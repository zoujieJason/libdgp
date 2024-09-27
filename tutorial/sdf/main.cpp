#include <igl/read_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/remove_unreferenced.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/orient_halfedges.h>
#include <igl/cr_vector_laplacian.h>
#include <igl/cr_vector_mass.h>
#include <igl/crouzeix_raviart_cotmatrix.h>
#include <igl/crouzeix_raviart_massmatrix.h>
#include <igl/edge_midpoints.h>
#include <igl/edge_vectors.h>
#include <igl/average_from_edges_onto_vertices.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/heat_geodesics.h>
#include <igl/principal_curvature.h>
#include <igl/isolines.h>
#include <igl/remesh_along_isoline.h>
#include <igl/slice.h>
#include <igl/barycentric_coordinates.h>
#include <igl/edge_lengths.h>

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/Geometry>

#include <iostream>
#include <set>
#include <limits>
#include <stdlib.h>

#include "dgp/read_labels.h"
#include "dgp/extract_kring.h"

#include <igl/remove_unreferenced.h>
#include <igl/boundary_loop.h>
#include <igl/grad.h>
#include <igl/intrinsic_delaunay_cotmatrix.h>
#include <igl/massmatrix_intrinsic.h>
#include <igl/invert_diag.h>
#include <igl/pathinfo.h>

#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

Eigen::RowVector3d edge_basis(int i)
{
  // f = (i,j,k);
  // index 
  //   0 -> e(j,k)
  //   1 -> e(k,i); 
  //   2 -> e(i,j);
  
  int j = (i + 1) % 3;
  int k = (i + 2) % 3;

  Eigen::RowVector3d e1;
  e1(i) = 0; e1(j) = -1; e1(k) = 1;
  return e1;
}

void crouzeix_raviart_connection_laplacian(
  int nE,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXi & E,
  const Eigen::MatrixXi & oE,
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXd & l_sq,
  const Eigen::MatrixXd & l,
  const Eigen::VectorXd & dblA,
  Eigen::SparseMatrix<std::complex<double>> & L)
{
  L = Eigen::SparseMatrix<std::complex<double>>(nE, nE);
  std::vector<Eigen::Triplet<std::complex<double>>> tripletList;
  for(int f = 0; f < F.rows(); ++f)
  {
    for(int i = 0; i < 3; ++i)
    {
      int j = (i + 1) % 3;
      int k = (i + 2) % 3;
      double weight = 4. * C(f, i);
      double sign_jk = static_cast<double>(oE(f, j) * oE(f, k));
      double cosTheta = (l_sq(f, j) + l_sq(f, k) - l_sq(f, i)) / (2. * l(f, j) * l(f, k));
      double sinTheta = dblA(f) / (l(f, j) * l(f, k));
      std::complex<double> rot_jk(-cosTheta, sinTheta);

      int ej = E(f, j);
      int ek = E(f, k);
      tripletList.emplace_back(ej, ej, weight);
      tripletList.emplace_back(ek, ek, weight);
      tripletList.emplace_back(ej, ek, -weight * rot_jk * sign_jk);
      tripletList.emplace_back(ek, ej, -weight * std::conj(rot_jk) * sign_jk);
    }
  }
  L.setFromTriplets(tripletList.begin(), tripletList.end());
}

void crouzeix_raviart_connection_mass(
  const std::vector<std::vector<int>> EF,
  const Eigen::VectorXd & dblA,
  Eigen::SparseMatrix<double> & M)
{
  size_t m = EF.size();
  M = Eigen::SparseMatrix<double>(m, m);
  std::vector<Eigen::Triplet<double>> tripletList;
  for(int i = 0; i < m; ++i)
  {
    for(int f: EF[i])
    {
      tripletList.emplace_back(i, i, dblA(f) / 6.);
    }
  }
  M.setFromTriplets(tripletList.begin(), tripletList.end());
}

void faces_vector_feild(
  int nf, 
  const Eigen::MatrixXi & E, 
  const Eigen::MatrixXi & oE,
  const Eigen::MatrixXd & l_sq, 
  const Eigen::MatrixXd & l,
  const Eigen::Matrix<std::complex<double>, -1, 1> &X,
  Eigen::MatrixXd &Xd)
{
  Xd = Eigen::MatrixXd::Zero(nf, 3);
  for(int f = 0; f < nf; ++f)
  {
    double l_jk_sq = l_sq(f, 0);
    double l_ki_sq = l_sq(f, 1);
    double l_ij_sq = l_sq(f, 2);
    double ai = l_ij_sq + l_ki_sq - l_jk_sq;
    double aj = l_jk_sq + l_ij_sq - l_ki_sq;
    double ak = l_ki_sq + l_jk_sq - l_ij_sq;
    Eigen::RowVector3d bc = {ai, aj, ak};
    double s = 0.5 * (l(f, 0) + l(f, 1) + l(f, 2));
    double A = std::sqrt(s * (s - l(f, 0)) * (s - l(f, 1)) * (s - l(f, 2)));
    if(A < 0.)
    {
      std::cout << "warning : A < 0. " << std::endl;
    }

    Eigen::RowVector3d fv = Eigen::RowVector3d::Zero();
    for(int i = 0; i < 3; ++i)
    {
      Eigen::RowVector3d e1 = edge_basis(i);
      e1 *= static_cast<double>(oE(f, i));

      Eigen::RowVector3d e2{e1(2) * bc(2) - e1(1) * bc(1),
                            e1(0) * bc(0) - e1(2) * bc(2),
                            e1(1) * bc(1) - e1(0) * bc(0)};
      e2 /= (4. * A);
      double norm_e1 = std::sqrt(-((l_ij_sq * e1(0) * e1(1) + 
                                    l_jk_sq * e1(1) * e1(2) +
                                    l_ki_sq * e1(2) * e1(0))));
      e1 /= norm_e1;

      double norm_e2 = std::sqrt(-((l_ij_sq * e2(0) * e2(1) + 
                                    l_jk_sq * e2(1) * e2(2) +
                                    l_ki_sq * e2(2) * e2(0))));
      e2 /= norm_e2;
      fv += std::real(X(E(f, i))) * e1;
      fv += std::imag(X(E(f, i))) * e2;
    }
    double norm_fv = std::sqrt(-((l_ij_sq * fv(0) * fv(1) + 
                                  l_jk_sq * fv(1) * fv(2) +
                                  l_ki_sq * fv(2) * fv(0))));
    if(norm_fv == 0)
    {
      Xd.row(f) = Eigen::RowVector3d::Zero();
    }
    else
    {
      Xd.row(f) = fv / norm_fv;
    }
  }
}

void remove_unref(
  const Eigen::MatrixXd & Vin,
  const Eigen::MatrixXi & Fin,
  const Eigen::VectorXi &K, 
  Eigen::MatrixXd & Vout,
  Eigen::MatrixXi & Fout,
  Eigen::VectorXi & Iout,
  Eigen::VectorXi & Jout,
  Eigen::VectorXi & FIout)
{
  Eigen::Index extract_faces = (K.array() == 1).count();
  FIout.setZero(extract_faces);
  Eigen::MatrixXi Fi(extract_faces, 3);
  for(int i = 0, f = 0; i < Fin.rows(); i++)
  {
    if(K(i) == 1)
    {
      FIout(f) = i;
      Fi.row(f++) = Fin.row(i);
    }
  }

  igl::remove_unreferenced(Vin, Fi, Vout, Fout, Iout, Jout);
}

int label = -1; 
Eigen::VectorXi labels;
double t = 1.;
Eigen::VectorXd Phi;
double lambda = 0.1;
double line_length = 0.1;
int skip = 6;

Eigen::RowVector3d red(igl::MAYA_RED(0), igl::MAYA_RED(1), igl::MAYA_RED(2));
Eigen::RowVector3d grey(igl::MAYA_GREY(0), igl::MAYA_GREY(1), igl::MAYA_GREY(2));
Eigen::RowVector3d green(igl::MAYA_GREEN(0), igl::MAYA_GREEN(1), igl::MAYA_GREEN(2));

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd centers;
Eigen::MatrixXd Xd;
Eigen::MatrixXd Xv; 
using VectorXComplex = Eigen::Matrix<std::complex<double>, -1, 1>;
std::vector<std::vector<int>> boundaries;

Eigen::VectorXi iI;
Eigen::MatrixXd iV;
Eigen::MatrixXi iE;
Eigen::VectorXd level_sets;
double level_set = 0.;
bool LTL = false;
bool hard_constriant = false;

void extract_tooth()
{
  std::cout << "extracted tooth " << label << std::endl;
  Eigen::MatrixXd Vt; 
  Eigen::MatrixXi Ft;
  Eigen::VectorXi It, Jt, FIt;
  remove_unref(V, F, labels.unaryExpr([&](int l) { return static_cast<int>(l == label); }).eval(), Vt, Ft, It, Jt, FIt);
  std::cout << "tooth " << label << " [#V" << V.rows() << " #F" << F.rows() << " ]." << std::endl;

  // this method might not get all boundary vertices 
  boundaries.clear();
  igl::boundary_loop(Ft, boundaries);
  for(auto &boundary: boundaries)
  {
    std::for_each(boundary.begin(), boundary.end(), [&Jt](int &v) { v = Jt(v); });
  }

  Eigen::MatrixXd C, l, l_sq;; 
  Eigen::VectorXd dblA;
  igl::doublearea(V, F, dblA);
  igl::squared_edge_lengths(V, F, l_sq);
  l = l_sq.array().sqrt();
  igl::cotmatrix_entries(l, C);

  Eigen::SparseMatrix<double> Grad; 
  igl::grad(V, F, Grad);

  igl::barycenter(V, F, centers);

  Eigen::MatrixXi E, oE;
  igl::orient_halfedges(F, E, oE);
  int nE = E.maxCoeff() + 1;

  Eigen::MatrixXi Elist, EMAP;
  std::vector<std::vector<int>> EF(nE); 
  Elist.setZero(nE,2);
  EMAP.setZero(3*F.rows(),1);
  for(int i=0; i<F.rows(); ++i) {
    for(int j=0; j<3; ++j) {
      const int e = E(i,j);
      EMAP(i+j*F.rows()) = e;
      EF[e].push_back(i);
      if(oE(i,j)>0) {
        Elist.row(e) << F(i, (j+1)%3), F(i, (j+2)%3);
      }
    }
  }

  Eigen::SparseMatrix<std::complex<double>> Lconn;
  Eigen::SparseMatrix<double> Mconn;
  crouzeix_raviart_connection_laplacian(nE, F, E, oE, C, l_sq, l, dblA, Lconn);
  crouzeix_raviart_connection_mass(EF, dblA, Mconn);

  std::vector<std::vector<int> > VF, VFi;
  igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);

  VectorXComplex X0 = VectorXComplex::Zero(nE);
  for(auto &boundary: boundaries)
  {
    for(int i = 0; i < boundary.size() - 1; i++)
    {
      int vi = boundary[i];
      int vj = boundary[(i + 1) % boundary.size()];

      int eid = -1;
      for(int f: VF[vi])
      {
        for(int j = 0; j < 3; ++j)
        {
          if(F(f, j) == vj)
          {
            for(int k = 0; k < 3; ++k)
            {
              if(F(f, k) != vi && F(f, k) != vj)
              {
                eid = E(f, k);
                std::complex<double> innerProd(0, Elist(eid, 0) == vi ? 1 : -1);
                std::complex<double> contrib = l(f, k) * innerProd;
                X0[eid] += contrib;
                break;
              }
            }
          }
        }

        if(eid != -1) break;
      }
    }
  }

  Eigen::SparseMatrix<std::complex<double>> vectorOp = Mconn.cast<std::complex<double>>() + t * Lconn;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> complex_solver;
  vectorOp.makeCompressed();
  complex_solver.compute(vectorOp);
  std::cout << "complex_solver factorization : " << complex_solver.info() << std::endl;
  VectorXComplex Xt = complex_solver.solve(X0);
  std::cout << "complex_solver solve : " << complex_solver.info() << std::endl;
  faces_vector_feild(F.rows(), E, oE, l_sq, l, Xt, Xd);

  Eigen::VectorXd Div = Eigen::VectorXd::Zero(V.rows());
  for(int v = 0; v < V.rows(); ++v)
  {
    for(int f: VF[v])
    {
      double l_jk_sq = l_sq(f, 0);
      double l_ki_sq = l_sq(f, 1);
      double l_ij_sq = l_sq(f, 2);
      for(int i = 0; i < 3; ++i)
      {
        if(F(f, i) == v)
        {
          Eigen::RowVector3d x = Xd.row(f);
          {
            // k = (i + 2) % 3 -> e(i, j);
            int vi = (i + 2) % 3;
            Eigen::RowVector3d eA = edge_basis(vi);
            double term1 = l_ij_sq * (eA(1) * x(0) + eA(0) * x(1));
            double term2 = l_jk_sq * (eA(1) * x(2) + eA(2) * x(1));
            double term3 = l_ki_sq * (eA(2) * x(0) + eA(0) * x(2));
            double dot = -0.5 * (term1 + term2 + term3);
            Div(v) += C(f, vi) * dot;
          }

          {
            // j = (i + 1) % 3 -> e(k, i);
            int vi = (i + 1) % 3;
            Eigen::RowVector3d eB = -edge_basis(vi); // careful with auto Eigen type
            double term1 = l_ij_sq * (eB(1) * x(0) + eB(0) * x(1));
            double term2 = l_jk_sq * (eB(1) * x(2) + eB(2) * x(1));
            double term3 = l_ki_sq * (eB(2) * x(0) + eB(0) * x(2));
            double dot = -0.5 * (term1 + term2 + term3);
            Div(v) += C(f, vi) * dot;
            }
          break;
        }
      }
    }
  }

  Xv.setZero(V.rows(), 3);
  for(int v = 0; v < V.rows(); ++v)
  {
    for(auto f: VF[v])
    {
      Xv.row(v) += Xd.row(f);
    }
    Xv.row(v) /= static_cast<double>(VF[v].size());
  }

  Eigen::SparseMatrix<double> A;
  igl::cotmatrix(V, F, A);
  Eigen::SparseMatrix<double> Lap = -A;
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > poisson_solver;
  if(LTL)
  {
    if(hard_constriant)
    {
      Eigen::SparseMatrix<double> Q = Lap.transpose() * Lap;
      Eigen::MatrixXd RHS = 2.0 * Lap * Div;
      int n_constraints = std::accumulate(boundaries.begin(), boundaries.end(), 0, [](int sum, const std::vector<int> &boundary) { return sum + boundary.size(); });
      std::cout << "n_constraints: " << n_constraints << std::endl;
      Eigen::VectorXi known = Eigen::VectorXi::Zero(n_constraints); 
      Eigen::VectorXd known_vals = Eigen::VectorXd::Zero(n_constraints);
      int count = 0; 
      for(const auto &boundary: boundaries)
      {
        for(int i = 0; i < boundary.size(); i++)
        {
          known(count) = boundary[i];
          known_vals(count++) = RHS(boundary[i]);
        }
      }
      Eigen::SparseMatrix<double> Aeq;
      Eigen::VectorXd Beq;
      igl::min_quad_with_fixed(Q, -RHS, known, known_vals, Aeq, Beq, false, Phi);
    }
    else 
    {
      Lap.makeCompressed();
      poisson_solver.compute(Lap.transpose() * Lap);
      std::cout << "poisson_solver compute info: " << poisson_solver.info() << std::endl;
      Phi = poisson_solver.solve(2.0 * Lap.transpose() * Div);
    }
  }
  else 
  {
    if(hard_constriant)
    {
      int n_constraints = std::accumulate(boundaries.begin(), boundaries.end(), 0, [](int sum, const std::vector<int> &boundary) { return sum + boundary.size(); });
      std::cout << "n_constraints: " << n_constraints << std::endl;
      Eigen::VectorXi known = Eigen::VectorXi::Zero(n_constraints); 
      Eigen::VectorXd known_vals = Eigen::VectorXd::Zero(n_constraints);
      int count = 0; 
      for(const auto &boundary: boundaries)
      {
        for(int i = 0; i < boundary.size(); i++)
        {
          known(count) = boundary[i];
          known_vals(count++) = Div(boundary[i]);
        }
      }
      Eigen::SparseMatrix<double> Aeq;
      Eigen::VectorXd Beq;
      std::cout << "igl::min_quad_with_fixed solve info: " <<
        igl::min_quad_with_fixed(Lap, Div, known, known_vals, Aeq, Beq, true, Phi) << 
      std::endl;
    }
    else 
    {
      Lap.makeCompressed();
      poisson_solver.compute(Lap);
      std::cout << "poisson_solver compute info: " << poisson_solver.info() << std::endl;
      Phi = poisson_solver.solve(Div);
    }
  }
  std::cout << "poisson_solver solve info: " << poisson_solver.info() << std::endl; 
  Phi *= -1.;
  // Compute the average value of phi along the input (integrate phi along the input geometry.)
  double shift = 0.;
  double normalization = 0.;
  for(const auto &boundary: boundaries)
  {
    for(int i = 0; i < boundary.size(); i++)
    {
      int vi = boundary[i];
      int vj = boundary[(i + 1) % boundary.size()];

      int eid = -1;
      for(int f: VF[vi])
      {
        int i_ = -1; 
        for(int j = 0; j < 3; ++j)
        {
          if(F(f, j) == vi)
          {
            i_ = j;
            break;
          }
        }

        int j_ = -1;
        for(int j = 0; j < 3; ++j)
        {
          if(F(f, j) == vj)
          {
            j_ = j;
            break;
          }
        }

        int k_ = -1; 
        for(int j = 0; j < 3; ++j)
        {
          if(F(f, j) != vi && F(f, j) != vj)
          {
            k_ = j;
            break;
          }
        }

        if(i_ == -1 || j_ == -1 || k_ == -1)
        {
          continue;
        }

        Eigen::RowVector3d vi_coords;
        vi_coords(i_) = 1; vi_coords(j_) = 0; vi_coords(k_) = 0;
        Eigen::RowVector3d vj_coords;
        vj_coords(i_) = 0; vj_coords(j_) = 1; vj_coords(k_) = 0;
        double l_ij = l(f, k_);
        Eigen::RowVector3d midp = 0.5 * (vi_coords + vj_coords);
        Eigen::RowVector3d phi; 
        phi(i_) = Phi(F(f, i_)); phi(j_) = Phi(F(f, j_)); phi(k_) = Phi(F(f, k_));
        shift += l_ij * midp.dot(phi);
        normalization += l_ij;
        break;
      }
    }
  }
  shift /= normalization;
  std::cout << "shjft: " << shift << std::endl;
  Phi -= (Eigen::VectorXd::Ones(V.rows()) * shift);
}

int main(int argc, char ** argv)
{
  std::string meshname(argv[1]);
  if(!igl::read_triangle_mesh(meshname,V,F)) 
  {
    std::cout << "Failed to load mesh " << meshname << std::endl;
    return 0;
  }

  std::string labelname(argv[2]);
  std::vector<int> veclabel;
  dgp::read_labels(labelname, veclabel);
  if(veclabel.size() != F.rows()) 
  {
    std::cout << "Failed to load label " << labelname << std::endl;
    return 0;
  }
  labels = Eigen::Map<Eigen::VectorXi>(veclabel.data(), veclabel.size());
  level_sets.resize(1);

  igl::opengl::glfw::Viewer viewer;
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  viewer.plugins.push_back(&plugin);
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  plugin.widgets.push_back(&menu);
  viewer.data().set_mesh(V,F);
  viewer.data().show_lines = false;
  Eigen::MatrixXd colors; 
  colors.setZero(F.rows(), 3);
  for(int f = 0; f < F.rows(); ++f)
  {
    colors.row(f) = labels(f) == 1 ? red : grey;
  }
  viewer.data().set_colors(colors);

  menu.callback_draw_viewer_menu = [&]()
  {
    menu.draw_viewer_menu();
    if (ImGui::CollapsingHeader("Parameters", ImGuiTreeNodeFlags_DefaultOpen))
    {
      static bool LTL_mode = false;
      if (ImGui::Checkbox("LTL mode", &LTL_mode))
      {
        LTL = LTL_mode;
      }

      static bool hard_constriant_mode = false;
      if (ImGui::Checkbox("hard constriant mode", &hard_constriant_mode))
      {
        hard_constriant = hard_constriant_mode;
      }
      
      ImGui::InputDouble("lambda", &lambda, 0, 0, "%.6f");
      ImGui::InputDouble("t", &t, 0, 0, "%.6f");
      ImGui::InputDouble("line length", &line_length, 0, 0, "%.4f");
      ImGui::InputDouble("level-set", &level_set, 0, 0, "%.4f");
      ImGui::InputInt("skip", &skip, 0, 0);

      static bool is_initialized = false;
      static int tooth_number = 1;
      static std::vector<std::string> tooth_numbers;

      if (!is_initialized) {
          std::set<int> label_set(veclabel.begin(), veclabel.end());
          std::for_each(label_set.begin(), label_set.end(), [&](int l) { tooth_numbers.push_back(std::to_string(l)); });
          is_initialized = true;
      }

      ImGui::Combo("tooth number", &tooth_number, tooth_numbers);
      label = std::stoi(tooth_numbers[tooth_number]);
    }
  };

  viewer.callback_key_pressed = [&](
    igl::opengl::glfw::Viewer& viewer,
    unsigned char key,
    int modifier) -> bool
  {
      switch (key)
      {
        case 'e':
          extract_tooth();
          viewer.data().clear();
          viewer.data().set_mesh(V,F);
          viewer.data().set_data(Phi);
          viewer.data().add_edges(centers, centers + line_length * Xd.rowwise().normalized(), red);
          return true;
        case '1':
          viewer.data().clear();
          viewer.data().set_mesh(V,F);
          viewer.data().show_lines = true;
          level_sets << level_set;
          igl::isolines(V, F, Phi, level_sets, iV, iE, iI);
          viewer.data().set_edges(iV, iE, green);
          viewer.data().set_colors(colors);
        return true;
        default:
          return false;
      }
  };
  viewer.launch();
  return 0;
}
