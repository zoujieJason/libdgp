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

#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

std::vector<Eigen::Triplet<double>> cot_triplets;
std::vector<int> boundary_vertices;

Eigen::SparseMatrix<double> Grad; 

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd BC; 
Eigen::MatrixXd X;
Eigen::VectorXd phi;

Eigen::VectorXi iI;
Eigen::MatrixXd iV;
Eigen::MatrixXi iE; 

Eigen::RowVector3d red(igl::MAYA_RED(0), igl::MAYA_RED(1), igl::MAYA_RED(2));
Eigen::RowVector3d grey(igl::MAYA_GREY(0), igl::MAYA_GREY(1), igl::MAYA_GREY(2));
Eigen::RowVector3d green(igl::MAYA_GREEN(0), igl::MAYA_GREEN(1), igl::MAYA_GREEN(2));

Eigen::MatrixXd FC;
Eigen::MatrixXd XC;

Eigen::VectorXd level_sets;

double lambda = 0.1;
double edge_length = 0.2; 
double level_set = 0;
double line_width = 0.5;
int skip = 6;

Eigen::MatrixXd VG;
Eigen::MatrixXi FG; 
Eigen::VectorXi labels; 
int label = -1;
int kring = 5;

std::vector<std::vector<int>> VGF;

void update_phi()
{
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(cot_triplets.size() + boundary_vertices.size());
  triplets.insert(triplets.end(), cot_triplets.begin(), cot_triplets.end());
  std::vector<int> bc; 
  for(int i = 0; i < boundary_vertices.size(); i++)
  {
    if(i % skip == 0)
    {
      bc.push_back(boundary_vertices[i]);
    }
  }
  for(int i = 0; i < bc.size(); i++)
  {
    triplets.emplace_back(V.rows() + i, bc[i], lambda);
  }

  Eigen::SparseMatrix<double> A(V.rows() + bc.size(), V.rows());
  A.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::VectorXd X_vec = Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(X.data(), X.size(), 1);
  Eigen::VectorXd rhs(V.rows() + bc.size());
  rhs.setZero();
  rhs.head(V.rows()) = Grad.transpose() * X_vec;
  Eigen::SparseMatrix<double> At  = A.transpose();
  Eigen::SparseMatrix<double> AtA = At * A;
  Eigen::VectorXd             Atb = At * rhs;
  Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > solver(AtA);
  std::cout << "Solver info: " << solver.info() << std::endl;
  phi = solver.solve(Atb).eval();
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
  auto extract_faces = (K.array() == 1).count();
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

void extract_tooth()
{
  std::cout << "extracted tooth " << label << std::endl;
  Eigen::MatrixXd Vt; 
  Eigen::MatrixXi Ft;
  Eigen::VectorXi It, Jt, FIt;
  Eigen::VectorXi tooth_label = labels.unaryExpr([&](int l) { return static_cast<int>(l == label); });
  remove_unref(VG, FG, tooth_label, Vt, Ft, It, Jt, FIt);
  
  Eigen::VectorXi tooth_loop; 
  igl::boundary_loop(Ft, tooth_loop);
  tooth_loop = tooth_loop.unaryExpr([&Jt](int v) { return Jt(v); });

  std::vector<int> fuzzy_vertices; 
  std::vector<int> fuzzy_triangle;
  std::for_each(tooth_loop.data(), tooth_loop.data() + tooth_loop.size(), [&](int v) { 
    std::vector<int> kring_vertices, kring_faces;
    dgp::extract_kring(FG, VGF, v, kring, kring_vertices, kring_faces); 
    fuzzy_triangle.insert(fuzzy_triangle.end(), kring_faces.begin(), kring_faces.end());
  });
  std::for_each(fuzzy_triangle.begin(), fuzzy_triangle.end(), [&tooth_label](auto f) { tooth_label(f) = 1; });

  Eigen::VectorXi If, Jf, FIf;
  remove_unref(VG, FG, tooth_label, V, F, If, Jf, FIf);
  std::cout << "tooth " << label << " [#V" << V.rows() << " #F" << F.rows() << " ]." << std::endl;
  
  Eigen::VectorXi Fft;
  Fft.setZero(F.rows());
  FC.setZero(F.rows(), 3);
  for(int f = 0; f < F.rows(); ++f)
  {
    if(labels(FIf(f)) == label)
    {
      FC.row(f) = red;
      Fft(f) = 1;
    }
    else 
    {
      FC.row(f) = grey;
    }
  }

  remove_unref(V, F, Fft, Vt, Ft, It, Jt, FIt);
  igl::boundary_loop(Ft, tooth_loop);
  boundary_vertices.clear();
  boundary_vertices.reserve(tooth_loop.size());
  std::for_each(tooth_loop.data(), tooth_loop.data() + tooth_loop.size(), [&](int v) { boundary_vertices.push_back(Jt(v)); });
  std::cout << "boundary_vertices " << boundary_vertices.size() << std::endl;
}

void prepare()
{
  if(V.rows() == 0) return;
  if(F.rows() == 0) return;

  Eigen::MatrixXd C; 
  igl::cotmatrix_entries(V, F, C);
  Eigen::MatrixXi edges_patten; edges_patten.resize(3, 2);
  edges_patten << 1, 2, 2, 0, 0, 1; 
  cot_triplets.clear();
  cot_triplets.reserve(F.rows() * edges_patten.rows() * 4);
  for(int f = 0; f < F.rows(); ++f)
  {
    for(int e = 0; e < edges_patten.rows(); ++e)
    {
      int source = F(f, edges_patten(e, 0));
      int dest = F(f, edges_patten(e, 1));
      cot_triplets.emplace_back(source, dest, -C(f, e));
      cot_triplets.emplace_back(dest, source, -C(f, e)); 
      cot_triplets.emplace_back(source, source, C(f, e));
      cot_triplets.emplace_back(dest, dest, C(f, e));
    }
  }
  igl::grad(V, F, Grad);
  igl::barycenter(V, F, BC);

  std::vector<int> bad_vertices; 
  Eigen::MatrixXd PD1, PD2; 
  Eigen::VectorXd PV1, PV2; 
  igl::principal_curvature(V, F, PD1, PD2, PV1, PV2, bad_vertices);
  
  //Compute vector Laplacian and mass matrix
  Eigen::MatrixXi E, oE;//Compute Laplacian and mass matrix
  Eigen::SparseMatrix<double> vecL, vecM;
  igl::cr_vector_mass(V, F, E, vecM);
  igl::cr_vector_laplacian(V, F, E, oE, vecL);
  const int m = vecL.rows()/2; //The number of edges in the mesh

  //Convert the E / oE matrix format to list of edges / EMAP format required
  // by the functions constructing scalar Crouzeix-Raviart functions
  Eigen::MatrixXi Elist(m,2), EMAP(3*F.rows(),1);
  for(int i=0; i<F.rows(); ++i) {
    for(int j=0; j<3; ++j) {
      const int e = E(i,j);
      EMAP(i+j*F.rows()) = e;
      if(oE(i,j)>0) {
        Elist.row(e) << F(i, (j+1)%3), F(i, (j+2)%3);
      }
    }
  }
  Eigen::SparseMatrix<double> scalarL, scalarM;
  igl::crouzeix_raviart_massmatrix(V, F, Elist, EMAP, scalarM);
  igl::crouzeix_raviart_cotmatrix(V, F, Elist, EMAP, scalarL);

  //Compute edge midpoints & edge vectors
  Eigen::MatrixXd edgeMps, parVec, perpVec;
  igl::edge_midpoints(V, F, E, oE, edgeMps);
  igl::edge_vectors(V, F, E, oE, parVec, perpVec);

  std::vector<std::vector<int> > VF, VFi;
  igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);

  std::vector<double> paraVals(boundary_vertices.size(), 0.08);
  std::vector<double> perpVals(boundary_vertices.size(), 0.95);
  std::vector<int> boundary_edges(boundary_vertices.size(), -1);
  for(int i = 0; i < boundary_vertices.size(); i++)
  {
    int vi = boundary_vertices[i];
    int vj = boundary_vertices[(i + 1) % boundary_vertices.size()];

    for(auto f: VF[vi])
    {
      for(int j = 0; j < 3; ++j)
      {
        if(F(f, j) == vj)
        {
          for(int k = 0; k < 3; ++k)
          {
            if(F(f, k) != vi && F(f, k) != vj)
            {
              boundary_edges[i] = E(f, k);
              break;
            }
          }
        }
      }

      if(boundary_edges[i] != -1) break;
    }

    if(boundary_edges[i] == -1)
    {
      std::cerr << "Could not find edge for boundary vertex "<< std::endl;
      continue;
    }

    auto e = boundary_edges[i];
    if(PD2.row(vi).dot(PD2.row(vj)) < 0)
    {
      PD2.row(vj) *= -1.0;
    }
    Eigen::Vector3d pd = 0.5 * (PD2.row(vi) + PD2.row(vj)); pd.normalize();
    Eigen::Vector3d par = parVec.row(e); par.normalize();
    Eigen::Vector3d perp = perpVec.row(e); perp.normalize();
    Eigen::Vector3d n = par.cross(perp); n.normalize();
    pd = pd - pd.dot(n) * n;
    
    paraVals[i] = pd.dot(par);
    perpVals[i] = pd.dot(perp);
  }

  //Perform the vector heat method
  Eigen::VectorXi known = Eigen::VectorXi::Zero(2 * boundary_edges.size());
  Eigen::VectorXd knownVals = Eigen::VectorXd::Zero(2 * boundary_edges.size());
  Eigen::VectorXd Y0 = Eigen::VectorXd::Zero(2 * m), Yt;
  for(int i = 0; i < boundary_edges.size(); ++i)
  {
    known(i) = boundary_edges[i]; 
    known(i + boundary_edges.size()) = boundary_edges[i] + m;
    knownVals(i) = paraVals[i];
    knownVals(i + boundary_edges.size()) = perpVals[i];
    Y0(boundary_edges[i]) = paraVals[i];
    Y0(boundary_edges[i] + m) = perpVals[i];
  }

  const double t = 0.01;
  Eigen::SparseMatrix<double> Aeq;
  Eigen::VectorXd Beq;
  igl::min_quad_with_fixed
  (Eigen::SparseMatrix<double>(vecM+t*vecL), Eigen::VectorXd(-vecM*Y0), known, knownVals,
   Aeq, Beq, false, Yt);
  
  Eigen::VectorXd u0 = Eigen::VectorXd::Zero(m), ut;
  Eigen::VectorXi knownScal = Eigen::VectorXi::Zero(boundary_edges.size());
  Eigen::VectorXd knownScalVals = Eigen::VectorXd::Zero(boundary_edges.size());
  for(int i = 0; i < boundary_edges.size(); ++i)
  {
    u0(boundary_edges[i]) = sqrt(paraVals[i]*paraVals[i] + perpVals[i]*perpVals[i]);
    knownScal(i) = boundary_edges[i];
    knownScalVals(i) = u0(boundary_edges[i]);
  }
  
  igl::min_quad_with_fixed
  (Eigen::SparseMatrix<double>(scalarM+t*scalarL), Eigen::VectorXd(-scalarM*u0), knownScal,
   knownScalVals, Aeq, Beq, false, ut);

  Eigen::VectorXd phi0 = Eigen::VectorXd::Zero(m), phit;
  Eigen::VectorXd knownScalValsPhi = Eigen::VectorXd::Zero(boundary_edges.size());
  for(int i = 0; i < boundary_edges.size(); ++i)
  {
    phi0(boundary_edges[i]) = 1;
    knownScalValsPhi(i) = 1;
  }

  igl::min_quad_with_fixed
  (Eigen::SparseMatrix<double>(scalarM+t*scalarL), Eigen::VectorXd(-scalarM*phi0), knownScal,
   knownScalValsPhi, Aeq, Beq, false, phit);

  Eigen::ArrayXd Xtfactor = ut.array() /
  (phit.array() * (Yt.array().segment(0,m)*Yt.array().segment(0,m)
                   + Yt.array().segment(m,m)*Yt.array().segment(m,m)).sqrt());
  Eigen::VectorXd Xt(2*m);
  Xt.segment(0,m) = Xtfactor * Yt.segment(0,m).array();
  Xt.segment(m,m) = Xtfactor * Yt.segment(m,m).array();

  //Convert vector field for plotting
  Eigen::MatrixXd vecs(m, 3);
  for(int i=0; i<edgeMps.rows(); ++i) {
    vecs.row(i) = Xt(i)*parVec.row(i) + Xt(i+edgeMps.rows())*perpVec.row(i);
  }

  XC.setZero(F.rows(), 3);
  X.setZero(F.rows(), 3);
  for(int i = 0; i < F.rows(); ++i)
  {
    XC.row(i) << 0.1, 0.1, 0.1;
    for(int j = 0; j < F.cols(); ++j)
    {
      X.row(i) += vecs.row(E(i, j));   
    }
    X.row(i) /= static_cast<double>(F.cols());
  }
  X.rowwise().normalize();
}

int main(int argc, char ** argv)
{
  const double howMuchToSmoothBy = 1e-1;
  const int howManySmoothingInterations = 50;

  std::string meshname(argv[1]);
  if(!igl::read_triangle_mesh(meshname,VG,FG)) 
  {
    std::cout << "Failed to load mesh " << meshname << std::endl;
    return 0;
  }
  std::vector<std::vector<int>> VGFi;
  igl::vertex_triangle_adjacency(VG, FG, VGF, VGFi);

  std::string labelname(argv[2]);
  std::vector<int> veclabel;
  dgp::read_labels(labelname, veclabel);
  if(veclabel.size() != FG.rows()) 
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
  viewer.data().set_mesh(VG,FG);
  viewer.data().show_lines = false;

  menu.callback_draw_viewer_menu = [&]()
  {
    menu.draw_viewer_menu();
    if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
    {
      ImGui::InputDouble("lambda", &lambda, 0, 0, "%.6f");
      ImGui::InputDouble("edge_length", &edge_length, 0, 0, "%.4f");
      ImGui::InputDouble("level-set", &level_set, 0, 0, "%.4f");
      ImGui::InputDouble("line_width", &line_width, 0, 0, "%.4f");
      ImGui::InputInt("skip", &skip, 0, 0);
      ImGui::InputInt("kring", &kring, 0, 0);
      static bool is_initialized = false;
      static int tooth_number = 1;
      static std::vector<std::string> tooth_numbers;

      if (!is_initialized) {
          std::set<int> label_set(veclabel.begin(), veclabel.end());
          std::for_each(label_set.begin(), label_set.end(), [&](int l) {
              tooth_numbers.push_back(std::to_string(l));
          });
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
          prepare();
          viewer.data().clear();
          viewer.data().set_mesh(V,F);
          viewer.core().align_camera_center(V, F);
          viewer.data().set_colors(FC);
          return true;
        case '1':
          viewer.data().clear();
          viewer.data().set_mesh(V,F);
          update_phi();
          level_sets << level_set;
          igl::isolines(V, F, phi, level_sets, iV, iE, iI);
          viewer.data().set_colors(FC);
          viewer.data().set_edges(iV, iE, green);
          return true;

        case '2':
          viewer.data().clear();
          viewer.data().set_mesh(V,F);
          update_phi();
          level_sets << level_set;
          igl::isolines(V, F, phi, level_sets, iV, iE, iI);
          viewer.data().set_data(phi);
          viewer.data().add_edges(BC, BC + edge_length*X, XC);
          return true;

        default:
          return false;
      }
  };
  viewer.launch();

  return 0;
}
