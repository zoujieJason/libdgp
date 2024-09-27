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

int count_bits(int n, int bits = sizeof(int) * 8)
{
  int count = 0;
  for(int i = 0; i < bits; i++)
  {
    count += n & 1;
    n >>= 1;
  }
  return count;
}

void remove_bit(int &type, int label)
{
  type &= ~(1 << label);
}

std::vector<Eigen::Triplet<double>> cot_triplets;
std::vector<int> boundary_vertices;

Eigen::SparseMatrix<double> Grad; 
Eigen::SparseMatrix<double> Div; 

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd BC; 
Eigen::MatrixXd X;
Eigen::VectorXd phi;
Eigen::VectorXd dblA;
Eigen::VectorXd VPV1;
Eigen::VectorXd VPV2;
Eigen::MatrixXd VPD1;
Eigen::MatrixXd VPD2;
Eigen::MatrixXd VPN;
Eigen::VectorXi FIt;

Eigen::VectorXi iI;
Eigen::MatrixXd iV;
Eigen::MatrixXi iE; 

Eigen::RowVector3d red(igl::MAYA_RED(0), igl::MAYA_RED(1), igl::MAYA_RED(2));
Eigen::RowVector3d grey(igl::MAYA_GREY(0), igl::MAYA_GREY(1), igl::MAYA_GREY(2));
Eigen::RowVector3d green(igl::MAYA_GREEN(0), igl::MAYA_GREEN(1), igl::MAYA_GREEN(2));

Eigen::MatrixXd FC;
Eigen::MatrixXd XC;

Eigen::VectorXd level_sets;

double t = 0.5;
double lambda = 0.5;
double edge_length = 0.2; 
double level_set = 0;
double line_width = 0.5;
int skip = 6;

Eigen::MatrixXd VG;
Eigen::MatrixXi FG; 
Eigen::VectorXi labels; 
Eigen::VectorXi vertices_type; 
int n_labels = 0;
int label = -1;
int kring = 6;
bool accumlative = true; 
bool intrinsic = false; 
bool edge_based_Laplace = false; 
bool preserve_source_normals = false;
bool min_quad_grad = true;

std::vector<std::vector<int>> VGF;

std::string dir;
std::string basename;
std::string ext;
std::string filename;

Eigen::MatrixXi E, oE;
Eigen::SparseMatrix<double> vecL, vecM;
Eigen::MatrixXi Elist, EMAP; 
Eigen::SparseMatrix<double> scalarL, scalarM;
Eigen::MatrixXd edgeMps, parVec, perpVec;
std::vector<double> paraVals, perpVals;
std::vector<int> boundary_edges;

Eigen::SparseMatrix<std::complex<double>> Lconn;
Eigen::SparseMatrix<double> Mconn;
Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> X0, Xt;
Eigen::MatrixXd Xd; 

double howMuchToSmoothBy = 0.1;
int howManySmoothingInterations = 50;

Eigen::MatrixXd l_sq;
Eigen::MatrixXd l;

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
  Eigen::VectorXi It, Jt;
  Eigen::VectorXi tooth_label = labels.unaryExpr([&](int l) { return static_cast<int>(l == label); });
  remove_unref(VG, FG, tooth_label, Vt, Ft, It, Jt, FIt);
  
  // this method might not get all boundary vertices 
  // Eigen::VectorXi tooth_loop; 
  // igl::boundary_loop(Ft, tooth_loop);
  // tooth_loop = tooth_loop.unaryExpr([&Jt](int v) { return Jt(v); });
  std::vector<std::vector<int>> loops;
  igl::boundary_loop(Ft, loops);
  int curr_label_type = 1 << label;
  std::vector<int> loop; 
  for(int i = 0; i < Vt.rows(); ++i)
  {
    auto type = vertices_type(Jt(i));
    if(count_bits(type, n_labels) == 1)
    {
      continue;
    }
    assert(type & curr_label_type);
    loop.push_back(Jt(i));
  }
  Eigen::VectorXi tooth_loop = Eigen::Map<Eigen::VectorXi>(loop.data(), loop.size());

  std::vector<int> fuzzy_vertices; 
  std::vector<int> fuzzy_triangle;
  std::for_each(tooth_loop.data(), tooth_loop.data() + tooth_loop.size(), [&](int v) { 
    std::vector<int> kring_vertices, kring_faces;
    dgp::extract_kring(FG, VGF, v, kring, kring_vertices, kring_faces); 
    fuzzy_triangle.insert(fuzzy_triangle.end(), kring_faces.begin(), kring_faces.end());
  });
  std::for_each(fuzzy_triangle.begin(), fuzzy_triangle.end(), [&tooth_label](auto f) { tooth_label(f) = 1; });

  Eigen::VectorXi If, Jf, FIf;
  remove_unref(VG, FG, tooth_label, V, F, If, Jf, FIt);
  std::cout << "tooth " << label << " [#V" << V.rows() << " #F" << F.rows() << " ]." << std::endl;
  
  Eigen::VectorXi Fft;
  Fft.setZero(F.rows());
  FC.setZero(F.rows(), 3);
  for(int f = 0; f < F.rows(); ++f)
  {
    if(labels(FIt(f)) == label)
    {
      FC.row(f) = red;
      Fft(f) = 1;
    }
    else 
    {
      FC.row(f) = grey;
    }
  }
   
  Eigen::VectorXi FIt_; 
  remove_unref(V, F, Fft, Vt, Ft, It, Jt, FIt_);
  igl::boundary_loop(Ft, tooth_loop);
  boundary_vertices.clear();
  boundary_vertices.reserve(tooth_loop.size());
  std::for_each(tooth_loop.data(), tooth_loop.data() + tooth_loop.size(), [&](int v) { boundary_vertices.push_back(Jt(v)); });
  std::cout << "boundary_vertices " << boundary_vertices.size() << std::endl;

  // prepare 
  igl::doublearea(V, F, dblA);
  igl::squared_edge_lengths(V, F, l_sq);
  l = l_sq.array().sqrt();

  Eigen::MatrixXd C; 
  igl::cotmatrix_entries(l, C);
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

  Eigen::MatrixXd VA, VB, VC; 
  igl::slice(V, F.col(0), Eigen::VectorXi::LinSpaced(3, 0, 3), VA);
  igl::slice(V, F.col(0), Eigen::VectorXi::LinSpaced(3, 0, 3), VB);
  igl::slice(V, F.col(0), Eigen::VectorXi::LinSpaced(3, 0, 3), VC);
  Eigen::MatrixXd FaceC; 
  igl::barycenter(V, F, FaceC);
  Eigen::MatrixXd BaryC;
  igl::barycentric_coordinates(FaceC, VA, VB, VC, BaryC);

  std::vector<int> bad_vertices; 
  Eigen::MatrixXd PD1, PD2; 
  Eigen::VectorXd PV1, PV2; 
  igl::principal_curvature(V, F, PD1, PD2, PV1, PV2, bad_vertices);
  VPD1 = PD1;
  VPD2 = PD2;
  VPV1 = PV1;
  VPV2 = PV2;
  VPN = PD2;
  for(int i = 0; i < PD1.rows(); ++i)
  {
    VPN.row(i) = Eigen::RowVector3d(VPD1.row(i)).cross(Eigen::RowVector3d(VPD2.row(i)));
  }
  std::cout << VPV1.minCoeff() << " " << VPV1.maxCoeff() << std::endl;  
  std::cout << VPV2.minCoeff() << " " << VPV2.maxCoeff() << std::endl;

  Eigen::MatrixXi E_, oE_;
  igl::cr_vector_mass(V, F, E_, vecM);
  igl::cr_vector_laplacian(V, F, E_, oE_, vecL);
  E = E_;
  oE = oE_;
  const int m = vecL.rows()/2;

  std::vector<std::vector<int>> EF(m); 
  Elist.setZero(m,2);
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

  igl::crouzeix_raviart_massmatrix(V, F, Elist, EMAP, scalarM);
  igl::crouzeix_raviart_cotmatrix(V, F, Elist, EMAP, scalarL);

  igl::edge_midpoints(V, F, E, oE, edgeMps);
  igl::edge_vectors(V, F, E, oE, parVec, perpVec);

  std::vector<std::vector<int> > VF, VFi;
  igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);

  paraVals = std::vector<double>(boundary_vertices.size(), 0.08);
  perpVals = std::vector<double>(boundary_vertices.size(), 0.95);
  boundary_edges = std::vector<int>(boundary_vertices.size(), -1);
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
    // auto e = boundary_edges[i];
    // if(PD2.row(vi).dot(PD2.row(vj)) < 0)
    // {
    //   PD2.row(vj) *= -1.0;
    // }
    // Eigen::Vector3d pd = 0.5 * (PD2.row(vi) + PD2.row(vj)); pd.normalize();
    // Eigen::Vector3d par = parVec.row(e); par.normalize();
    // Eigen::Vector3d perp = perpVec.row(e); perp.normalize();
    // Eigen::Vector3d n = par.cross(perp); n.normalize();
    // pd = pd - pd.dot(n) * n;
    
    // paraVals[i] = pd.dot(par);
    // perpVals[i] = pd.dot(perp);
  }
}

void calculateX()
{
  const int m = vecL.rows()/2;
  Eigen::VectorXi known;
  Eigen::VectorXd knownVals;
  if(preserve_source_normals)
  {
    known = Eigen::VectorXi::Zero(2 * boundary_edges.size());
    knownVals = Eigen::VectorXd::Zero(2 * boundary_edges.size());
  }

  Eigen::VectorXd Y0 = Eigen::VectorXd::Zero(2 * m), Yt;
  for(int i = 0; i < boundary_edges.size(); ++i)
  {
    if(preserve_source_normals) known(i) = boundary_edges[i]; 
    if(preserve_source_normals) known(i + boundary_edges.size()) = boundary_edges[i] + m;
    if(preserve_source_normals) knownVals(i) = paraVals[i];
    if(preserve_source_normals) knownVals(i + boundary_edges.size()) = perpVals[i];
    Y0(boundary_edges[i]) = paraVals[i];
    Y0(boundary_edges[i] + m) = perpVals[i];
  }

  Eigen::SparseMatrix<double> Aeq;
  Eigen::VectorXd Beq;
  for(int i=0; i<howManySmoothingInterations; ++i) 
  {
    igl::min_quad_with_fixed
    (Eigen::SparseMatrix<double>(vecM+t*vecL), Eigen::VectorXd(-vecM*Y0), known, knownVals,
    Aeq, Beq, false, Yt);
    Y0 = Yt;
  }

  Eigen::VectorXi knownScal;
  Eigen::VectorXd knownScalVals;
  if(preserve_source_normals)
  {
    knownScal = Eigen::VectorXi::Zero(boundary_edges.size());
    knownScalVals = Eigen::VectorXd::Zero(boundary_edges.size());
  }

  Eigen::VectorXd u0 = Eigen::VectorXd::Zero(m), ut;
  for(int i = 0; i < boundary_edges.size(); ++i)
  {
    u0(boundary_edges[i]) = sqrt(paraVals[i]*paraVals[i] + perpVals[i]*perpVals[i]);
    if(preserve_source_normals) knownScal(i) = boundary_edges[i];
    if(preserve_source_normals) knownScalVals(i) = u0(boundary_edges[i]);
  }

  for(int i=0; i<howManySmoothingInterations; ++i) 
  {
    igl::min_quad_with_fixed
    (Eigen::SparseMatrix<double>(scalarM+t*scalarL), Eigen::VectorXd(-scalarM*u0), knownScal,
    knownScalVals, Aeq, Beq, false, ut);
    u0 = ut;
  }

  Eigen::VectorXd knownScalValsPhi;
  if(preserve_source_normals)
  {
    knownScalValsPhi = Eigen::VectorXd::Zero(boundary_edges.size());
  }

  Eigen::VectorXd phi0 = Eigen::VectorXd::Zero(m), phit;
  for(int i = 0; i < boundary_edges.size(); ++i)
  {
    phi0(boundary_edges[i]) = 1;
    if(preserve_source_normals) knownScalValsPhi(i) = 1;
  }

  for(int i=0; i<howManySmoothingInterations; ++i) 
  {
    igl::min_quad_with_fixed
    (Eigen::SparseMatrix<double>(scalarM+t*scalarL), Eigen::VectorXd(-scalarM*phi0), knownScal,
    knownScalValsPhi, Aeq, Beq, false, phit);
    phi0 = phit;
  }

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

void update_phi_accumulatively()
{
  Eigen::SparseMatrix<double> L;
  if(intrinsic)
  {
    Eigen::MatrixXd l; 
    igl::edge_lengths(V, F, l);
    igl::cotmatrix_intrinsic(l, F, L);
  }
  else
  {
    igl::cotmatrix(V, F, L);
  }

  Eigen::SparseMatrix<double> A = -L;
  Eigen::SparseMatrix<double> At = A.transpose();
  Eigen::SparseMatrix<double> AtA = At * A;

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(boundary_vertices.size());
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
    triplets.emplace_back(i, bc[i], 1.);
  }
  Eigen::SparseMatrix<double> B(bc.size(), V.rows());
  B.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::SparseMatrix<double> BtB = B.transpose() * B;

  Div = -0.25*Grad.transpose()*dblA.colwise().replicate(3).asDiagonal();
  Eigen::VectorXd X_vec = Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(X.data(), X.size(), 1);
  Eigen::VectorXd b = -Div * X_vec;

  Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > solver(AtA + lambda * BtB);
  std::cout << "\tsolver info: " << solver.info() << std::endl;
  phi = solver.solve(2.f * At * b).eval();
}

void update_phi()
{
  std::cout << "solving parameters: " << std::endl;
  std::cout << "\tt: " << t << std::endl;
  std::cout << "\tlambda: " << lambda << std::endl;
  std::cout << "\taccumlative: " << accumlative << std::endl;
  std::cout << "\tintrinsic: " << intrinsic << std::endl;
  std::cout << "\tpreserve: " << preserve_source_normals << std::endl;
  std::cout << "\tmin_quad_grad: " << min_quad_grad << std::endl;
  calculateX();
  if (accumlative) return update_phi_accumulatively();

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
  Div = -0.25*Grad.transpose()*dblA.colwise().replicate(3).asDiagonal();
  rhs.head(V.rows()) = Grad.transpose() * X_vec;
  if(min_quad_grad)
    rhs.head(V.rows()) = -Div * X_vec;
  Eigen::SparseMatrix<double> At  = A.transpose();
  Eigen::SparseMatrix<double> AtA = At * A;
  Eigen::VectorXd             Atb = At * rhs;
  Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > solver(AtA);
  std::cout << "\tsolver info: " << solver.info() << std::endl;
  phi = solver.solve(Atb).eval();
}

void save_curve(std::string curve_file)
{
  std::ofstream out(curve_file);
  out << "signed_curve 1\nl ";
  std::for_each(boundary_vertices.begin(), boundary_vertices.end(), [&](int v) {out << v+1 << " ";});
  out.close();
}

void save_label(std::string label_file)
{
  std::ofstream out(label_file);
  out << "something\n";
  for(int f = 0; f < F.rows(); ++f)
  {
    out << (labels(FIt(f)) == label ? 1 : 0) << std::endl;
  }
  out.close();
}

Eigen::MatrixXd U; 
Eigen::MatrixXi G;
Eigen::VectorXd SU;
Eigen::VectorXi J;
Eigen::SparseMatrix<double> UBC;
Eigen::VectorXi GI;

Eigen::MatrixXd phiFC;
void update_phiFC(Eigen::MatrixXd &U, Eigen::MatrixXi &G, Eigen::VectorXd &SU)
{
  phiFC.setZero(G.rows(), 3);
  for(int f = 0; f < G.rows(); ++f)
  {
    double vphi = 0.;
    for(int i = 0; i < 3; ++i)
    {
      vphi += SU(G(f, i));
    }
    vphi /= 3.;
    phiFC.row(f) = vphi > 0 ? red: grey;
  }
}

int main(int argc, char ** argv)
{
  std::string meshname(argv[1]);
  if(!igl::read_triangle_mesh(meshname,VG,FG)) 
  {
    std::cout << "Failed to load mesh " << meshname << std::endl;
    return 0;
  }
  std::vector<std::vector<int>> VGFi;
  igl::vertex_triangle_adjacency(VG, FG, VGF, VGFi);
  igl::pathinfo(meshname, dir, basename, ext, filename);

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

  n_labels = labels.maxCoeff() + 1;
  vertices_type.setZero(VG.rows());
  for(int i = 0; i < VG.rows(); i++)
  {
    if(VGF[i].empty()) continue;
    auto &type = vertices_type(i);
    for(auto f: VGF[i])
    {
      type |= 1 << labels(f);
    }
  }

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
    if (ImGui::CollapsingHeader("Parameters", ImGuiTreeNodeFlags_DefaultOpen))
    {
      static bool accumlative_mode = true;
      if (ImGui::Checkbox("accumlative mode", &accumlative_mode))
      {
        accumlative = accumlative_mode;
      }
      
      static bool intrinsic_mode = false;
      if (ImGui::Checkbox("intrinsic mode", &intrinsic_mode))
      {
        intrinsic = intrinsic_mode;
      }

      static bool grad_mode = true;
      if (ImGui::Checkbox("grad mode", &grad_mode))
      {
        min_quad_grad = grad_mode;
      }
       
      static bool edge_based = false;
      if (ImGui::Checkbox("edge-based Laplace", &edge_based))
      {
        edge_based_Laplace = edge_based;
      }

      static bool preserve = false;
      if (ImGui::Checkbox("preserve source normals", &preserve))
      {
        preserve_source_normals = preserve;
      }

      ImGui::InputDouble("lambda", &lambda, 0, 0, "%.6f");
      ImGui::InputDouble("t", &t, 0, 0, "%.6f");
      ImGui::InputDouble("edge length", &edge_length, 0, 0, "%.4f");
      ImGui::InputDouble("level-set", &level_set, 0, 0, "%.4f");
      ImGui::InputDouble("line width", &line_width, 0, 0, "%.4f");
      ImGui::InputDouble("howMuchToSmoothBy", &howMuchToSmoothBy, 0, 0, "%.4f");
      ImGui::InputInt("howManySmoothingInterations", &howManySmoothingInterations, 0, 0);
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

      static std::string str = filename;
      ImGui::InputText("saving name", str);
      filename = str;
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
          viewer.core().align_camera_center(V, F);
          viewer.data().set_colors(FC);
          viewer.data().set_points(V(boundary_vertices, Eigen::all), green);
          return true;
        case 'x':
          extract_tooth();
          viewer.data().clear();
          viewer.data().set_mesh(V,F);
          viewer.data().set_data(VPV1);
          viewer.data().add_edges(V(boundary_vertices, Eigen::all), V(boundary_vertices, Eigen::all) + VPN(boundary_vertices, Eigen::all), red);
          return true;
        case 'y':
          extract_tooth();
          viewer.data().clear();
          viewer.data().set_mesh(V,F);
          viewer.data().set_data(VPV2);
          viewer.data().add_edges(V(boundary_vertices, Eigen::all), V(boundary_vertices, Eigen::all) + VPD2(boundary_vertices, Eigen::all), red);
          return true;
        case 'z':
          extract_tooth();
          viewer.data().clear();
          viewer.data().set_mesh(V,F);
          viewer.data().set_data(0.5 * VPV1 + 0.5 * VPV2);
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
          viewer.data().set_data(phi);
          viewer.data().add_edges(BC, BC + edge_length*X, XC);
          return true;
        case '3':
          viewer.data().clear();
          update_phi();
          level_sets << level_set;
          igl::remesh_along_isoline(V, F, phi, level_set, U, G, SU, J, UBC, GI);
          viewer.data().set_mesh(U, G);
          update_phiFC(U, G, SU);
          viewer.data().set_colors(phiFC);
          return true;
        case 's':
          if(V.rows() != 0 && F.rows() != 0 && filename != "")
          {
            std::cout << "saving " << dir + "/" + filename + "_" + std::to_string(label) + ".obj" << std::endl;
            igl::writeOBJ(dir + "/" + filename + "_" + std::to_string(label) + ".obj", V, F);
            std::cout << "saving " << dir + "/" + filename + "_" + std::to_string(label) +".txt" << std::endl;
            save_curve(dir + "/" + filename + "_" + std::to_string(label) +".txt");
            std::cout << "saving " << dir + "/" + filename + "_" + std::to_string(label) +".label" << std::endl;
            save_label(dir + "/" + filename + "_" + std::to_string(label) +".label");
          }
          return true;

        default:
          return false;
      }
  };
  viewer.launch();

  return 0;
}
