#include "smoother.h"
#include "remove_unref.h"
#include "extract_kring.h"

#define NOMINMAX

#include <igl/boundary_loop.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/cr_vector_mass.h>
#include <igl/cr_vector_laplacian.h>
#include <igl/crouzeix_raviart_massmatrix.h>
#include <igl/crouzeix_raviart_cotmatrix.h>
#include <igl/edge_midpoints.h>
#include <igl/edge_vectors.h>
#include <igl/principal_curvature.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix_intrinsic.h>
#include <igl/grad.h>
#include <igl/remesh_along_isoline.h>
#include <igl/isolines.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/facet_adjacency_matrix.h>
#include <igl/sum.h>

#include <unordered_map>

using VectorXComplex = Eigen::Matrix<std::complex<double>, -1, 1>;

namespace dgp
{
	namespace alg
	{
		int count_bits(int n, int bits = sizeof(int) * 8)
		{
			assert(n < bits); 
			int count = 0;
			for(int i = 0; i < bits; i++)
			{
				count += n & 1;
				n >>= 1;
			}
			return count; 
		}

		void remove_bit(int &v_label, int label)
		{
			v_label &= ~(1 << label);
		}

		struct SmootherDATA
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXi F;
			std::vector<int> all_bVertices;
			std::vector<std::vector<int>> all_bEdges;
			Eigen::VectorXi faces_label;
			std::vector<Eigen::VectorX<bool>> vertices_label;
			std::vector<std::vector<int> > VF, VFi;
			Eigen::MatrixXi E, oE;
			Eigen::MatrixXi Elist, EMAP;

			Eigen::SparseMatrix<double> vecL, vecM;
			Eigen::SparseMatrix<double> scalarL, scalarM;
			Eigen::MatrixXd edgeMps, paraVec, perpVec;

			bool X_ready{ false };
			bool PHI_ready{ false };
			bool SDF_ready{ false };
			double SDF_levelSet{ 0.2 };
			double SDF_consistencyRatio{ 0.8 };
			Eigen::MatrixXd X;
			Eigen::VectorXd PHI;
			Eigen::VectorXd SDF;

			Eigen::SparseMatrix<double> L;
			Eigen::SparseMatrix<double> Lt, LtL;
			Eigen::SparseMatrix<double> Grad; 
			Eigen::SparseMatrix<double> Div;

			Eigen::MatrixXd C;
			Eigen::MatrixXd lsq, l; 
			Eigen::VectorXd dblA;

			Eigen::VectorXi VJ;

			Eigen::VectorXi BadVertex;

			Eigen::MatrixXd SDFD;
			Eigen::MatrixXd SDFPD; 

			void precompute(
				const std::vector<size_t> &block_vertices,
				const std::vector<std::vector<int>>& local_boundaries,
				const std::vector<bool>& is_border_vertex,
				int block_kring)
			{
				std::set<size_t> blocker(block_vertices.begin(), block_vertices.end());
				const auto nE = E.maxCoeff() + 1;
				Eigen::VectorXi edge_count = Eigen::VectorXi::Zero(nE);
				Eigen::VectorXi vertex_count = Eigen::VectorXi::Zero(V.rows());
				Eigen::VectorXi block = Eigen::VectorXi::Zero(V.rows());
				for (const auto& boundary : local_boundaries)
				{
					std::vector<bool> padding_vertex(boundary.size(), false);
					for (int i = 0; i < boundary.size(); i++)
					{
						int inext = (i + 1) % boundary.size();
						int iprev = (i + boundary.size() - 1) % boundary.size();
						int vip = boundary[iprev];
						int vic = boundary[i];
						int vin = boundary[inext];

						std::unordered_map<int, int> count;
						for (int f : VF[vip]) count[f]++;
						for (int f : VF[vic]) count[f]++;
						for (int f : VF[vin]) 
						{
							if (++count[f] == 3) 
							{
								padding_vertex[i] = true; 
							}
						}
					}

					for (int i = 0; i < boundary.size(); i++)
					{
						int inext = (i + 1) % boundary.size();
						if (inext == 0) 
						{
							continue;
						}
						int vi = boundary[i];
						int vj = boundary[inext];

						if (is_border_vertex[VJ(vi)])
						{
							int prev = i;
							int next = i;
							for(int k = 0; k < block_kring; ++k)
							{
								block(boundary[prev]) += 1;
								block(boundary[next]) += 1;
								prev = (prev + boundary.size() - 1) % boundary.size();
								next = (next + 1) % boundary.size();
							}
						}

						vertex_count(vi)++;
						if (padding_vertex[i] || padding_vertex[inext]) 
						{
							continue;
						}
						if (blocker.find(vi) != blocker.end())
						{
							continue; 
						}

						int eid = -1;
						for (int f : VF[vi])
						{
							for (int j = 0; j < 3; ++j)
							{
								if (F(f, j) == vj)
								{
									for (int k = 0; k < 3; ++k)
									{
										if (F(f, k) != vi && F(f, k) != vj)
										{ 
											eid = E(f, k);
											if (edge_count(eid)++ > 0) 
											{
												continue;
											}
											all_bEdges.push_back(std::vector<int>{f, k, vi, vj});			
											all_bVertices.push_back(vi);
											all_bVertices.push_back(vj);
											break;
										}
									}
								}
							}
							if (eid != -1) break;
						} 
					} 
				}
				std::sort(all_bVertices.begin(), all_bVertices.end()); 
				all_bVertices.erase(std::unique(all_bVertices.begin(), all_bVertices.end()), all_bVertices.end());
				BadVertex = vertex_count.unaryExpr([](int count) { return count > 1 ? 1 : 0; });
				BadVertex += block; 
			} 

			void precompute(bool use_intrinsic)
			{
				// Warning: E and oE might not be recompute
				igl::cr_vector_mass(V, F, E, vecM);
				igl::cr_vector_laplacian(V, F, E, oE, vecL);

				const int nE = E.maxCoeff() + 1;
				assert(nE == vecL.rows() / 2);

				Elist.setZero(nE, 2);
				EMAP.setZero(3 * F.rows(), 1);
				for (int i = 0; i < F.rows(); ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						const int e = E(i, j);
						EMAP(i + j * F.rows()) = e;
						if (oE(i, j) > 0)
						{
							Elist.row(e) << F(i, (j + 1) % 3), F(i, (j + 2) % 3);
						}
					}
				}

				igl::crouzeix_raviart_massmatrix(V, F, Elist, EMAP, scalarM);
				igl::crouzeix_raviart_cotmatrix(V, F, Elist, EMAP, scalarL);

				igl::edge_midpoints(V, F, E, oE, edgeMps);
				igl::edge_vectors(V, F, E, oE, paraVec, perpVec);

				igl::vertex_triangle_adjacency(V, F, VF, VFi);

				igl::squared_edge_lengths(V, F, lsq);
				l = lsq.array().sqrt().eval();

				if (use_intrinsic)
				{
					igl::cotmatrix_intrinsic(l, F, L);
				}
				else
				{
					igl::cotmatrix(V, F, L);
				}
				igl::cotmatrix_entries(l, C);
				L *= -1.;
				Lt = L.transpose();
				LtL = Lt * L;
				igl::grad(V, F, Grad);
				igl::doublearea(V, F, dblA);
				Div = -0.25 * Grad.transpose() * dblA.colwise().replicate(3).asDiagonal();
			}

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

			void complex_map(const VectorXComplex &Xt, Eigen::MatrixXd &Xf)
			{
				Xf = Eigen::MatrixXd::Zero(F.rows(), 3);
				for (int f = 0; f < F.rows(); ++f)
				{
					double l_jk_sq = lsq(f, 0);
					double l_ki_sq = lsq(f, 1);
					double l_ij_sq = lsq(f, 2);
					double ai = l_ij_sq + l_ki_sq - l_jk_sq;
					double aj = l_jk_sq + l_ij_sq - l_ki_sq;
					double ak = l_ki_sq + l_jk_sq - l_ij_sq;
					Eigen::RowVector3d bc = { ai, aj, ak };
					double s = 0.5 * (l(f, 0) + l(f, 1) + l(f, 2));
					double ar = std::sqrt(s * (s - l(f, 0)) * (s - l(f, 1)) * (s - l(f, 2)));
					if(std::isnan(ar) || std::isinf(ar) || ar < 0.)
					{
						ar = 0.5 * dblA(f);
					}

					Eigen::RowVector3d fv = Eigen::RowVector3d::Zero();
					for (int i = 0; i < 3; ++i)
					{
						Eigen::RowVector3d e1 = edge_basis(i);
						e1 *= static_cast<double>(oE(f, i));

						Eigen::RowVector3d e2{ e1(2) * bc(2) - e1(1) * bc(1),
											  e1(0) * bc(0) - e1(2) * bc(2),
											  e1(1) * bc(1) - e1(0) * bc(0) };
						e2 /= (4. * ar);
						double norm_e1 = std::sqrt(-((l_ij_sq * e1(0) * e1(1) +
							l_jk_sq * e1(1) * e1(2) +
							l_ki_sq * e1(2) * e1(0))));
						e1 /= norm_e1;

						double norm_e2 = std::sqrt(-((l_ij_sq * e2(0) * e2(1) +
							l_jk_sq * e2(1) * e2(2) +
							l_ki_sq * e2(2) * e2(0))));
						e2 /= norm_e2;
						fv += std::real(Xt(E(f, i))) * e1;
						fv += std::imag(Xt(E(f, i))) * e2;
					}
					double norm_fv = std::sqrt(-((l_ij_sq * fv(0) * fv(1) +
						l_jk_sq * fv(1) * fv(2) +
						l_ki_sq * fv(2) * fv(0))));
					if (norm_fv == 0)
					{
						Xf.row(f) = Eigen::RowVector3d::Zero();
					}
					else
					{
						Xf.row(f) = fv / norm_fv;
					}
				}
			}

			void vector_field_smoothing(const Eigen::MatrixXd& Xf, Eigen::MatrixXd& Xv)
			{
				Xv.setZero(V.rows(), 3);
				for (int v = 0; v < V.rows(); ++v)
				{
					for (auto f : VF[v])
					{
						Xv.row(v) += Xf.row(f);
					}
					Xv.row(v) /= static_cast<double>(VF[v].size());
				}
			}

			void compute_Div(const Eigen::MatrixXd& Xd, Eigen::VectorXd &Div)
			{
				Div = Eigen::VectorXd::Zero(V.rows());
				for (int v = 0; v < V.rows(); ++v)
				{
					for (int f : VF[v])
					{
						double l_jk_sq = lsq(f, 0);
						double l_ki_sq = lsq(f, 1);
						double l_ij_sq = lsq(f, 2);
						for (int i = 0; i < 3; ++i)
						{
							if (F(f, i) == v)
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
			}

			void shift_sdf()
			{
				for (int i = 0; i < all_bVertices.size(); ++i)
				{
					SDF(all_bVertices[i]) = 0.;
				}

				double shift = 0.;
				double normalization = 0.;
				for (const auto& edge : all_bEdges)
				{
					int vi = edge[2];
					int vj = edge[3];
					for (int f : VF[vi])
					{
						int i_ = -1;
						for (int j = 0; j < 3; ++j)
						{
							if (F(f, j) == vi)
							{
								i_ = j;
								break;
							}
						}

						int j_ = -1;
						for (int j = 0; j < 3; ++j)
						{
							if (F(f, j) == vj)
							{
								j_ = j;
								break;
							}
						}

						int k_ = -1;
						for (int j = 0; j < 3; ++j)
						{
							if (F(f, j) != vi && F(f, j) != vj)
							{
								k_ = j;
								break;
							}
						}

						if (i_ == -1 || j_ == -1 || k_ == -1)
						{
							continue;
						}

						Eigen::RowVector3d vi_coords;
						vi_coords(i_) = 1; vi_coords(j_) = 0; vi_coords(k_) = 0;
						Eigen::RowVector3d vj_coords;
						vj_coords(i_) = 0; vj_coords(j_) = 1; vj_coords(k_) = 0;
						double l_ij = l(f, k_);
						Eigen::RowVector3d midp = 0.5 * (vi_coords + vj_coords);
						Eigen::RowVector3d fsdf;
						fsdf(i_) = SDF(F(f, i_)); fsdf(j_) = SDF(F(f, j_)); fsdf(k_) = SDF(F(f, k_));
						shift += l_ij * midp.dot(fsdf);
						normalization += l_ij;
						break;
					}
				}

				shift /= normalization;
				SDF -= (Eigen::VectorXd::Ones(V.rows()) * shift);
			}

			void solve_sdf(SDFParam param)
			{
				SDF_levelSet = param.level_set;
				SDF_consistencyRatio = param.consistency_ratio;
				SDF_ready = true; 
				const auto nE = E.maxCoeff() + 1; 
				VectorXComplex X0 = VectorXComplex::Zero(nE);
				std::vector<int> constriants; constriants.reserve(all_bEdges.size());
				for(const auto &edge: all_bEdges)
				{
					int f = edge[0];
					int k = edge[1];
					int vi = edge[2];
					int eid = E(f, k);
					std::complex<double> innerProd(0, Elist(eid, 0) == vi ? 1 : -1);
					std::complex<double> contrib = l(f, k) * innerProd;
					X0[eid] += contrib;
					constriants.push_back(eid);
				}

				VectorXComplex Xt;
				Eigen::VectorXd Y0 = Eigen::VectorXd::Zero(2 * nE), Yt;
				for (int i = 0; i < nE; ++i)
				{
					Y0(i) = std::real(X0(i));
					Y0(i + nE) = std::imag(X0(i));
				}

				Eigen::VectorXi known = Eigen::VectorXi::Zero(2 * constriants.size());
				Eigen::VectorXd known_vals = Eigen::VectorXd::Zero(2 * constriants.size());
				for (int i = 0; i < constriants.size(); ++i)
				{
					known(i) = constriants[i];
					known(i + constriants.size()) = constriants[i] + nE;
					known_vals(i) = std::real(X0(constriants[i]));
					known_vals(i + constriants.size()) = std::imag(X0(constriants[i]));
				}

				Eigen::SparseMatrix<double> Aeq;
				Eigen::VectorXd Beq;
				SDF_ready &= igl::min_quad_with_fixed(Eigen::SparseMatrix<double>(vecM + param.t * vecL),
					Eigen::VectorXd(-vecM * Y0), known, known_vals, Aeq, Beq, false, Yt);
				Xt.setZero(nE);
				for (int i = 0; i < nE; ++i)
				{
					Xt(i) = std::complex<double>(Yt(i), Yt(i + nE));
				}

				Eigen::MatrixXd Xf;
				complex_map(Xt, Xf);

				Eigen::VectorXd Div;
				compute_Div(Xf, Div);

				Eigen::MatrixXd Xv;
				vector_field_smoothing(Xf, Xv);

				known = Eigen::VectorXi::Zero(all_bVertices.size());
				known_vals = Eigen::VectorXd::Zero(all_bVertices.size());
				for(int i = 0; i < all_bVertices.size(); ++i)
				{
					known(i) = all_bVertices[i];
					known_vals(i) = Div(all_bVertices[i]);
				}
				SDF_ready &= igl::min_quad_with_fixed(L, -Div, known, known_vals, Aeq, Beq, false, SDF);
				shift_sdf();
			}

			bool isolines_test(
				const Eigen::MatrixXd &iV, 
				const Eigen::MatrixXi &iE)
			{
				if (iV.rows() == 0 || iE.rows() == 0) return false; 
				Eigen::VectorXi vnv_count = Eigen::VectorXi::Zero(iV.rows());
				std::for_each(iE.data(), iE.data() + iE.size(), [&vnv_count](int v)
				{
					vnv_count(v)++;
				});
				return (vnv_count.array() == 2).all();
			}

			bool isolines_test(const Eigen::VectorXd& VALUE, double level_set)
			{
				Eigen::VectorXd level_sets; level_sets.resize(1);
				level_sets << level_set;
				Eigen::MatrixXd iV;
				Eigen::MatrixXi iE;
				Eigen::VectorXi iI;
				igl::isolines(V, F, VALUE, level_sets, iV, iE, iI);
				return isolines_test(iV, iE);
			}

			bool sdf_boundary_condition(
				std::vector<int> &edge_constraints,
				std::vector<double> &paraVals,
				std::vector<double> &perpVals,
				int min_spacing_source_normals,
				const std::vector<size_t> &block_vertices)
			{
				if (!SDF_ready) return false;
				edge_constraints.clear();
				paraVals.clear();
				perpVals.clear();

				bool success_isolines_test = false;
				Eigen::MatrixXd iV;
				Eigen::MatrixXi iE;
				Eigen::VectorXi iI;
				for(int i = 0; i < 3; ++i)
				{
					Eigen::VectorXd level_sets; level_sets.resize(1);
					level_sets << (average_value_of_labeled_faces(SDF) > 0. ? 1. : -1.) * (SDF_levelSet + static_cast<double>(i) * 0.1);
					igl::isolines(V, F, SDF, level_sets, iV, iE, iI);
					if (isolines_test(iV, iE))
					{
						success_isolines_test = true;
						break;
					}
				}

				if (!success_isolines_test) return false; 

				std::vector<int> all_eVertices;
				for (const auto& edge : all_bEdges)
				{
					all_eVertices.push_back(edge[2]);
					all_eVertices.push_back(edge[3]);
				}

				Eigen::VectorXd sqD;
				Eigen::VectorXi sqDI;
				Eigen::MatrixXd sqDC;
				igl::point_mesh_squared_distance(V(all_eVertices, Eigen::all).eval(), iV, iE, sqD, sqDI, sqDC);

				std::vector<int> bad_vertices;
				Eigen::MatrixXd PD1, PD2;
				Eigen::VectorXd PV1, PV2;
				igl::principal_curvature(V, F, PD1, PD2, PV1, PV2, bad_vertices);
				std::for_each(bad_vertices.begin(), bad_vertices.end(), [&](int v) { BadVertex(v) += 1; });

				Eigen::VectorXi Blocker = Eigen::VectorXi::Zero(BadVertex.size());
				std::for_each(block_vertices.begin(), block_vertices.end(), [&](auto v){ Blocker(v) += 1; });

				SDFD.setZero(V.rows(), 3);
				SDFPD.setZero(V.rows(), 3);
				std::vector<bool> is_sampled(all_bEdges.size(), false);
				for(int i = 0; i < all_bEdges.size(); ++i)
				{
					const auto& edge = all_bEdges[i];
					int e = E(edge[0], edge[1]);
					int vi = edge[2];
					int vj = edge[3];
					if (BadVertex(vi) > 0 || BadVertex(vj) > 0)
					{
						is_sampled[i] = true; 
						continue;
					}
					if (Blocker(vi) > 0 || Blocker(vj) > 0)
					{
						is_sampled[i] = false; 
						continue;
					}

					Eigen::Vector3d pd;
					Eigen::Vector3d sdfd;
					double dot;
					{
						Eigen::Vector3d vi_pd;
						Eigen::Vector3d vi_sdfd;
						double vi_dot;
						{
							Eigen::Vector3d sdf_dir = sqDC.row(2 * i + 0) - V.row(vi); sdf_dir.normalize();
							Eigen::Vector3d maxpd = PD1.row(vi); maxpd.normalize();
							Eigen::Vector3d minpd = PD2.row(vi); minpd.normalize();
							Eigen::Vector3d maxpd_verify = sdf_dir.cross(minpd); maxpd_verify.normalize();
							vi_dot = std::abs(std::max(std::min(maxpd_verify.dot(maxpd), 1.), -1.));
							vi_pd = minpd;
							vi_sdfd = sdf_dir;
						}

						Eigen::Vector3d vj_pd;
						Eigen::Vector3d vj_sdfd;
						double vj_dot;
						{
							Eigen::Vector3d sdf_dir = sqDC.row(2 * i + 1) - V.row(vj); sdf_dir.normalize();
							Eigen::Vector3d maxpd = PD1.row(vj); maxpd.normalize();
							Eigen::Vector3d minpd = PD2.row(vj); minpd.normalize();
							Eigen::Vector3d maxpd_verify = sdf_dir.cross(minpd); maxpd_verify.normalize();
							vj_dot = std::abs(std::max(std::min(maxpd_verify.dot(maxpd), 1.), -1.));
							vj_pd = minpd;
							vj_sdfd = sdf_dir;

						}

						if(vi_dot > vj_dot)
						{
							dot = vi_dot;
							pd = vi_pd;
							sdfd = vi_sdfd;
						}
						else
						{
							dot = vj_dot;
							pd = vj_pd;
							sdfd = vj_sdfd;
						}

						SDFD.row(vi) = vi_sdfd;
						SDFD.row(vj) = vj_sdfd;

						if (dot < SDF_consistencyRatio)
						{
							continue;
						}

						SDFPD.row(vi) = vi_dot * vi_pd;
						SDFPD.row(vj) = vj_dot * vj_pd;
					}

					is_sampled[i] = true;
					pd *= (pd.dot(sdfd) > 0. ? 1. : -1.);

					Eigen::Vector3d para = paraVec.row(e); para.normalize();
					Eigen::Vector3d perp = perpVec.row(e); perp.normalize();
					Eigen::Vector3d n = para.cross(perp); n.normalize();
					pd = pd - pd.dot(n) * n;

					paraVals.push_back(pd.dot(para));
					perpVals.push_back(pd.dot(perp));
					edge_constraints.push_back(e);
				}
				
				for(size_t i = 0; i < is_sampled.size(); )
				{
					if(is_sampled[i])
					{
						size_t j = (i + 1) % is_sampled.size();
						int count = 0;
						while (!is_sampled[j]) 
						{
							j = (j + 1) % is_sampled.size();
							count++; 
						}

						int n = count / min_spacing_source_normals;
						int interval = count / (n + 1) + (count % 2 == 0 ? 0 : 1);
						for(size_t k = 1; k < n + 1; ++k)
						{
							size_t index = (i + k * interval) % is_sampled.size();
							const auto& edge = all_bEdges[index];
							int e = E(edge[0], edge[1]);
							int vi = edge[2];
							int vj = edge[3];
							//TODO: sdf_dir = arg min E(vi, vj)
							Eigen::Vector3d sdf_dir = sqDC.row(2 * index + 0) - V.row(vi); sdf_dir.normalize();

							Eigen::Vector3d para = paraVec.row(e); para.normalize();
							Eigen::Vector3d perp = perpVec.row(e); perp.normalize();
							Eigen::Vector3d n = para.cross(perp); n.normalize();
							sdf_dir = sdf_dir - sdf_dir.dot(n) * n;
							
							SDFPD.row(vi) = sdf_dir;

							paraVals.push_back(sdf_dir.dot(para));
							perpVals.push_back(sdf_dir.dot(perp));
							edge_constraints.push_back(e);
						}

						if(j <= i)
						{
							break;
						}
						i = j;
					}
					else
					{
						++i;
					}
				}

				return true;
			}

			void solve_X(const VectorTransportParam &param)
			{
				X_ready = false;
				std::vector<int> edge_constraints; edge_constraints.reserve(all_bEdges.size());
				std::for_each(all_bEdges.begin(), all_bEdges.end(), [&](const auto& edge)
				{
					edge_constraints.push_back(E(edge[0], edge[1]));
				});
				std::vector<double> paraVals = std::vector<double>(edge_constraints.size(), 0.08);
				std::vector<double> perpVals = std::vector<double>(edge_constraints.size(), 0.98);
				if (!sdf_boundary_condition(edge_constraints, paraVals, perpVals, param.min_spacing_source_normals, param.block_vertices))
				{
					return;
				}
				X_ready = true;

				const int nE = E.maxCoeff() + 1; 
				Eigen::VectorXi known;
				Eigen::VectorXd knownVals;
				if (param.preserve_source_normals)
				{
					known = Eigen::VectorXi::Zero(2 * edge_constraints.size());
					knownVals = Eigen::VectorXd::Zero(2 * edge_constraints.size());
				}

				Eigen::VectorXd Y0 = Eigen::VectorXd::Zero(2 * nE), Yt;
				for (int i = 0; i < edge_constraints.size(); ++i)
				{
					if (param.preserve_source_normals) known(i) = edge_constraints[i];
					if (param.preserve_source_normals) known(i + edge_constraints.size()) = edge_constraints[i] + nE;
					if (param.preserve_source_normals) knownVals(i) = paraVals[i];
					if (param.preserve_source_normals) knownVals(i + edge_constraints.size()) = perpVals[i];
					Y0(edge_constraints[i]) = paraVals[i];
					Y0(edge_constraints[i] + nE) = perpVals[i];
				}

				Eigen::SparseMatrix<double> Aeq;
				Eigen::VectorXd Beq;
				X_ready &= igl::min_quad_with_fixed
				(Eigen::SparseMatrix<double>(vecM + param.t * vecL), Eigen::VectorXd(-vecM * Y0), known, knownVals,
					Aeq, Beq, false, Yt);

				Eigen::VectorXi knownScal;
				Eigen::VectorXd knownScalVals;
				if (param.preserve_source_normals)
				{
					knownScal = Eigen::VectorXi::Zero(edge_constraints.size());
					knownScalVals = Eigen::VectorXd::Zero(edge_constraints.size());
				}

				Eigen::VectorXd u0 = Eigen::VectorXd::Zero(nE), ut;
				for (int i = 0; i < edge_constraints.size(); ++i)
				{
					u0(edge_constraints[i]) = sqrt(paraVals[i] * paraVals[i] + perpVals[i] * perpVals[i]);
					if (param.preserve_source_normals) knownScal(i) = edge_constraints[i];
					if (param.preserve_source_normals) knownScalVals(i) = u0(edge_constraints[i]);
				}

				X_ready &= igl::min_quad_with_fixed
				(Eigen::SparseMatrix<double>(scalarM + param.t * scalarL), Eigen::VectorXd(-scalarM * u0), knownScal,
					knownScalVals, Aeq, Beq, false, ut);

				Eigen::VectorXd knownScalValsPhi;
				if (param.preserve_source_normals)
				{
					knownScalValsPhi = Eigen::VectorXd::Zero(edge_constraints.size());
				}

				Eigen::VectorXd phi0 = Eigen::VectorXd::Zero(nE), phit;
				for (int i = 0; i < edge_constraints.size(); ++i)
				{
					phi0(edge_constraints[i]) = 1;
					if (param.preserve_source_normals) knownScalValsPhi(i) = 1;
				}

				X_ready &= igl::min_quad_with_fixed
				(Eigen::SparseMatrix<double>(scalarM + param.t * scalarL), Eigen::VectorXd(-scalarM * phi0), knownScal,
					knownScalValsPhi, Aeq, Beq, false, phit);

				Eigen::ArrayXd Xtfactor = ut.array() /
					(phit.array() * (Yt.array().segment(0, nE) * Yt.array().segment(0, nE)
						+ Yt.array().segment(nE, nE) * Yt.array().segment(nE, nE)).sqrt());
				Eigen::VectorXd Xt(2 * nE);
				Xt.segment(0, nE) = Xtfactor * Yt.segment(0, nE).array();
				Xt.segment(nE, nE) = Xtfactor * Yt.segment(nE, nE).array();

				Eigen::MatrixXd vecs(nE, 3);
				for (int i = 0; i < edgeMps.rows(); ++i)
				{
					vecs.row(i) = Xt(i) * paraVec.row(i) + Xt(i + edgeMps.rows()) * perpVec.row(i);
				}

				X.setZero(F.rows(), 3);
				for (int i = 0; i < F.rows(); ++i)
				{
					for (int j = 0; j < F.cols(); ++j)
					{
						X.row(i) += vecs.row(E(i, j));
					}
					X.row(i) /= static_cast<double>(F.cols());
				}
				X.rowwise().normalize();
			}

			void solve_PHI(const DiffusionParam& param)
			{
				PHI_ready = false;
				if (!X_ready) return;

				std::set<size_t> blocker;
				std::for_each(param.block_vertices.begin(), param.block_vertices.end(), [&blocker](auto pair){ blocker.insert(pair.first); });

				std::vector<std::pair<int, double> > constraints;
				constraints.reserve(all_bVertices.size() + param.potential_vertices.size()); 
				for (int i = 0; i < all_bVertices.size(); i++)
				{
					if (i % param.min_spacing_constriants == 0 &&
						!BadVertex(all_bVertices[i]) && 
						blocker.find(all_bVertices[i]) == blocker.end())
					{
						auto vertex = all_bVertices[i];
						auto vertex_label = vertices_label[vertex];
						vertex_label(0) = false; 
						constraints.emplace_back(std::make_pair(vertex, vertex_label.count() > 1? param.adjacent_point_potential: 0.));
					}
				}
				constraints.insert(constraints.end(), param.potential_vertices.begin(), param.potential_vertices.end());

				std::vector<Eigen::Triplet<double>> triplets;
				triplets.reserve(constraints.size());
				for (int i = 0; i < constraints.size(); i++)
				{
					triplets.emplace_back(i, constraints[i].first, 1.);
				}

				Eigen::SparseMatrix<double> B(constraints.size(), V.rows());
				B.setFromTriplets(triplets.begin(), triplets.end());
				Eigen::SparseMatrix<double> Bt = B.transpose();
				Eigen::SparseMatrix<double> BtB = Bt * B;

				Eigen::MatrixXd Dense_uL = Eigen::MatrixXd::Zero(V.rows(), V.rows());
				for(auto block_vertex: param.block_vertices)
				{
					auto v = block_vertex.first;
					auto nvf = VF[v].size(); 
					for(int i = 0; i < nvf; ++i)
					{
						int f = VF[v][i];
						int _1 = VFi[v][i]; 
						int _0 = (_1 - 1 + 3) % 3; 
						int _2 = (_1 + 1) % 3; 
						Dense_uL(v, F(f, _0)) = 1.; 
						Dense_uL(v, F(f, _2)) = 1.; 
					}
					Dense_uL(v, v) = -(double)nvf;
				}

				Eigen::SparseMatrix<double> uL = Dense_uL.sparseView();

				Eigen::VectorXd DenseW = Eigen::VectorXd::Ones(V.rows());
				for(auto pair: param.block_vertices)
				{
					DenseW(pair.first) = pair.second * pair.second;
				}
				Eigen::SparseMatrix<double> W = DenseW.asDiagonal().toDenseMatrix().sparseView();

				Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > solver(Lt * W * L + param.lambda_regularization * uL.transpose() * uL + param.lambda * BtB);
				if (solver.info() != Eigen::Success) return;

				Eigen::VectorXd X_vec = Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(X.data(), X.size(), 1);
				Eigen::VectorXd l = -Div * X_vec;

				Eigen::VectorXd b = Eigen::VectorXd::Zero(constraints.size());
				for (int i = 0; i < constraints.size(); i++)
				{
					b(i) = constraints[i].second;
				}
				PHI = solver.solve(2. * Lt * W * l + 2. * Bt * b).eval();
				if (solver.info() != Eigen::Success) return;
				PHI_ready = true;
			}

			Eigen::VectorXd get_PHI(Eigen::Index number_of_vertices)
			{
				//TODO: 
				Eigen::VectorXi boundary;
				igl::boundary_loop(F, boundary);
				double phi = PHI(boundary).mean();
				Eigen::VectorXd Phi_out; Phi_out.setConstant(number_of_vertices, phi);
				for(int v = 0; v < V.rows(); ++v)
				{
					Phi_out(VJ(v)) = PHI(v);
				}
				return  Phi_out;
			}

			double average_value_of_labeled_faces(const Eigen::VectorXd &VALUE)
			{
				double phi = 0.;
				for(int f = 0; f < F.rows(); ++f)
				{
					double f_phi = 0.;
					if (static_cast<bool>(faces_label(f)))
					{
						for (int i = 0; i < 3; ++i)
						{
							f_phi += VALUE(F(f, i));
						}
					}
					phi += f_phi / 3.;
				}
				return phi;
			}
		};

		class SmootherImpl
		{
		public:
			SmootherImpl() = default;
			SmootherImpl(const SmootherImpl&) = delete;
			SmootherImpl& operator=(const SmootherImpl&) = delete;
			~SmootherImpl() = default;

			void precompute(const std::vector<double>& v_buffer, const std::vector<int>& f_buffer)
			{
				V_ = Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(v_buffer.data(), v_buffer.size() / 3, 3);
				F_ = Eigen::Map<const Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(f_buffer.data(), f_buffer.size() / 3, 3);

				igl::vertex_triangle_adjacency(V_, F_, VF_, VFi_);
			}
			
			void precompute(const std::vector<int>& labels, PrecomputeParam& param)
			{
				is_border_vertex_ = igl::is_border_vertex(F_);
				faces_label_ = Eigen::Map<const Eigen::VectorXi>(labels.data(), labels.size());
				vertices_label_.resize(V_.rows(), Eigen::VectorX<bool>::Zero(faces_label_.maxCoeff() + 1));
				for(int v = 0; v < V_.rows(); ++v)
				{
					auto &vextex_label_ = vertices_label_[v];
					for (auto f : VF_[v]) 
					{
						vextex_label_(faces_label_(f)) = true; 
					}
				}
				data_ptrs_.clear();
				data_ptrs_ = std::vector<std::unique_ptr<SmootherDATA>>(faces_label_.maxCoeff() + 1);
				std::set<int> label_set(faces_label_.begin(), faces_label_.end());
				for(auto label: label_set)
				{
					if (label == 0) continue;
					data_ptrs_[label] = std::make_unique<SmootherDATA>();
					std::vector<std::vector<int>> label_boundaries;
					extracting(label, param.fuzzy_kring, data_ptrs_[label]->V, data_ptrs_[label]->F, label_boundaries, data_ptrs_[label]->faces_label, data_ptrs_[label]->vertices_label, data_ptrs_[label]->VJ);
					data_ptrs_[label]->precompute(param.use_intrinsic);
					data_ptrs_[label]->precompute(param.block_vertices, label_boundaries, is_border_vertex_, param.block_kring);
					if (data_ptrs_[label]->all_bEdges.empty() || data_ptrs_[label]->all_bVertices.empty()) 
					{
						std::cout << "[SmootherImpl::precompute unexpected error: extracting label " << label << "  failed]\n";
						data_ptrs_[label].reset(nullptr);
					}
				}
			}
			
			bool valid_labels(const std::vector<int>& labels)
			{
				return labels.size() == F_.rows();
			}

			void extension(const Eigen::VectorXi& v_candidates, int kring, Eigen::VectorXi& K)
			{
				std::vector<int> rings;
				for (auto v : v_candidates)
				{
					std::vector<int> vrings, frings;
					extract_kring(F_, VF_, v, kring, vrings, frings);
					rings.insert(rings.end(), frings.begin(), frings.end());
				}
				std::for_each(rings.begin(), rings.end(), [&K](auto f) { K(f) = 1; });
			}

			bool extracting(
				int label,
				int kring,
				Eigen::MatrixXd& Vo, Eigen::MatrixXi& Fo,
				std::vector<std::vector<int>>& boundaries,
				Eigen::VectorXi& faces_label,
				std::vector<Eigen::VectorX<bool>> &vertices_label,
				Eigen::VectorXi& VJ)
			{
				Eigen::VectorXi Ki = faces_label_.unaryExpr([&label](auto fl) { return fl == label ? 1 : 0; });

				Eigen::MatrixXd Vi;
				Eigen::MatrixXi Fi;
				Eigen::VectorXi Ii, Ji, FiMAP;
				remove_unref(V_, F_, Ki, Vi, Fi, Ii, Ji, FiMAP);
				if(Fi.rows() == 0 || Vi.rows() == 0)
				{
					return false; 
				}

				// this method might not get all boundaries
				Eigen::VectorXi all_bVertices;
				{
					std::vector<std::vector<int>> tmp_boundaries;
					igl::boundary_loop(Fi, tmp_boundaries);
					std::vector<int> tmp_all_bVertices;
					for (const auto& tmp_boundary : tmp_boundaries)
					{
						for (int v : tmp_boundary)
						{
							tmp_all_bVertices.push_back(Ji(v));
						}
					}
					all_bVertices = Eigen::Map<Eigen::VectorXi>(tmp_all_bVertices.data(), tmp_all_bVertices.size());
				}
				extension(all_bVertices, kring, Ki);

				Eigen::VectorXi Ii_ex, Ji_ex, Fi_exMAP;
				remove_unref(V_, F_, Ki, Vo, Fo, Ii_ex, Ji_ex, Fi_exMAP);
				faces_label = Eigen::VectorXi::LinSpaced(Fo.rows(), 0, Fo.rows()).unaryExpr([&](int f) {return static_cast<int>(faces_label_(Fi_exMAP(f)) == label); });
				VJ = Ji_ex;

				Eigen::VectorXi Ki_ex = Eigen::VectorXi::LinSpaced(Fo.rows(), 0, Fo.rows()).unaryExpr([&](int f) {return faces_label_(Fi_exMAP(f)) == label ? 1 : 0; });
				remove_unref(Vo, Fo, Ki_ex, Vi, Fi, Ii, Ji, FiMAP);

				igl::boundary_loop(Fi, boundaries);
				for(auto &boundary: boundaries)
				{
					for(auto &v: boundary)
					{
						v = Ji(v);
					}
				}
				
				vertices_label.clear(); 
				vertices_label.reserve(Vo.rows());
				for(Eigen::Index v = 0; v < Vo.rows(); ++v)
				{
					vertices_label.emplace_back(vertices_label_[VJ(v)]);
				}
				return true;
			}

			void sdf_required(SDFParam param)
			{
				for (int label = 0; label < data_ptrs_.size(); ++label)
				{
					if (data_ptrs_[label] == nullptr) continue;
					auto& data_ptr = data_ptrs_[label];
					data_ptr->solve_sdf(param);
				}
			}

			void vector_field_solve(const VectorTransportParam& param)
			{
				for (int label = 0; label < data_ptrs_.size(); ++label)
				{
					if (data_ptrs_[label] == nullptr) continue;
					auto& data_ptr = data_ptrs_[label];
					data_ptr->solve_X(param);
				}
			}

			void scalar_field_solve(const DiffusionParam& param)
			{
				for(int label = 0; label < data_ptrs_.size(); ++label)
				{
					if (data_ptrs_[label] == nullptr) continue;
					auto& data_ptr = data_ptrs_[label];
					data_ptr->solve_PHI(param);
				}
			}

			void relaxing(
				double collapse_lenRatio,
				Eigen::MatrixXd& V,
				Eigen::MatrixXi& F,
				Eigen::VectorXi& faces_label,
				Eigen::VectorXi& J)
			{
				const auto isCommonVertex = [](
					const Eigen::VectorXi& vertices_label,
					int vi)->bool
				{
					return count_bits(vertices_label(vi)) > 2;
				};

				const auto commonValues = [](
					const std::vector<int>& v0,
					const std::vector<int>& v1)->std::vector<int>
				{
					std::set<int> set0(v0.begin(), v0.end());
					std::set<int> result;
					for (int val : v1)
					{
						if (set0.find(val) != set0.end())
						{
							result.insert(val);
						}
					}
					return std::vector<int>(result.begin(), result.end());
				};

				const auto commonVertices = [&commonValues](
					const Eigen::MatrixXi& F,
					const Eigen::VectorXi& fKeep,
					const std::vector<int>& vf1,
					const std::vector<int>& vf2)->std::vector<int>
				{
					auto vf1_ = vf1;
					auto vf2_ = vf2;
					vf1_.erase(
						std::remove_if(
							vf1_.begin(),
							vf1_.end(),
							[&fKeep](int f) { return !fKeep(f); }),
						vf1_.end());
					vf2_.erase(
						std::remove_if(
							vf2_.begin(),
							vf2_.end(),
							[&fKeep](int f) { return !fKeep(f); }),
						vf2_.end());
					std::vector<int> v1, v2;
					for (int f : vf1_)
					{
						for (int i = 0; i < 3; ++i) v1.push_back(F(f, i));
					}
					for (int f : vf2_)
					{
						for (int i = 0; i < 3; ++i) v2.push_back(F(f, i));
					}
					return commonValues(v1, v2);
				};

				const auto topoTest = [&commonVertices](
					const Eigen::MatrixXi& F,
					const std::vector<std::vector<int>> &VF,
					const Eigen::VectorXi& fKeep,
					int v0,
					int v1)->bool
				{
					const auto cmv = commonVertices(F, fKeep, VF[v0], VF[v1]);
					if (cmv.size() != 4) return false;
					int cnt_false = 0;
					int cnt_true = 0;
					for(int v: cmv)
					{
						if (v == v0 || v == v1) cnt_true++;
						else cnt_false++;
					}
					return cnt_true == 2 && cnt_true == cnt_false;
				};

				const auto filppedTest = [](
					const Eigen::Vector3d& q,
					const Eigen::MatrixXd& G,
					const Eigen::MatrixXi& T,
					const Eigen::VectorXi& fKeep,
					const std::set<int>& common_faces,
					const std::vector<int>& vf,
					const std::vector<int>& vfi)->bool
				{
					for (int iter = 0; iter < vf.size(); ++iter)
					{
						const int f = vf[iter];
						const int i = vfi[iter];
						if (!fKeep(f)) continue;
						if (common_faces.find(f) != common_faces.end()) continue;
						const int j = (i + 1) % 3;
						const int k = (i + 2) % 3;
						const Eigen::Vector3d u = (Eigen::Vector3d(G.row(T(f, j))) - q).normalized();
						const Eigen::Vector3d v = (Eigen::Vector3d(G.row(T(f, k))) - q).normalized();
						if (std::fabs(u.dot(v)) > 0.999) return false;
						const Eigen::Vector3d n = u.cross(v).normalized();
						const Eigen::Vector3d u_ = (Eigen::Vector3d(G.row(T(f, j))) - Eigen::Vector3d(G.row(T(f, i)))).normalized();
						const Eigen::Vector3d v_ = (Eigen::Vector3d(G.row(T(f, k))) - Eigen::Vector3d(G.row(T(f, i)))).normalized();
						const Eigen::Vector3d tn = u_.cross(v_).normalized();
						if (n.dot(tn) < 0.2) return false;
					}
					return true;
				};

				const auto geoTest = [&filppedTest](
					int vi,
					int vj,
					const Eigen::MatrixXd& G,
					const Eigen::MatrixXi& T,
					const Eigen::VectorXi& fKeep,
					const std::vector<int>& common_faces,
					const std::vector<std::vector<int>>& VF,
					const std::vector<std::vector<int>>& VFi,
					Eigen::Vector3d &bP)->bool
				{
					std::set<int> cmf_set(common_faces.begin(), common_faces.end());
					double u;
					for (int i = 0; i < 10; i++)
					{
						u = static_cast<double>(i) / 2.;
						bP = u * Eigen::Vector3d(G.row(vi)) + (1. - u) * Eigen::Vector3d(G.row(vj));
						if(!filppedTest(bP, G, T, fKeep, cmf_set, VF[vi], VFi[vi]) ||
							!filppedTest(bP,G, T, fKeep, cmf_set, VF[vj], VFi[vj]))
						{
							continue;
						}
						return true;
					}
					return false; 
				};

				struct RelaxingNode
				{
					int vi{ -1 };
					int vj{ -1 };
					double sq_len{ std::numeric_limits<double>::max() };
				};
				struct RelaxingNodeCmp
				{
					bool operator()(const RelaxingNode& lhs, const RelaxingNode& rhs) { return lhs.sq_len > rhs.sq_len; }
				};
				using pQueue = std::priority_queue<RelaxingNode, std::vector<RelaxingNode>, RelaxingNodeCmp>;

				J = Eigen::VectorXi::LinSpaced(F.rows(), 0, F.rows() - 1);
				double avg_len = igl::avg_edge_length(V, F);
				double avg_sq_len = collapse_lenRatio * avg_len * avg_len;
				std::set<int> label_set(faces_label.begin(), faces_label.end());
				bool init_set = true;
				double collapse_ratio = 0.;
				int num_collapse = 0;
				int n_candidate = 1; 
				int stop_cond = 0;
				int iteration = 0; 
				while (stop_cond < 2)
				{
					pQueue pQ;
					for (auto label : label_set)
					{
						if (label == 0) continue;
						Eigen::MatrixXd Vi;
						Eigen::MatrixXi Fi;
						Eigen::VectorXi Ii, Ji, FiMAP;
						remove_unref(V, F, faces_label.unaryExpr([&label](auto fl) { return fl == label ? 1 : 0; }).eval(), Vi, Fi, Ii, Ji, FiMAP);

						std::vector<std::vector<int>> boundaries;
						igl::boundary_loop(Fi, boundaries);
						for (size_t i = 0; i < boundaries.size(); ++i)
						{
							auto& boundary = boundaries[i];
							std::for_each(boundary.begin(), boundary.end(), [&Ji](int& v) { v = Ji(v); });

							const size_t n = boundary.size();
							for (size_t j = 0; j < boundary.size(); ++j)
							{
								RelaxingNode node;
								node.vi = boundary[j];
								node.vj = boundary[(j + 1) % n];
								node.sq_len = (V.row(boundary[j]) - V.row(boundary[(j + 1) % n])).squaredNorm();
								pQ.push(node);
							}
						}
					}
					if (init_set)
					{
						init_set = false;
						n_candidate = static_cast<int>(pQ.size());
					}

					std::vector<std::vector<int>> VF, VFi;
					igl::vertex_triangle_adjacency(V, F, VF, VFi);
					Eigen::VectorXi vertices_label = Eigen::VectorXi::Zero(V.rows());
					for (int v = 0; v < V.rows(); ++v)
					{
						for (auto f : VF[v]) vertices_label(v) |= 1 << faces_label(f);
					}
					Eigen::VectorXi vmap = Eigen::VectorXi::LinSpaced(V.rows(), 0, V.rows() - 1);
					Eigen::VectorXi vmark = Eigen::VectorXi::Zero(V.rows());
					Eigen::VectorXi fKeep = Eigen::VectorXi::Ones(F.rows());
					collapse_ratio = static_cast<double>(num_collapse) / static_cast<double>(n_candidate);
					int local_num_collapse = 0;
					while (!pQ.empty())
					{
						const auto node = pQ.top();
						pQ.pop();
						auto vi = node.vi;
						auto vj = node.vj;
						if (vmap(vi) != vi || vmap(vj) != vj) continue;
						if (vmark(vi) || vmark(vj)) continue;
						if (isCommonVertex(vertices_label, vi) || isCommonVertex(vertices_label, vj)) continue;
						if (node.sq_len > avg_sq_len)
						{
							stop_cond = 1; 
							break;
						}

						if(!topoTest(F, VF, fKeep, vi, vj)) continue;

						std::vector<int> cmf = commonValues(VF[vi], VF[vj]);
						Eigen::Vector3d bP;
						if (!geoTest(vi, vj, V, F, fKeep, cmf, VF, VFi, bP)) continue;

						std::for_each(cmf.begin(), cmf.end(), [&fKeep](int f) { fKeep(f) = 0; });
						vmap[vj] = vi;
						V.row(vi) = bP;
						for(size_t i = 0; i < VF[vj].size(); ++i)
						{
							const int vf  = VF[vj][i];
							const int vfi = VFi[vj][i];
							while (vmap[F(vf, vfi)] != F(vf, vfi)) F(vf, vfi) = vmap[F(vf, vfi)];
						}
						VF[vi].insert(VF[vi].end(), VF[vj].begin(), VF[vj].end());
						VFi[vi].insert(VFi[vi].end(), VFi[vj].begin(), VFi[vj].end());
						vmark(vi) = 1;
						vmark(vj) = 1;
						local_num_collapse++;
						num_collapse++;
						collapse_ratio = static_cast<double>(num_collapse) / static_cast<double>(n_candidate);
						if(collapse_ratio >= 0.5)
						{
							stop_cond = 2;
							break;
						}
					}
					if (local_num_collapse == 0)
					{
						stop_cond = 3;
					}
					std::cout << "relaxing[iter(" << iteration++ << ") "
					<< "stop_cond(" << stop_cond << ") "
					<< "local_num_collapse(" << local_num_collapse << ") "
					<< "collapse_ratio(" << collapse_ratio << ").]\n";

					for (int& v : F.reshaped())
					{
						while (vmap[v] != v) v = vmap[v];
					}

					Eigen::MatrixXd Vr;
					Eigen::MatrixXi Fr;
					Eigen::VectorXi rI, rJ, rFMAP;
					remove_unref(V, F, fKeep, Vr, Fr, rI, rJ, rFMAP);
					std::vector<int> labels;
					for (int f = 0; f < faces_label.size(); ++f)
					{
						if (fKeep(f) > 0) labels.push_back(faces_label(f));
					}

					if (Fr.rows() != labels.size())
					{
						return;
					}

					V = Vr;
					F = Fr;
					faces_label = Eigen::Map<Eigen::VectorXi>(labels.data(), labels.size());
					J = rFMAP.unaryExpr([&J](int f) { return J(f); }).eval();
				}

				switch (stop_cond)
				{
				case 1:
					std::cout << "relaxing[stop condition: avg_sq_len(" << static_cast<double>(num_collapse) / static_cast<double>(n_candidate) << ").]" << std::endl;
					break;
				case 2:
					std::cout << "relaxing[stop condition: collapse_ratio(" << static_cast<double>(num_collapse) / static_cast<double>(n_candidate) << ").]" << std::endl;
					break;
				case 3:
					std::cout << "relaxing[stop condition: no edge can be collapsed(" << static_cast<double>(num_collapse) / static_cast<double>(n_candidate) << ").]" << std::endl;
					break;
				default:
					std::cout << "relaxing[stop condition: error(" << stop_cond << ").]" << std::endl;
					break;
				}
			}
			
			void remesh_along_isoline(
				Eigen::MatrixXd &V, 
				Eigen::MatrixXi &F, 
				Eigen::VectorXi &faces_label,
				Eigen::VectorXi& J,
				double phi)
			{
				V = V_;
				F = F_;
				faces_label = faces_label_;
				Eigen::SparseMatrix<double> BC0(V_.rows(), V_.rows()); BC0.setIdentity();
				J = Eigen::VectorXi::LinSpaced(F.rows(), 0, F.rows() - 1);
				for (int label = 0; label < data_ptrs_.size(); ++label)
				{
					if (data_ptrs_[label] == nullptr || !data_ptrs_[label]->PHI_ready) continue;
					auto& data_ptr = data_ptrs_[label];
					if (!data_ptr->isolines_test(data_ptr->PHI, phi)) continue;
					Eigen::VectorXd Phi = BC0 * data_ptr->get_PHI(V_.rows());

					Eigen::MatrixXd U;
					Eigen::MatrixXi G;
					Eigen::VectorXd SU;
					Eigen::VectorXi Jlocal;
					Eigen::SparseMatrix<double> BC;
					Eigen::VectorXi GI;
					igl::remesh_along_isoline(V, F, Phi, phi, U, G, SU, Jlocal, BC, GI);
					BC0 = BC * BC0;

					faces_label = Jlocal.unaryExpr([&faces_label](int j) { return faces_label(j); }).eval();
					J = Jlocal.unaryExpr([&J](int j) { return J(j); }).eval();
					V = U;
					F = G;
					auto lphi = data_ptr->average_value_of_labeled_faces(data_ptr->PHI);
					for(int f = 0; f < G.rows(); ++f)
					{
						double f_phi = 0.;
						for (int i = 0; i < 3; ++i)
						{
							f_phi += SU(F(f, i));
						}
						f_phi /= 3.;
						if(lphi * f_phi >= 0.)
						{
							faces_label(f) = label;
						}
						else
						{
							faces_label(f) = faces_label(f) == label ? 0 : faces_label(f);
						}
					}
				}
			}

			void mesh_connections(
				Eigen::MatrixXi& F,
				Eigen::VectorXi& MC)
			{
				const auto connectedComponent = [](
					const Eigen::SparseMatrix<int> &A,
					int seed,
					int cid,
					Eigen::VectorXi &FI)
				{
					std::queue<int> Q;
					Q.push(seed);
					while (!Q.empty())
					{
						const auto fid = Q.front();
						Q.pop();

						if (FI(fid) != -1)
						{
							continue;
						}
						FI(fid) = cid;

						for (Eigen::SparseMatrix<int>::InnerIterator it(A, fid); it; ++it)
						{
							const Eigen::SparseMatrix<int>::Index adj_fid = it.row();
							if (FI(adj_fid) != -1)
							{
								continue;
							}

							Q.push(adj_fid);
						}
					}
				};
				
				Eigen::SparseMatrix<int> A;
				igl::facet_adjacency_matrix(F, A);

				MC.setConstant(F.rows(), -1);
				for (int fid = 0, cid = 0; fid < F.rows(); ++fid)
				{
					if (MC(fid) != -1)
					{
						continue;
					}

					connectedComponent(A, fid, cid++, MC);
				}
			}

			void optimize_adjacency(
				Eigen::MatrixXd& V,
				Eigen::MatrixXi& F,
				Eigen::VectorXi& faces_label)
			{
				std::vector<std::vector<int> > VF, VFi;
				igl::vertex_triangle_adjacency(V, F, VF, VFi);

				Eigen::MatrixXd V0;
				Eigen::MatrixXi F0;
				Eigen::VectorXi I0, J0, F0MAP;
				remove_unref(V, F, faces_label.unaryExpr([](int fl) { return static_cast<int>(fl == 0); }).eval(), V0, F0, I0, J0, F0MAP);

				if (F0.rows() == 0)
				{
					return;
				}

				Eigen::VectorXi MC;
				mesh_connections(F0, MC);

				int nc = MC.maxCoeff() + 1;
				Eigen::VectorXi count = Eigen::VectorXi::Zero(nc);
				std::for_each(MC.data(), MC.data()+ MC.size(), [&count](int c){ count(c)++; });
				int c;
				count.maxCoeff(&c);
				count(c) = 0;
				const int nl = faces_label.maxCoeff() + 1;
				Eigen::VectorXi tmp_faces_label = faces_label;
				while(count.maxCoeff(&c) > 0)
				{
					Eigen::MatrixXd Vc;
					Eigen::MatrixXi Fc;
					Eigen::VectorXi Ic, Jc, FcMAP;
					remove_unref(V0, F0, MC.unaryExpr([&c](int mc) { return static_cast<int>(mc == c); }).eval(), Vc, Fc, Ic, Jc, FcMAP);

					std::vector<std::vector<int>> boundaries;
					igl::boundary_loop(Fc, boundaries);

					Eigen::VectorXi label_count = Eigen::VectorXi::Zero(nl);
					bool is_target = true; 
					for (const auto& boundary : boundaries)
					{
						for(auto v: boundary)
						{
							int bit = 0;
							for (int vf : VF[J0(Jc(v))])
							{
								bit |= 1 << faces_label(vf);
								label_count(faces_label(vf))++;
							}
							remove_bit(bit, 0);

							if (bit == 0)
							{
								is_target = false;
								break;
							}
						}

						if (!is_target)
						{
							break;
						}
					}

					label_count(0) = 0;

					if (is_target)
					{
						int label;
						if (label_count.maxCoeff(&label) > 0)
						{
							for(int f = 0; f < Fc.rows(); ++f)
							{
								tmp_faces_label(F0MAP(FcMAP(f))) = label;
							}
						}
					}
					count(c) = 0;
				}

				faces_label = tmp_faces_label;
			}

			void get_mesh(
				size_t label,
				Eigen::MatrixXd &V, 
				Eigen::MatrixXi &F)
			{
				if (label < data_ptrs_.size() && data_ptrs_[label] != nullptr)
				{
					V = data_ptrs_[label]->V; 
					F = data_ptrs_[label]->F; 
				}
			}

			void get_sdf(
				size_t label,
				Eigen::VectorXd &sdf)
			{
				if (label < data_ptrs_.size() && data_ptrs_[label] != nullptr && data_ptrs_[label]->SDF_ready)
				{
					sdf = data_ptrs_[label]->SDF;
				}
			}

			void get_flow(
				size_t label,
				Eigen::MatrixXd &flow)
			{
				if (label < data_ptrs_.size() && data_ptrs_[label] != nullptr && data_ptrs_[label]->X_ready)
				{
					flow = data_ptrs_[label]->X;
				}
			}

			void get_phi(
				size_t label,
				Eigen::VectorXd &phi)
			{
				if (label < data_ptrs_.size() && data_ptrs_[label] != nullptr && data_ptrs_[label]->PHI_ready)
				{
					phi = data_ptrs_[label]->PHI;
				}
			}

			void get_faces_label(
				size_t label,
				Eigen::VectorXi &faces_label)
			{
				if (label < data_ptrs_.size() && data_ptrs_[label] != nullptr)
				{
					faces_label = data_ptrs_[label]->faces_label;
				}
			}

			void update_vertices_label()
			{
				vertices_label_.resize(V_.rows(), Eigen::VectorX<bool>::Zero(faces_label_.maxCoeff() + 1));
				for(int v = 0; v < V_.rows(); ++v)
				{
					auto &vextex_label_ = vertices_label_[v];
					for (auto f : VF_[v]) 
					{
						vextex_label_(faces_label_(f)) = true; 
					}
				}
			}

			void load(
				const Eigen::MatrixXd& V,
				const Eigen::MatrixXi& F,
				const Eigen::VectorXi& faces_label)
			{
				V_ = V; 
				F_ = F; 
				faces_label_ = faces_label;
				is_border_vertex_ = igl::is_border_vertex(F_);
				igl::vertex_triangle_adjacency(V_, F_, VF_, VFi_);
				update_vertices_label();

				data_ptrs_.clear();
				data_ptrs_ = std::vector<std::unique_ptr<SmootherDATA>>(faces_label_.maxCoeff() + 1);
			}

			bool precompute(
				size_t label,
				const PrecomputeParam &param)
			{
				if(label >= data_ptrs_.size())
				{
					return false;  
				}
				
				if(data_ptrs_[label] != nullptr)
				{
					for(auto v: param.remove_vertices)
					{
						for(auto f: VF_[data_ptrs_[label]->VJ(v)])
						{			
							if(faces_label_(f) == label)
							{
								faces_label_(f) = 0;
							}
						}
					}
					update_vertices_label();
				}

				Eigen::MatrixXd V; 
				Eigen::MatrixXi F;
				Eigen::VectorXi faces_label, VJ;
				std::vector<Eigen::VectorX<bool>> vertices_label;
				std::vector<std::vector<int>> label_boundaries;
				if(!extracting(
					label, 
					param.fuzzy_kring,
					V,
					F, 
					label_boundaries,
					faces_label,
					vertices_label,
					VJ))
				{
					return false; 
				}

				data_ptrs_[label] = std::make_unique<SmootherDATA>();
				data_ptrs_[label]->V = std::move(V);
				data_ptrs_[label]->F = std::move(F);
				data_ptrs_[label]->faces_label = std::move(faces_label);
				data_ptrs_[label]->vertices_label = std::move(vertices_label);
				data_ptrs_[label]->VJ = std::move(VJ);
				std::cout << "current precompute label(" << label << ") #V(" << data_ptrs_[label]->V.rows() << ") #F(" <<  data_ptrs_[label]->F.rows() << ").\n";
				data_ptrs_[label]->precompute(param.use_intrinsic);
				data_ptrs_[label]->precompute(param.block_vertices, label_boundaries, is_border_vertex_, param.block_kring);
				if (data_ptrs_[label]->all_bEdges.empty() || data_ptrs_[label]->all_bVertices.empty()) 
				{
					std::cout << "[SmootherImpl::precompute unexpected error: extracting label " << label << "  failed]\n";
					data_ptrs_[label].reset(nullptr);
				} 
				return true; 
			}

			bool sdf_solve(
				size_t label,
				const SDFParam &param)
			{
				if(label >= data_ptrs_.size() || data_ptrs_[label] == nullptr)
				{
					return false;  
				}

				data_ptrs_[label]->solve_sdf(param);
				return data_ptrs_[label]->SDF_ready;
			}

			bool vector_field_solve(
				size_t label,
				const VectorTransportParam& param)
			{
				if(label >= data_ptrs_.size() || data_ptrs_[label] == nullptr)
				{
					return false;  
				}

				data_ptrs_[label]->solve_X(param);
				return data_ptrs_[label]->X_ready;
			}

			bool scalar_field_solve(
				size_t label,
				const DiffusionParam& param)
			{
				if(label >= data_ptrs_.size() || data_ptrs_[label] == nullptr)
				{
					return false;  
				}

				data_ptrs_[label]->solve_PHI(param);
				return data_ptrs_[label]->PHI_ready; 
			}

		private:
			Eigen::MatrixXd V_;
			Eigen::MatrixXi F_;
			Eigen::VectorXi faces_label_;
			std::vector<Eigen::VectorX<bool>> vertices_label_;
			std::vector<bool> is_border_vertex_;
			std::vector<std::vector<int>> VF_, VFi_;
			std::vector<std::unique_ptr<SmootherDATA>> data_ptrs_;
		};

		Smoother::Smoother() {}

		Smoother::~Smoother() {}

		void Smoother::precompute(
			const std::vector<double>& vbuffer,
			const std::vector<int>& fbuffer,
			const std::vector<int>& faces_label,
			PrecomputeParam& param)
		{
			assert(vbuffer.size() % 3 == 0);
			assert(fbuffer.size() % 3 == 0);

			impl_ = std::make_unique<SmootherImpl>();
			impl_->precompute(vbuffer, fbuffer);
			impl_->precompute(faces_label, param);
		}

		void Smoother::precompute(
			const Eigen::MatrixXd &V, 
			const Eigen::MatrixXi &F, 
			const Eigen::VectorXi &faces_label,
			PrecomputeParam &param)
		{

			std::vector<double> vbuffer(V.size());
			Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(vbuffer.data(), V.rows(), V.cols()) = V;
			std::vector<int> fbuffer(F.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(fbuffer.data(), F.rows(), F.cols()) = F;
			std::vector<int> labels(faces_label.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(labels.data(), faces_label.rows(), faces_label.cols()) = faces_label;
		
			assert(vbuffer.size() % 3 == 0);
			assert(fbuffer.size() % 3 == 0);

			impl_ = std::make_unique<SmootherImpl>();
			impl_->precompute(vbuffer, fbuffer);
			impl_->precompute(labels, param);
		}

		void Smoother::sdf_required(const SDFParam& param)
		{
			impl_->sdf_required(param);
		}

		void Smoother::vector_field_solve(const VectorTransportParam& param)
		{
			impl_->vector_field_solve(param);
		}

		void Smoother::scalar_field_solve(const DiffusionParam& param)
		{
			impl_->scalar_field_solve(param);
		}

		void Smoother::remesh_along_isoline(
			double phi,
			std::vector<double>& vbuffer,
			std::vector<int>& fbuffer,
			std::vector<int>& labels,
			std::vector<int>& fmap)
		{
			Eigen::MatrixXd V;
			Eigen::MatrixXi F;
			Eigen::VectorXi faces_label;
			Eigen::VectorXi J;
			impl_->remesh_along_isoline(V, F, faces_label, J, phi);

			vbuffer.resize(V.size());
			Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(vbuffer.data(), V.rows(), V.cols()) = V;
			fbuffer.resize(F.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(fbuffer.data(), F.rows(), F.cols()) = F;
			labels.resize(faces_label.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(labels.data(), faces_label.rows(), faces_label.cols()) = faces_label;
			fmap.resize(J.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(fmap.data(), J.rows(), J.cols()) = J;
		}

		void Smoother::relaxing(
			double collapse_lenRatio,
			std::vector<double>& vbuffer,
			std::vector<int>& fbuffer,
			std::vector<int>& labels,
			std::vector<int>& fmap)
		{
			Eigen::MatrixXd V = Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(vbuffer.data(), vbuffer.size() / 3, 3); 
			Eigen::MatrixXi F = Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(fbuffer.data(), fbuffer.size() / 3, 3);
			Eigen::VectorXi faces_label = Eigen::Map<Eigen::VectorXi>(labels.data(), labels.size());
			Eigen::VectorXi J = Eigen::Map<Eigen::VectorXi>(fmap.data(), fmap.size());

			impl_->relaxing(collapse_lenRatio, V, F, faces_label, J);

			vbuffer.resize(V.size());
			Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(vbuffer.data(), V.rows(), V.cols()) = V;
			fbuffer.resize(F.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(fbuffer.data(), F.rows(), F.cols()) = F;
			labels.resize(faces_label.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(labels.data(), faces_label.rows(), faces_label.cols()) = faces_label;
			fmap.resize(J.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(fmap.data(), J.rows(), J.cols()) = J;
		}

		void Smoother::optimize_adjacency(
			std::vector<double>& vbuffer,
			std::vector<int>& fbuffer,
			std::vector<int>& labels)
		{
			Eigen::MatrixXd V = Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(vbuffer.data(), vbuffer.size() / 3, 3);
			Eigen::MatrixXi F = Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(fbuffer.data(), fbuffer.size() / 3, 3);
			Eigen::VectorXi faces_label = Eigen::Map<Eigen::VectorXi>(labels.data(), labels.size());

			impl_->optimize_adjacency(V, F, faces_label);

			vbuffer.resize(V.size());
			Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(vbuffer.data(), V.rows(), V.cols()) = V;
			fbuffer.resize(F.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(fbuffer.data(), F.rows(), F.cols()) = F;
			labels.resize(faces_label.size());
			Eigen::Map<Eigen::Matrix<int, -1, -1, Eigen::RowMajor>>(labels.data(), faces_label.rows(), faces_label.cols()) = faces_label;
		}

		void Smoother::get_mesh(
			size_t label,
			Eigen::MatrixXd &V, 
			Eigen::MatrixXi &F)
		{
			impl_->get_mesh(label, V, F);
		}

		void Smoother::get_sdf(
			size_t label,
			Eigen::VectorXd &sdf)
		{
			impl_->get_sdf(label, sdf);
		}

		void Smoother::get_flow(
			size_t label,
			Eigen::MatrixXd &flow)
		{
			impl_->get_flow(label, flow);
		}

		void Smoother::get_phi(
			size_t label,
			Eigen::VectorXd &phi)
		{
			impl_->get_phi(label, phi);
		}

		void Smoother::get_faces_label(
			size_t label,
			Eigen::VectorXi &faces_label)
		{
			impl_->get_faces_label(label, faces_label);
		}

		void Smoother::load(
			const Eigen::MatrixXd &V, 
			const Eigen::MatrixXi &F, 
			const Eigen::VectorXi &faces_label)
		{
			impl_ = std::make_unique<SmootherImpl>();
			impl_->load(V, F, faces_label);
		}

		bool Smoother::precompute(
			size_t label,
			const PrecomputeParam &param)
		{
			return impl_->precompute(label, param);
		}

		bool Smoother::sdf_solve(
			size_t label,
			const SDFParam &param)
		{
			return impl_->sdf_solve(label, param);
		}
		
		bool Smoother::vector_field_solve(
			size_t label,
			const VectorTransportParam& param)
		{
			return impl_->vector_field_solve(label, param);
		}

		bool Smoother::scalar_field_solve(
			size_t label,
			const DiffusionParam& param)
		{
			return impl_->scalar_field_solve(label, param);
		}
	}
}
