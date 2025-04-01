#include "vector_transport.h"
#include "halfedge_opposite.h"
#include "vertex_halfedge_adjacency.h"

#include <igl/internal_angles.h>
#include <igl/is_border_vertex.h>

#include <vector>

namespace dgp
{
    bool vector_transport_precompute(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &E,
        const std::vector<std::vector<int> > &VF,
        const std::vector<std::vector<int> > &VFi, 
        const std::vector<std::vector<int> > &EF,
        const std::vector<std::vector<int> > &EFi, 
        const Eigen::VectorXd &l,
        const Eigen::MatrixXd &C)
    {
        Eigen::MatrixXd corner_angles; 
        igl::internal_angles(V, F, corner_angles);
        Eigen::VectorXd vectices_angle = Eigen::VectorXd::Zero(V.rows());
        for(int f = 0; f < F.rows(); f++)
        {
            for(int i = 0; i < 3; i++)
            {
                vectices_angle(F(f, i)) += corner_angles(f, i);
            }
        }

        const std::vector<bool> is_border_vertex = igl::is_border_vertex(F);
        for(int f = 0; f < F.rows(); f++)
        {
            for(int i = 0; i < 3; i++)
            {
                double scale = (is_border_vertex[F(f, i)]? 1.0 : 2.0) * 180. / vectices_angle(F(f, i));
                corner_angles(f, i) *= scale;
            }
        }

        Eigen::MatrixXi HeOpp;
        halfedge_opposite(E, HeOpp);

        std::vector<std::vector<std::pair<int, int>>> VHe;
        if(!vertex_halfedge_adjacency(VF, VFi, F, HeOpp, is_border_vertex, VHe)) return false;
        
        Eigen::Matrix<std::complex<double>, -1, -1> HeCplx, HeOppCplx;
        HeCplx.setZero(F.rows(), 3);
        HeOppCplx.setZero(F.rows(), 3);
        for(size_t v = 0; v < VHe.size(); ++v)
        {
            const auto &helist = VHe[v];
            double angle = 0.; 
            for(const std::pair<int,int> &he : helist)
            {
                int f = he.first;
                int k = he.second;   // k -> e(i, j);
                int i = (k + 1) % 3; // i -> e(j, k);
                std::complex<double> cpl = l(f, k) * std::complex(std::cos(angle), std::sin(angle));
                HeCplx(f, k) = cpl;
                auto he_opp = HeOpp(f, k);
                if(he_opp != -1)
                {
                    int f_opp = he_opp % F.rows();
                    int j_opp = he_opp / F.rows();
                    HeOppCplx(f, j_opp) = cpl;
                }
                angle += corner_angles(f, i);
            }
            if(!is_border_vertex[v]) continue;
            auto he_back = helist.back();
            int f = he_back.first;
            int k = he_back.second;
            int j = (k + 2) % 3;
            auto he_opp = HeOpp(f, j);
            assert(he_opp == -1);
            std::complex<double> cpl = l(f, j) * std::complex(std::cos(angle), std::sin(angle));
            HeOppCplx(f, j) = cpl;
        }

        std::vector<Eigen::Triplet<std::complex<double>>> triplets;
        for(size_t e = 0; e < EF.size(); ++e)
        {
            const auto &ef = EF[e];
            const auto &efi = EFi[e];
            assert(ef.size() > 0 && ef.size() < 3); 
            int f0 = ef.front(); 
            int i0 = efi.front();
            std::complex<double> cpl = -HeCplx(f0, i0) / HeOppCplx(f0, i0);
            cpl /= std::abs(cpl);
            std::complex<double> cpl_inv = std::conj(cpl) / std::norm(cpl);

            int j0 = (i0 + 1) % 3;
            int k0 = (i0 + 2) % 3;
            int vj = F(f0, j0);
            int vk = F(f0, k0);

            // he -> (vj, vk)
            double c = C(f0, i0);
            triplets.emplace_back(vj, vj, c);
            triplets.emplace_back(vj, vk, -c * cpl_inv);

            triplets.emplace_back(vk, vk, c);
            triplets.emplace_back(vk, vj, -c * cpl);
        }
        
        using SpMat = Eigen::SparseMatrix<std::complex<double>>;
        SpMat Lconn = SpMat(V.rows(), V.rows());
        Lconn.setFromTriplets(triplets.begin(), triplets.end());
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif
