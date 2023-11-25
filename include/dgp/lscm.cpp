#include "lscm.h"
#include "local_isometric_parameterization.h"

#include <igl/boundary_loop.h>
#include <igl/doublearea.h>

#include <igl/lscm.h>

namespace dgp
{
    DGP_INLINE void lscm(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::MatrixXd &U)
    {
        Eigen::MatrixXd LV; 
        local_isometric_parameterization(V, F, LV);
        
        Eigen::VectorXd dlbA;
        igl::doublearea(V, F, dlbA);

        const auto nv = V.rows();
        const auto nf = F.rows();
        const auto niv = nv - 2;

        std::vector<size_t> mesh_boundary;
        igl::boundary_loop(F, mesh_boundary);

        const auto v0 = std::min(mesh_boundary[0], mesh_boundary[1]); 
        const auto v1 = std::max(mesh_boundary[0], mesh_boundary[1]);
        const double scale = (V.row(v0)-V.row(v1)).norm();
        Eigen::VectorXd UV;
        UV.setZero(2 * nv);
        UV(0) = 0.0;
        UV(2) = 0.0;
        UV(1) = scale;
        UV(3) = 0.0;
        
        std::vector<Eigen::Triplet<double>> triplets;
        for(int fid = 0; fid < nf; ++fid)
        {
            const auto &f = F.row(fid);
            
            const auto &lv1 = LV.row(fid * 3 + 0);
            const auto &lv2 = LV.row(fid * 3 + 1);
            const auto &lv3 = LV.row(fid * 3 + 2);
            const double area = 0.5 * dlbA(fid);

            Eigen::MatrixXd coeffs;
            coeffs.resize(3,2);
            coeffs.row(0) = Eigen::Vector2d(lv3(0) - lv2(0), lv3(1) - lv2(1));
            coeffs.row(1) = Eigen::Vector2d(lv1(0) - lv3(0) ,lv1(1) - lv3(1));
            coeffs.row(2) = Eigen::Vector2d(lv2(0) - lv1(0), lv2(1) - lv1(1));

            for(int i = 0; i < 3; ++i)
            {
                const auto& vid = f(i);
                auto curr_coeff = coeffs.row(i);
                if (vid == v0) 
                {
                    curr_coeff = -1.0 * curr_coeff;
                    //real part
                    triplets.emplace_back(fid, 0, curr_coeff(0) / area);
                    triplets.emplace_back(fid + nf, 2, curr_coeff(0) / area);

                    //image part
                    triplets.emplace_back(fid, 2, -curr_coeff(1) / area);
                    triplets.emplace_back(fid + nf, 0, curr_coeff(1) / area);
                } 
                else if (vid == v1) 
                {
                    curr_coeff = -1.0 * curr_coeff;
                    //real part
                    triplets.emplace_back(fid, 1, curr_coeff(0) / area);
                    triplets.emplace_back(fid + nf, 3, curr_coeff(0) / area);

                    //image part
                    triplets.emplace_back(fid, 3, -curr_coeff(1) / area);
                    triplets.emplace_back(fid + nf, 1, curr_coeff(1) / area);
                } 
                else 
                {
                    const int offset = vid < v0 ? 0 : (vid < v1 ? 1 : 2); 
                    //real part
                    triplets.emplace_back(fid, vid + 4 - offset, curr_coeff(0) / area);
                    triplets.emplace_back(fid + nf, niv + vid + 4 - offset, curr_coeff(0) / area);

                    //image part
                    triplets.emplace_back(fid, niv + vid + 4 - offset, -curr_coeff(1) / area);
                    triplets.emplace_back(fid + nf, vid + 4 - offset, curr_coeff(1) / area);
                }
            }
        }

        Eigen::SparseMatrix<double> M(2*nf,2*nv);
        M.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver(M.block(0,4,2*nf,2*(nv-2)));
        const Eigen::VectorXd b = M.block(0,0,2*nf,4) * UV.block(0,0,4,1);
        UV.block(4,0,2*(nv-2),1) = solver.solve(b);

        U.setZero(nv, 2);
        U.row(v0) = Eigen::Vector2d(0.0,0.0);
        U.row(v1) = Eigen::Vector2d(scale,0.0);
        for(int vid = 0; vid < nv; ++vid) 
        {
            if(vid == v0 || vid == v1)
            {
                continue;
            }
            const int offset = vid < v0 ? 0 : (vid < v1 ? 1 : 2); 
            U.row(vid) = Eigen::Vector2d(UV(vid+4-offset), UV(vid+4+nv-2-offset));
        }
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif