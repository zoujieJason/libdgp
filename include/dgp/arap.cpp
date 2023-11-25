#include "arap.h"
#include "local_isometric_parameterization.h"

#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/polar_svd.h>
#include <igl/vertex_triangle_adjacency.h>

#include <numeric>

#define LIBDGP_ARAP_DEBUG

namespace dgp
{
    DGP_INLINE void arap(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXd &U0,
        Eigen::MatrixXd &U,
        int max_iteration,
        double eps)
    {
        assert(U0.rows() == V.rows());
        assert(F.cols() == 3);

        const auto nv = V.rows();
        const auto nf = F.rows();

        Eigen::MatrixXd LV;
        local_isometric_parameterization(V, F, LV);

        Eigen::MatrixXd C;
        igl::cotmatrix_entries(V, F, C);

        const auto local_phase = [&C, &F, &LV, &nf](const Eigen::MatrixXd &UV, Eigen::MatrixXd &Rs) 
        {
            Rs.setZero(nf * 2, 2);
            for(int fid = 0; fid < nf; ++fid)
            {
                Eigen::Matrix2d si = Eigen::Matrix2d::Zero();
                for (size_t i = 0; i < 3; ++i)
                {
                    const double cotan = 2.0 * C(fid, (i + 2) % 3);
                    const Eigen::Vector2d u = UV.row(F(fid, i))   - UV.row(F(fid, (i + 1) % 3));
                    const Eigen::Vector2d x = LV.row(3 * fid + i) - LV.row(3 * fid + (i + 1) % 3);
                    si = si + cotan * u * x.transpose();
                }

				Eigen::Matrix2d ri, ti, ui, vi;
				Eigen::Vector2d _;
				igl::polar_svd(si, ri, ti, ui, _, vi);

                Rs.block(fid * 2, 0, 2, 2) = ri;
            }
        };

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
		{
            Eigen::SparseMatrix<double> L;
            igl::cotmatrix(V, F, L);
            solver.compute(-L);
            if (solver.info() != Eigen::Success) 
            {
                std::cerr << "matrix decomposition failed." << std::endl;
                return;
            }
        }
    
        std::vector<std::vector<int> > VF, VFI; 
		igl::vertex_triangle_adjacency(V, F, VF, VFI); 

        const auto global_phase = [&C, &F, &solver, &VF, &VFI, &LV, &nv](const Eigen::MatrixXd &Rs, Eigen::MatrixXd &B)
        {
            B.setZero(nv, 2);
            for(int vid = 0; vid < nv; ++vid)
            {
                Eigen::Vector2d b = Eigen::Vector2d::Zero();
                const auto &vf = VF[vid];
                const auto &vfi = VFI[vid];
                assert(vf.size() == vfi.size());
                for(int i = 0; i < vf.size(); ++i)
                {
                    const auto &fid = vf[i];
                    const auto &f = F.row(fid);

                    const auto &i1 = vfi[i];
                    const auto i2 = (i1+1)%3;
                    const auto i3 = (i2+1)%3;

                    const Eigen::Matrix2d R = Rs.block(fid * 2, 0, 2, 2);

                    const Eigen::Vector2d x_ij = LV.row(3 * fid + i1) - LV.row(3 * fid + i2);
                    b += C(fid, i3) * R * x_ij;

                    const Eigen::Vector2d x_ik = LV.row(3 * fid + i1) - LV.row(3 * fid + i3);
                    b += C(fid, i2) * R * x_ik;
                }
                B.row(vid) = b;
            }
        };

        const auto arap_energy = [&C, &nf, &LV, &F](const Eigen::MatrixXd &UV, const Eigen::MatrixXd &Rs)->double
        {
            double e = 0.0; 
            for(int fid = 0; fid < nf; ++fid)
            {
                const auto &f = F.row(fid);
                const Eigen::Matrix2d R = Rs.block(fid * 2, 0, 2, 2);
                for (size_t i = 0; i < 3; ++i) 
                {
					const Eigen::Vector2d u = UV.row(f(i)) - UV.row(f((i + 1) % 3));
					const Eigen::Vector2d x = LV.row(3 * fid + i) - LV.row(3 * fid + (i + 1) % 3);
					e += C(fid, (i + 2) % 3) * (u - R * x).norm(); 
				}
            }
            return e;
        };

        U = U0;
        double eprev = std::numeric_limits<double>::max(); 
        for(int i = 0; i < max_iteration; ++i)
        {
            Eigen::MatrixXd Rs; 
            local_phase(U, Rs);
            const auto ecurr = arap_energy(U, Rs);
#ifdef LIBDGP_ARAP_DEBUG
            std::cout << "Iteration " << i << ": eprev " << eprev << " " << "| ecurr " << ecurr << " | " <<  eprev - ecurr << std::endl;
#endif
            if(eprev - ecurr < eps)
            {
                std::cout << "Converged." << std::endl;
                return;
            }
            eprev = ecurr;

            Eigen::MatrixXd B;
            global_phase(Rs, B);

            U = solver.solve(B);
			if (solver.info() != Eigen::Success) 
            {
				std::cerr << "solving failed." << std::endl;
				return;
			}
        }
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif
