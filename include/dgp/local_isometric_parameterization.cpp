#include "local_isometric_parameterization.h"

namespace dgp
{
    DGP_INLINE void local_isometric_parameterization(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::MatrixXd &U)
    {
        assert(F.cols() == 3);
        const auto nf = F.rows();
        
        U.resize(nf * 3, 2);
        for(int fid = 0; fid < nf; ++fid)
        {
            const auto &f = F.row(fid);
            const Eigen::Vector3d p1 = V.row(f(0));
            const Eigen::Vector3d p2 = V.row(f(1));
            const Eigen::Vector3d p3 = V.row(f(2));

            const Eigen::Vector3d e1 = p2 - p1;
            const Eigen::Vector3d e2 = p3 - p1;
            
            const Eigen::Vector3d n = e1.cross(e2);
            const Eigen::Vector3d y = n.cross(e1).normalized();
            const Eigen::Vector3d x = e1.normalized();

            U.row(fid * 3 + 0) = Eigen::Vector2d(0.0, 0.0);
            U.row(fid * 3 + 1) = Eigen::Vector2d(e1.norm(), 0.0);
            U.row(fid * 3 + 2) = Eigen::Vector2d(e2.dot(x), e2.dot(y));
        }
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif
