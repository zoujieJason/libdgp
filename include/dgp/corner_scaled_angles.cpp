#include "corner_scaled_angles.h"

#include <igl/internal_angles.h>
#include <igl/is_border_vertex.h>

#include <vector>

namespace dgp
{
    DGP_INLINE void CornerScaledAngles(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        std::vector<bool> &is_border_vertex,
        Eigen::MatrixXd &corner_angles)
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

        is_border_vertex = igl::is_border_vertex(F);
        for(int f = 0; f < F.rows(); f++)
        {
            for(int i = 0; i < 3; i++)
            {
                double scale = (is_border_vertex[F(f, i)]? 1.0 : 2.0) * 180. / vectices_angle(F(f, i));
                corner_angles(f, i) *= scale;
            }
        }
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif
