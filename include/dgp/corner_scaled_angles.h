#ifndef LIBDGP_CORNER_SCALED_ANGLES_H
#define LIBDGP_CORNER_SCALED_ANGLES_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    DGP_INLINE void CornerScaledAngles(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        std::vector<bool> &is_border_vertex,
        Eigen::MatrixXd &corner_angles);
}

#ifndef DGP_STATIC_LIBRARY
#include   "corner_scaled_angles.cpp"
#endif

#endif //LIBDGP_CORNER_SCALED_ANGLES_H