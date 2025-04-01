#ifndef LIBDGP_VECTOR_TRANSPORT_H
#define LIBDGP_VECTOR_TRANSPORT_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    DGP_INLINE void CornerScaledAngles(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const std::vector<bool> &is_border_vertex,
        Eigen::MatrixXd &corner_angles);

    DGP_INLINE void vector_transport_precompute(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F);
}

#ifndef DGP_STATIC_LIBRARY
#include   "vector_transport.cpp"
#endif

#endif //LIBDGP_VECTOR_TRANSPORT_H