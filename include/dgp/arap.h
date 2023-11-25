#ifndef LIBDGP_ARAP_H
#define LIBDGP_ARAP_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    DGP_INLINE void arap(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXd &U0,
        Eigen::MatrixXd &U,
        int max_iteration = 20,
        double eps = 1e-6);
}

#ifndef DGP_STATIC_LIBRARY
#include   "arap.cpp"
#endif

#endif //LIBDGP_ARAP_H