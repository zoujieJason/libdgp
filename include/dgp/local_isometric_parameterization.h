#ifndef LIBDGP_LOCAL_ISOMETRIC_PARAMETERIZATION_H
#define LIBDGP_LOCAL_ISOMETRIC_PARAMETERIZATION_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    DGP_INLINE void local_isometric_parameterization(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::MatrixXd &U);
}

#ifndef DGP_STATIC_LIBRARY
#include   "local_isometric_parameterization.cpp"
#endif

#endif //LIBDGP_LOCAL_ISOMETRIC_PARAMETERIZATION_H