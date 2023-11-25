#ifndef LIBDGP_LSCM_H
#define LIBDGP_LSCM_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    DGP_INLINE void lscm(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::MatrixXd &U);
}

#ifndef DGP_STATIC_LIBRARY
#include   "lscm.cpp"
#endif

#endif //LIBDGP_LSCM_H