#ifndef LIBDGP_INTERSECTION2D_H
#define LIBDGP_INTERSECTION2D_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    DGP_INLINE bool intersection2d(
        const Eigen::Vector3d& O1, 
        const Eigen::Vector3d& E1,
        const Eigen::Vector3d& O2,
        const Eigen::Vector3d& E2,
        Eigen::Vector3d& pt, 
        double& t1, 
        double& t2);
}

#ifndef DGP_STATIC_LIBRARY
#include   "intersection2d.cpp"
#endif

#endif //LIBDGP_INTERSECTION2D_H