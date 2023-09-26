#ifndef LIBDGP_UNIT_CIRCLE_POINTS_H
#define LIBDGP_UNIT_CIRCLE_POINTS_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    DGP_INLINE void unit_circle_points(size_t n, double radius, Eigen::MatrixXd &P);
}

#ifndef DGP_STATIC_LIBRARY
#include   "unit_circle_points.cpp"
#endif

#endif //LIBDGP_UNIT_CIRCLE_POINTS_H