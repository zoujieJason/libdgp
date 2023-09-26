#include "unit_circle_points.h"

namespace dgp
{
    DGP_INLINE void unit_circle_points(size_t n, double radius, Eigen::MatrixXd &P)
    {
        P.resize(n, 2);
        const auto per_radian = 2.0 * M_PI / static_cast<double>(n);
        for (size_t i = 0; i < n; ++i)
        {
            P(i, 0) = radius * std::cos(static_cast<double>(i) * per_radian);
            P(i, 1) = radius * std::sin(static_cast<double>(i) * per_radian);
        }
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif