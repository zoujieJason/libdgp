#ifndef LIBDGP_REMOVE_UNREF_H
#define LIBDGP_REMOVE_UNREF_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    void remove_unref(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXi& K,
        Eigen::MatrixXd& Vo,
        Eigen::MatrixXi& Fo,
        Eigen::VectorXi& Io,
        Eigen::VectorXi& Jo,
        Eigen::VectorXi& FoMAP);
}

#ifndef DGP_STATIC_LIBRARY
#include   "remove_unref.cpp"
#endif

#endif //LIBDGP_REMOVE_UNREF_H