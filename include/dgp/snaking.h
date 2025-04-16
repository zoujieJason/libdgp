#ifndef LIBDGP_SNAKING_H
#define LIBDGP_SNAKING_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include "Snaxl.h"

#include <list>

namespace dgp
{
    struct SnakingData
    {
        Eigen::MatrixXi E, uE; 
        Eigen::VectorXi EMAP; 
        std::vector<std::vector<int> > uE2E; 
        std::vector<std::vector<int> > uE2Fsorted;

        std::vector<std::vector<int> > VF, VFi;
        Eigen::MatrixXi HeOpp;
        std::vector<std::pair<bool, std::vector<std::tuple<int, int, bool>>>> VHe;

        double threshold_t{1e-3};
    };

    DGP_INLINE bool snaking_precompute(
        const Eigen::MatrixXi &F,
        SnakingData &snaking_data);

    template <typename stdit>
    DGP_INLINE void extract_snaxls(
        const Eigen::MatrixXd &V, 
        stdit begin, 
        stdit end,
        Eigen::MatrixXd &iV,
        Eigen::MatrixXi &iE,
        bool loop = false);

    DGP_INLINE void handle_critical_vertices(
        std::list<Snaxl> &snaxls);

    DGP_INLINE bool snaking(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const SnakingData &snaking_data,
        std::list<Snaxl> &snaxls);
}

#ifndef DGP_STATIC_LIBRARY
#include   "snaking.cpp"
#endif

#endif //LIBDGP_SNAKING_H