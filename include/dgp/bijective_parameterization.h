#ifndef LIBDGP_BIJECTIVE_PARAMETERIZATION_H
#define LIBDGP_BIJECTIVE_PARAMETERIZATION_H

#include "dgp_inline.h"

#include <Eigen/Dense>

namespace dgp
{
    class BijectiveParameterization
    {
public:
        BijectiveParameterization(){};
        
        DGP_INLINE void initialize(
            const Eigen::MatrixXd &V,
            const Eigen::MatrixXi &F,
            const Eigen::MatrixXd &U);

        DGP_INLINE void per_iteration(int iteration, double &g_norm, double &lambda);

        DGP_INLINE void output(Eigen::MatrixXd &U) const;

private:
        DGP_INLINE double energy_dirichlet_part1(const Eigen::VectorXd &dlbA_U, int fid) const;

        DGP_INLINE double energy_dirichlet_part2(const Eigen::MatrixXd &U, int fid) const;

        DGP_INLINE double energy_dirichlet(const Eigen::MatrixXd &U) const;

        DGP_INLINE void calculate_boundary_gradient(const Eigen::MatrixXd &U, Eigen::VectorXd &g) const;

        DGP_INLINE void calculate_gradient(const Eigen::MatrixXd &U, Eigen::VectorXd &g) const;

        DGP_INLINE bool non_degenerate(const Eigen::MatrixXd &U) const; 

        DGP_INLINE double maximum_step_length(const Eigen::MatrixXd &U, const Eigen::MatrixXd &Delta) const;

        DGP_INLINE double calculate_step_length(const Eigen::MatrixXd &U, const Eigen::MatrixXd &Delta, double max_lambda) const;

private:
        Eigen::MatrixXd V_;
        Eigen::MatrixXi F_;
        Eigen::MatrixXd C_;
        Eigen::VectorXd U_;

        Eigen::VectorXd dlbA_V_; 

        double plus_zero_ {0.0};

        int nv_;
        int nf_;

        // L_BFGS variables
        int number_of_LBFGS_variables_ {7};
        Eigen::VectorXd alpha_, rho_;
        Eigen::MatrixXd S_k_, Y_k_; 
    };
}

#ifndef DGP_STATIC_LIBRARY
#include   "bijective_parameterization.cpp"
#endif

#endif //LIBDGP_BIJECTIVE_PARAMETERIZATION_H