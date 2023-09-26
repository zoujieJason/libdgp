#include "bijective_parameterization.h"

#include <igl/cotmatrix_entries.h>
#include <igl/doublearea.h>

namespace dgp
{
    DGP_INLINE void BijectiveParameterization::initialize(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXd &U)
    {
        nv_ = V.rows();
        nf_ = F.rows();
        V_ = V;
        F_ = F;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> U_R(U);
        U_ = Eigen::Map<Eigen::VectorXd>(U_R.data(), U_R.size());
        igl::cotmatrix_entries(V, F, C_);
        igl::doublearea(V, F, dlbA_V_);

        alpha_.setZero(number_of_LBFGS_variables_); 
        rho_.setZero(number_of_LBFGS_variables_);
        S_k_.setZero(number_of_LBFGS_variables_, U_.size());
        Y_k_.setZero(number_of_LBFGS_variables_, U_.size());
    }

    DGP_INLINE double BijectiveParameterization::energy_dirichlet_part1(const Eigen::VectorXd &dlbA_U, int fid) const 
    {
        const auto sq_area_p = 0.25 * dlbA_V_(fid) * dlbA_V_(fid);
        const auto sq_area_u = 0.25 * dlbA_U(fid) * dlbA_U(fid);
        return 1.0 + sq_area_p / (sq_area_u + plus_zero_);
    }

    DGP_INLINE double BijectiveParameterization::energy_dirichlet_part2(const Eigen::MatrixXd &U, int fid) const 
    {
        const auto &f = F_.row(fid);
        const auto &u1 = U.row(f(0));
        const auto &u2 = U.row(f(1));
        const auto &u3 = U.row(f(2));

        const auto &p1 = V_.row(f(0));
        const auto &p2 = V_.row(f(1));
        const auto &p3 = V_.row(f(2));

        return (((u3 - u1).squaredNorm() * (p2 - p1).squaredNorm() + (u2 - u1).squaredNorm() * (p3 - p1).squaredNorm()) - 2.0 * ((u3 - u1).dot(u2 - u1) * (p3 - p1).dot(p2 - p1))) / (2.0 * dlbA_V_(fid) + plus_zero_);
    }

    DGP_INLINE double BijectiveParameterization::energy_dirichlet(const Eigen::MatrixXd &U) const
    {
        Eigen::VectorXd dlbA_U;
        igl::doublearea(U, F_, dlbA_U);
        double e = 0.0;
        for(int fid = 0; fid < nf_; ++fid)
        {
            e += energy_dirichlet_part1(dlbA_U, fid) * energy_dirichlet_part2(U, fid);
        }
        return e;
    }

    DGP_INLINE void BijectiveParameterization::calculate_boundary_gradient(const Eigen::MatrixXd &U, Eigen::VectorXd &g) const
    {
        
    }
    
    DGP_INLINE void BijectiveParameterization::calculate_gradient(const Eigen::MatrixXd &U, Eigen::VectorXd &g) const
    {
        Eigen::VectorXd dlbA_U;
        igl::doublearea(U, F_, dlbA_U);
        Eigen::Matrix2d R; 
        R << 0.0, -1.0, 1.0, 0.0;
        g.setZero(nv_ * 2);
        for(int fid = 0; fid < nf_; ++fid)
        {
            const auto &f = F_.row(fid);
            const auto _1 = energy_dirichlet_part1(dlbA_U, fid);
            const auto _2 = energy_dirichlet_part2(U, fid);
            const auto const_term = - (0.25 * dlbA_V_(fid) * dlbA_V_(fid)) / (0.125 * dlbA_U(fid) * dlbA_U(fid) * dlbA_U(fid) + plus_zero_); 
            for(int i = 0; i < 3; ++i)
            {
                const auto v1 = f(i);
                const auto v2 = f((i + 1) % 3);
                const auto v3 = f((i + 2) % 3);

                const Eigen::Vector2d ui1 = U.row(v1);
                const Eigen::Vector2d ui2 = U.row(v2);
                const Eigen::Vector2d ui3 = U.row(v3);

                const auto cot2 = 2.0 * C_(fid, (i + 1) % 3);
                const auto cot3 = 2.0 * C_(fid, (i + 2) % 3);

                const Eigen::Vector2d u23_bot = R * (ui3 - ui2);

                const Eigen::Vector2d gi = (const_term * u23_bot) * _2 + _1 * (-cot2 * ui3 - cot3 * ui2 + (cot2 + cot3) * ui1);
                g(v1 * 2 + 0) += gi(0);
                g(v1 * 2 + 1) += gi(1);
            }
        }
    }

    DGP_INLINE bool BijectiveParameterization::non_degenerate(const Eigen::MatrixXd &U) const
    {
        Eigen::VectorXd dlbA; 
        igl::doublearea(U, F_, dlbA);
        return dlbA.minCoeff() > 0.0;
    }

    DGP_INLINE double BijectiveParameterization::maximum_step_length(const Eigen::MatrixXd &U, const Eigen::MatrixXd &Delta) const 
    {
        const auto quadraticEquationSolver = [](double a, double b, double c) -> double
        {
            const double delta = b * b - 4.0 * a * c;
            if(delta < 0.0)
            {
                return 0.0;
            }

            if(delta < 1e-8)
            {
                return -b / (2.0 * a);
            }

            const double sqrt_delta = std::sqrt(delta);
            const double x1 = (-b + sqrt_delta) / (2.0 * a);
            const double x2 = (-b - sqrt_delta) / (2.0 * a);
            const double min_one = std::min(x1, x2);
            return min_one < 0.0 ? std::max(x1, x2) : min_one;
        };

        double max_lambda = std::numeric_limits<double>::lowest();
        for(int fid = 0; fid < nf_; ++fid)
        {
            const auto &f = F_.row(fid);
            const Eigen::Vector2d u1 = U.row(f(0));
            const Eigen::Vector2d u2 = U.row(f(1));
            const Eigen::Vector2d u3 = U.row(f(2));
            const Eigen::Vector2d u21 = u2 - u1;
            const Eigen::Vector2d u31 = u3 - u1;

            const Eigen::Vector2d du1 = Delta.row(f(0));
            const Eigen::Vector2d du2 = Delta.row(f(1));
            const Eigen::Vector2d du3 = Delta.row(f(2));
            const Eigen::Vector2d du21 = du2 - du1;
            const Eigen::Vector2d du31 = du3 - du1;

            const auto lambda = quadraticEquationSolver(
                du21(0) * du31(1) - du21(1) * du31(0), 
                u21(0) * du31(1) + du21(0) * u31(1) - u31(0) * du21(1) - du31(0) * u21(1), 
                u21(0) * u31(1) - u21(1) * u31(0));
            if(lambda <= 0.0)
            {
                std::cout << lambda << std::endl;
            }
            max_lambda = std::max(max_lambda, lambda); 
        }
        return max_lambda;
    }

    DGP_INLINE double BijectiveParameterization::calculate_step_length(const Eigen::MatrixXd &U, const Eigen::MatrixXd &Delta, double max_lambda) const
    {
        const double c = 1e-4; 
        const double phi = 1;
        const auto e0 = energy_dirichlet(U);
        while(max_lambda > 1e-8)
        {
            const auto U1 = U - max_lambda * Delta;
            const auto e1 = energy_dirichlet(U1);
            if(e1 < e0 + c * max_lambda && non_degenerate(U1))
            {
                return max_lambda;
            }
            max_lambda *= 0.5;
        }
        return 0.0;
    }

    DGP_INLINE void BijectiveParameterization::per_iteration(int iteration, double &g_norm, double &lambda)
    {
        Eigen::MatrixXd U = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(U_.data(), nv_, 2);
        Eigen::VectorXd g0;
        calculate_gradient(U, g0);
        
        Eigen::VectorXd r;
        const auto &k = iteration;
        const auto &n = number_of_LBFGS_variables_;
        auto q = g0; 
        for(int i = k - 1; i >= std::max(0, k - n); --i)
        {
            alpha_(i % n) = q.dot(S_k_.row(i % n)) * rho_(i % n);
            q -= alpha_(i % n) * Y_k_.row(i % n);
        }

        if(k == 0)
        {
            r = q; 
        }
        else 
        {
            r = S_k_.row((k - 1) % n).dot(Y_k_.row((k - 1) % n)) / Y_k_.row((k - 1) % n).dot(Y_k_.row((k - 1) % n)) * q; 
        }

        for(int i = std::max(0, k - n); i < k; ++i)
        {
            const auto beta = rho_(i % n) * Y_k_.row(i % n).dot(r);
            r += (alpha_(i % n) - beta) * S_k_.row(i % n);
        }

        Eigen::MatrixXd Delta = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(r.data(), nv_, 2);
        // double max_lambda = k==0 ? 1e-5 : maximum_step_length(U, Delta);
        double max_lambda = maximum_step_length(U, Delta);
        lambda = calculate_step_length(U, Delta, max_lambda);

        S_k_.row(iteration % n) = lambda * -r; 
        U_ += S_k_.row(iteration % n);
        Eigen::VectorXd g1;
        U = Eigen::Map<Eigen::MatrixXd>(U_.data(), nv_, 2);
        calculate_gradient(U, g1);
        Y_k_.row(iteration % n) = g1 - g0;
        rho_(iteration % n) = 1.0 / S_k_.row(iteration % n).dot(Y_k_.row(iteration % n));
        g_norm = g1.norm();
    }

    DGP_INLINE void BijectiveParameterization::output(Eigen::MatrixXd &U) const
    {
        U = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(U_.data(), nv_, 2);
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif