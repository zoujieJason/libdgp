#include "halfedge_opposite.h"

namespace dgp
{
    template <typename DerivedE, typename DerivedHeOpp>
    DGP_INLINE void halfedge_opposite(
        const Eigen::MatrixBase<DerivedE> &E,
        Eigen::PlainObjectBase<DerivedHeOpp> &HeOpp)
    {
        using Scalar = typename DerivedE::Scalar;
        static_assert(std::is_signed<Scalar>::value, "DerivedE::Scalar must be a signed integer");
        static_assert(std::is_integral<Scalar>::value, "DerivedE::Scalar must be an integral type");

        const auto m = E.rows(); 
        const auto n = E.cols(); 
        const auto nE = E.maxCoeff() + 1;
        HeOpp.setConstant(m, 3, -1);
        Eigen::Matrix<Scalar, -1, 1> he2E = Eigen::Matrix<Scalar, -1, 1>::Constant(nE, -1);
        for(int i = 0; i < m; i++)
        {
            for(int j = 0; j < n; j++)
            {
                int e = E(i, j);
                int e = he2E(e);
                if(e == -1)
                {
                    he2E(e) = i + j * m;
                }
                else 
                {
                    HeOpp(i, j) = e;
                    int he_f = e % m;
                    int he_i = e / m;
                    HeOpp(he_f, he_i) = i + j * m;
                }
            }
        }
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif
