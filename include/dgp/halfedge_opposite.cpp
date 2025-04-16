#include "halfedge_opposite.h"

#include <igl/orient_halfedges.h>

namespace dgp
{
    template <typename DerivedF, typename DerivedHeOpp>
    DGP_INLINE void halfedge_opposite(
        const Eigen::MatrixBase<DerivedF> &F,
        Eigen::PlainObjectBase<DerivedHeOpp> &HeOpp)
    {
        using Scalar = typename DerivedF::Scalar;
        static_assert(std::is_signed<Scalar>::value, "DerivedF::Scalar must be a signed integer");
        static_assert(std::is_integral<Scalar>::value, "DerivedF::Scalar must be an integral type");

        DerivedF E, oE; 
        igl::orient_halfedges(F, E, oE);

        const auto nf = E.rows(); 
        const auto dim = E.cols(); 
        const auto nE = E.maxCoeff() + 1; // Unique
        HeOpp.setConstant(nf, 3, -1);
        Eigen::Matrix<Scalar, -1, 1> E2He = Eigen::Matrix<Scalar, -1, 1>::Constant(nE, -1);
        for(int i = 0; i < nf; i++)
        {
            for(int j = 0; j < dim; j++)
            {
                int he = i + j * nf;;
                int e = E(i, j);
                if(E2He(e) == -1)
                {
                    E2He(e) = he;
                }
                else 
                {
                    HeOpp(i, j) = E2He(e);
                    int he_f = E2He(e) % nf;
                    int he_i = E2He(e) / nf;
                    HeOpp(he_f, he_i) = he;
                }
            }
        }
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif
