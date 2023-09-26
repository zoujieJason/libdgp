#include "slice_laplacian.h"

namespace dgp
{
    DGP_INLINE bool slice_laplacian(
        const Eigen::SparseMatrix<double> &L,
        const std::vector<size_t> &boundary,
        Eigen::SparseMatrix<double> &LA,
        Eigen::SparseMatrix<double> &LB,
        std::vector<bool> &is_boundary, 
        std::vector<size_t> &vertices_map)
    {
        const auto number_of_vertices = L.rows();
        if(number_of_vertices != L.cols())
        {
            return false; 
        }

        const auto number_of_boundary = boundary.size();
        if(number_of_vertices <= number_of_boundary)
        {
            return false; 
        }

       is_boundary.resize(number_of_vertices, false);
        for(const auto &i: boundary)
        {
            if(i >= number_of_vertices)
            {
                return false; 
            }

            is_boundary[i] = true; 
        }

        vertices_map.resize(number_of_vertices, 0);
        {
            size_t fixed_vid = 0; 
            size_t free_vid = 0; 
            for(size_t i = 0; i < number_of_vertices; ++i)
            {
                if(is_boundary[i])
                {
                    vertices_map[i] = fixed_vid++; 
                }
                else
                {
                    vertices_map[i] = free_vid++; 
                }
            }
        }

        const auto number_of_nonzeros = L.nonZeros();
        const auto number_of_fixed_vertices = number_of_boundary;
		const auto number_of_free_vertices = number_of_vertices - number_of_fixed_vertices;
        std::vector<Eigen::Triplet<double> > LA_triplets, LB_triplets; 
        LA_triplets.reserve(number_of_nonzeros * number_of_free_vertices / number_of_vertices);
		LB_triplets.reserve(number_of_nonzeros * number_of_fixed_vertices / number_of_vertices);
        for (int k = 0; k < L.outerSize(); ++k) 
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) 
			{
				if (is_boundary[it.row()])
				{
					continue;
				}

				if (is_boundary[it.col()])
				{
					LB_triplets.emplace_back(it.row(), vertices_map[it.col()], it.value());
				}
				else
				{
					LA_triplets.emplace_back(it.row(), vertices_map[it.col()], it.value());
				}
			}
		}

		LA.resize(number_of_vertices, number_of_free_vertices);
		LA.setFromTriplets(LA_triplets.begin(), LA_triplets.end());
		LB.resize(number_of_vertices, number_of_fixed_vertices);
		LB.setFromTriplets(LB_triplets.begin(), LB_triplets.end());
        return true;
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif