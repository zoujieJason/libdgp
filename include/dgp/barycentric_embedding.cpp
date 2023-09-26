#include "barycentric_embedding.h"
#include "slice_laplacian.h"

#include <igl/cotmatrix.h>

namespace dgp
{
    DGP_INLINE bool barycentric_embedding(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const std::vector<size_t> &mesh_boundary,
        const Eigen::MatrixXd &embedding_boundary,
        Eigen::MatrixXd &UV)
    {
        const auto number_of_boundary = mesh_boundary.size();
        if(number_of_boundary < 3 || number_of_boundary != embedding_boundary.rows() || embedding_boundary.cols() != 2)
        {
            return false;
        }

        Eigen::SparseMatrix<double> L; 
        igl::cotmatrix(V, F, L);

        Eigen::SparseMatrix<double> LA, LB;
        std::vector<bool> is_boundary;
        std::vector<size_t> vertices_map; 
        if(!slice_laplacian(L, mesh_boundary, LA, LB, is_boundary, vertices_map))
        {
            return false; 
        }

        Eigen::MatrixXd b;
        b.resize(number_of_boundary, 2);
        for(size_t i = 0; i < number_of_boundary; ++i)
        {
            b.row(vertices_map[mesh_boundary[i]]) = embedding_boundary.row(i); 
        }

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(LA.transpose() * LA);
		if(solver.info() != Eigen::Success) 
		{
			std::cerr << "LeastSquareProblemSolver error: matrix decomposition failed." << std::endl;
			return false;
		}

		const auto result = solver.solve(LA.transpose() * (-1.0 * LB * b));
		if(solver.info() != Eigen::Success) 
		{
			std::cerr << "LeastSquareProblemSolver error: solving failed." << std::endl;
			return false;
		}

        const auto number_of_vertices = V.rows();
		UV.resize(number_of_vertices, 2);
		for(size_t vid = 0; vid < number_of_vertices; ++vid) 
		{
			if (is_boundary[vid])
			{
				UV.row(vid) = b.row(vertices_map[vid]);
			}
			else
			{
				UV.row(vid) = result.row(vertices_map[vid]);
			}
		}
        return true;
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif