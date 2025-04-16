#ifndef LIBDGP_SMOOTHER_H
#define LIBDGP_SMOOTHER_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <vector>
#include <memory>
#include <string>

namespace dgp
{
	namespace alg
	{
		struct PrecomputeParam
		{
			int fuzzy_kring{ 6 };
			int block_kring{ 3 };
			bool use_intrinsic{ false };
			std::vector<size_t> remove_vertices; 
			std::vector<size_t> block_vertices;
		};

		struct SDFParam
		{
			double t{ 0.01 };
			bool preserve_source_normals{ true };
			bool hard_constriant{ true };
			double level_set{ 0.2 };
			double consistency_ratio{ 0.8 };
		};

		struct VectorTransportParam
		{
			double t{ 0.5 };
			bool preserve_source_normals{ false };
			int min_spacing_source_normals{ 4 };
			std::vector<size_t> block_vertices;
		};

		struct DiffusionParam
		{
			double lambda{ 0.5 };
			double adjacent_point_potential{ 0.01 };
			int min_spacing_constriants{ 3 };
			std::vector<std::pair<size_t, double>> block_vertices;
			std::vector<std::pair<size_t, double>> potential_vertices;
			double lambda_regularization{0.1};
		};

		class SmootherImpl;

		class Smoother
		{
		public:
			Smoother();

			~Smoother();

			//Warning: max(labels) < sizeof(int) * 8;
			void precompute(
				const std::vector<double>& vbuffer,
				const std::vector<int>& fbuffer,
				const std::vector<int>& faces_label,
				PrecomputeParam &param);
			
			void precompute(
				const Eigen::MatrixXd &V, 
				const Eigen::MatrixXi &F, 
				const Eigen::VectorXi &faces_label,
				PrecomputeParam &param);

			void sdf_required(const SDFParam& param = SDFParam());

			void vector_field_solve(const VectorTransportParam& param = VectorTransportParam());

			void scalar_field_solve(const DiffusionParam& param = DiffusionParam());

			void remesh_along_isoline(
				double phi,
				std::vector<double>& vbuffer,
				std::vector<int>& fbuffer,
				std::vector<int>& labels,
				std::vector<int>& fmap);

			void relaxing(
				double collapse_lenRatio,
				std::vector<double>& vbuffer,
				std::vector<int>& fbuffer,
				std::vector<int>& labels,
				std::vector<int>& fmap);

			void optimize_adjacency(
				std::vector<double>& vbuffer,
				std::vector<int>& fbuffer,
				std::vector<int>& labels);
				
			void get_mesh(size_t label, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

			void get_sdf(size_t label, Eigen::VectorXd &sdf); 

			void get_flow(size_t label, Eigen::MatrixXd &flow);

			void get_phi(size_t label, Eigen::VectorXd &phi); 

			void get_faces_label(size_t label, Eigen::VectorXi &labels); 

			void load(
				const Eigen::MatrixXd &V, 
				const Eigen::MatrixXi &F, 
				const Eigen::VectorXi &faces_label);

			bool precompute(
				size_t label,
				const PrecomputeParam &param = PrecomputeParam());
			
			bool sdf_solve(
				size_t label,
				const SDFParam &param = SDFParam());

			bool vector_field_solve(
				size_t label,
				const VectorTransportParam& param = VectorTransportParam());

			bool scalar_field_solve(
				size_t label,
				const DiffusionParam& param = DiffusionParam());

			bool local_optimize(
				size_t label, 
				const std::vector<size_t> &vertices);

		private:
			std::unique_ptr<SmootherImpl> impl_;
		};
	}
}

#ifndef DGP_STATIC_LIBRARY
#include   "smoother.cpp"
#endif

#endif //LIBDGP_SMOOTHER_H