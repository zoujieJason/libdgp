#include "filling_holes.h"
#include "vertex_edge_adjacency.h"

#include <igl/boundary_loop.h>

#include <limits>
#include <vector>

//#define FILLING_HOLES_DEBUG

namespace dgp
{
	template<
		typename DerivedV,
		typename DerivedF,
		typename DerivedVr,
		typename DerivedFr>
	bool filling_holes(
		const Eigen::PlainObjectBase<DerivedV>& V,
		const Eigen::PlainObjectBase<DerivedF>& F,
		Eigen::PlainObjectBase<DerivedVr>& Vr, 
		Eigen::PlainObjectBase<DerivedFr>& Fr)
	{
		assert(V.cols() == 3); 
		assert(F.cols() == 3);

		const size_t nv = V.rows();
		const size_t nf = F.rows();

		using Scalar = typename DerivedV::Scalar;

		using LoopType  = std::vector<int>;
		using LoopsType = std::vector<LoopType>;
		using MatrixXS  = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
		using MatrixXI  = Eigen::MatrixXi;
		using Vector3S  = Eigen::Matrix<Scalar, 3, 1>;
		using Vector2I  = Eigen::Matrix<int, 2, 1>;

		using FrScalar   = typename DerivedFr::Scalar;
		using FacesCache = std::vector<FrScalar>;
		using FrMapType  = Eigen::Matrix<FrScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

		LoopsType L;
		igl::boundary_loop(F, L);

		Eigen::MatrixXi E, uE;
		Eigen::VectorXi EMAP;
		std::vector<std::vector<size_t> > uE2E;
		igl::unique_edge_map(F, E, uE, EMAP, uE2E);

		std::vector<std::vector<size_t> > VuE;
		vertex_edge_adjacency(E, EMAP, nv, uE.rows(), VuE);

		const Scalar SMAX = std::numeric_limits<Scalar>::max();

		const auto e2f = [&](const int &vi, const int &vj)->int
		{
			for (size_t i = 0; i < E.rows(); ++i)
			{
				if ((E(i, 0) == vi && E(i, 1) == vj) || (E(i, 0) == vj && E(i, 1) == vi))
				{
					return i % nf;
				}
			}
			return -1;
		};

		const auto findTriangleVertex = [&](const int &fid, const int &vi, const int &vj)->int
		{
			if (fid == -1)
			{
				std::cerr << "findTriangleVertex: fid == -1" << std::endl;
				return -1;
			}

			size_t cnt = 0;
			int id = -1;
			for (size_t i = 0; i < 3; ++i)
			{
				if (F(fid, i) == vi || F(fid, i) == vj)
				{
					cnt++;
				}
				else
				{
					id = F(fid, i);
				}
			}
			return id;
		};

		const auto isInteriorEdge = [&](const int &vi, const int &vj)
		{
			for (const auto &ueid : VuE[vi])
			{
				if ((uE(ueid, 0) == vi && uE(ueid, 1) == vj) || (uE(ueid, 0) == vj && uE(ueid, 1) == vi))
				{
					const auto edges = uE2E[ueid];
					if (edges.size() > 1)
					{
						return true;
					}
				}
			}
			return false;
		};

		const auto initializeDP = [&](const LoopType &loop, MatrixXS &AN, MatrixXS &AR, MatrixXI &I)->bool
		{
			const size_t n = loop.size();
			AN.setConstant(n, n, SMAX);
			AR.setConstant(n, n, SMAX);
			I.setConstant(n, n, -1);
			for (size_t i = 0; i < n - 1; ++i)
			{
				AN(i, i + 1) = 0;
				AR(i, i + 1) = 0;
			}
			for (size_t i = 0; i < n; ++i)
			{
				const size_t _0 = i;
				const size_t _1 = (i + 1) % n;

				const int vi = loop[_0];
				const int vj = loop[_1];

				const auto triangle_vid = findTriangleVertex(e2f(vi, vj), vi, vj);
				if (triangle_vid == -1)
				{
					return false;
				}
				I(_0, _1) = nv + triangle_vid;
			}
			return true;
		};

		const auto calculateNormal = [&](const Vector3S &v0, const Vector3S &v1, const Vector3S &v2, Vector3S &n)-> bool
		{
			n = (v1 - v0).cross(v2 - v0);
			n.normalize();
			return true;
		};

		const auto calculateArea = [&](const Vector3S &v0, const Vector3S &v1, const Vector3S &v2, Scalar &area)-> bool
		{
			const Vector3S n = (v1 - v0).cross(v2 - v0);
			area = 0.5 * n.norm();
			return true;
		};

		const auto triangleNormal = [&](const int &_0, const int &_1, const LoopType &loop, const MatrixXI &I, Vector3S &n)->bool
		{
			int _s = I(_0, _1);
			if (_s == -1)
			{
				return false;
			}
			int _v0, _v1, _v2;
			if (_s < nv)
			{
				_v0 = loop[_0];
				_v1 = loop[_s];
				_v2 = loop[_1];
			}
			else
			{
				_v0 = loop[_0];
				_v1 = _s % nv;
				_v2 = loop[_1];
			}
			const Vector3S v0 = V.row(_v0);
			const Vector3S v1 = V.row(_v1);
			const Vector3S v2 = V.row(_v2);
#ifdef FILLING_HOLES_DEBUG
			std::cout << "Calculate triangle normal: " << _0 << " " << _1 << " I(_0, _1) = " << I(_0, _1) << " | " << _v0 << " " << _v1 << " " << _v2 << std::endl;
#endif
			return calculateNormal(v0, v1, v2, n);
		};

		const auto calculateDihedralAngle = [&](const Vector3S &n0, const Vector3S &n1, Scalar &angle)
		{
			//make sure all angle measures are positive.
			angle = static_cast<Scalar>(1.0) - std::max(static_cast<Scalar>(-1.0), std::min(static_cast<Scalar>(1.0), n0.dot(n1)));
		};

		const auto calculateTriangleWeight = [&](
			const int &i, const int &m, const int &j,
			const LoopType &loop,
			const MatrixXI &I,
			Scalar &angle, Scalar& area)->bool
		{
			const Vector3S v0 = V.row(loop[i]);
			const Vector3S v1 = V.row(loop[m]);
			const Vector3S v2 = V.row(loop[j]);

			calculateArea(v0, v1, v2, area);

			Vector3S n;
			calculateNormal(v0, v1, v2, n);

			Vector3S n_im, n_mj;
			triangleNormal(i, m, loop, I, n_im);
			triangleNormal(m, j, loop, I, n_mj);

			Scalar angle_im, angle_mj;
			calculateDihedralAngle(n, n_im, angle_im);
			calculateDihedralAngle(n, n_mj, angle_mj);

			angle = std::max(angle, angle_im);
			angle = std::max(angle, angle_mj);

			if (i == 0 && j == loop.size() - 1)
			{
				Vector3S n_ij;
				triangleNormal(j, i, loop, I, n_ij);
				Scalar angle_ij;
				calculateDihedralAngle(n, n_ij, angle_ij);
				angle = std::max(angle, angle_ij);
			}
			return true;
		};

		const auto cacheDP = [&](const LoopType &loop, MatrixXS &ANgle, MatrixXS &ARea, MatrixXI &I)->bool
		{
			const size_t n = loop.size();
			if (n < 3)
			{
				return false;
			}
			
			if (!initializeDP(loop, ANgle, ARea, I))
			{
				std::cerr << "!initializeDP(loop, ANgle, ARea, I)" << std::endl;
				return false; 
			}

			//for all vertexes but last two. 
			for (int k = 2; k < n; ++k)
			{
				//for all polygons belong to [i,i+k],
				//e.g. for all single triangle when i == 2.
				for (int i = 0; i < n - k; ++i)
				{
					int j = i + k;
					Scalar min_angle = SMAX;
					Scalar min_area = SMAX;
					int index = -1;

					for (int m = i + 1; m < j; ++m)
					{
						Scalar angle = 0;
						Scalar area = 0;

						if (isInteriorEdge(loop[i], loop[m]) || isInteriorEdge(loop[m], loop[j]) || isInteriorEdge(loop[j], loop[i]))
						{
							angle = SMAX;
							area = SMAX;
						}
						else
						{
							calculateTriangleWeight(i, m, j, loop, I, angle, area);
						}
#ifdef FILLING_HOLES_DEBUG
						std::cout << "Current process: " << i << " " << m << " " << j << std::endl;
						std::cout << "Triangle weight: " << angle << " " << area << std::endl;
						std::cout << "w" << i << m << ": " << ANgle(i, m) << ", w" << m << j << ": " << ANgle(m, j) << std::endl;
#endif
						angle = std::max(angle, ANgle(i, m));
						angle = std::max(angle, ANgle(m, j));
						area = area + ARea(i, m);
						area = area + ARea(m, j);

						if (angle < min_angle || (angle == min_angle && area < min_area))
						{
							min_angle = angle;
							min_area = area;
							index = m;
						}
					}
					ANgle(i, j) = min_angle;
					ARea(i, j) = min_area;
					I(i, j) = index;
#ifdef FILLING_HOLES_DEBUG
					std::cout << "Result: w" << i << j << ", I(i, j) = " << index << std::endl;
					std::cout << std::endl;
#endif
				}
			}

			return true;
		};

		
		const auto fillingHole = [&](const LoopType &loop, const MatrixXI &I, FacesCache &cache)
		{
			Vector2I e0;
			e0 << 0, loop.size() - 1;
			std::vector<Vector2I> edges;
			edges.push_back(e0);
			while (!edges.empty())
			{
				Vector2I edge = edges.back();
				edges.pop_back();

				const int i = edge(0);
				const int k = edge(1);
				const int j = I(i, k);
				if (j == -1 || j >= nv)
				{
					continue;
				}
				if (k - i < 2)
				{
					continue;
				}
				cache.push_back(loop[i]);
				cache.push_back(loop[j]);
				cache.push_back(loop[k]);
				{
					Vector2I edge;
					edge << i, j;
					edges.push_back(edge);
				}

				{
					Vector2I edge;
					edge << j, k;
					edges.push_back(edge);
				}
			}
		};

		FacesCache faces_cache;
		for (auto &loop : L)
		{
			std::reverse(loop.begin(), loop.end());

#ifdef FILLING_HOLES_DEBUG
			std::cout << "Boundary loop: ";
			for (auto val : loop)
			{
				std::cout << val << " ";
			}
			std::cout << std::endl;
#endif
			if (loop.size() < 3)
			{
				std::cerr << "loop.size() < 3 && Impossible error." << std::endl;
				return false;
			}

			MatrixXS ANgle;
			MatrixXS ARea;
			MatrixXI I;
			if (!cacheDP(loop, ANgle, ARea, I))
			{
				return false; 
			}

#ifdef FILLING_HOLES_DEBUG
			std::cout << "ANgle: " << std::endl;
			std::cout << ANgle << std::endl;
			std::cout << std::endl;

			std::cout << "ARea: " << std::endl;
			std::cout << ARea << std::endl;
			std::cout << std::endl;

			std::cout << "I: " << std::endl;
			std::cout << I << std::endl;
#endif 
			fillingHole(loop, I, faces_cache);
		}

		Vr = V;

		const size_t nnf = faces_cache.size() / 3;
		Fr.resize(nf + nnf, 3);
		Fr.block(0, 0, nf, 3) = F;
		Fr.block(nf, 0, nnf, 3) = Eigen::Map<FrMapType>(faces_cache.data(), nnf, 3);
		return true;
	}

	template<
		typename DerivedV,
		typename DerivedF,
		typename DerivedFr>
	bool filling_holes(
		const Eigen::PlainObjectBase<DerivedV>& V,
		const Eigen::PlainObjectBase<DerivedF>& F,
		std::vector<int> &polyline_of_hole,
		Eigen::PlainObjectBase<DerivedFr>& Fr)
	{
		assert(V.cols() == 3); 
		assert(F.cols() == 3);

		const size_t nv = V.rows();
		const size_t nf = F.rows();

		using Scalar = typename DerivedV::Scalar;

		using LoopType = std::vector<int>;
		using LoopsType = std::vector<LoopType>;
		using MatrixXS = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
		using MatrixXI = Eigen::MatrixXi;
		using Vector3S = Eigen::Matrix<Scalar, 3, 1>;
		using Vector2I = Eigen::Matrix<int, 2, 1>;

		using FrScalar = typename DerivedFr::Scalar;
		using FacesCache = std::vector<FrScalar>;
		using FrMapType = Eigen::Matrix<FrScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

		Eigen::MatrixXi E, uE;
		Eigen::VectorXi EMAP;
		std::vector<std::vector<size_t> > uE2E;
		igl::unique_edge_map(F, E, uE, EMAP, uE2E);

		std::vector<std::vector<size_t> > VuE;
		vertex_edge_adjacency(E, EMAP, nv, uE.rows(), VuE);

		const Scalar SMAX = std::numeric_limits<Scalar>::max();

		const auto e2f = [&](const int &vi, const int &vj)->int
		{
			for (size_t i = 0; i < E.rows(); ++i)
			{
				if ((E(i, 0) == vi && E(i, 1) == vj) || (E(i, 0) == vj && E(i, 1) == vi))
				{
					return i % nf;
				}
			}
			return -1;
		};

		const auto findTriangleVertex = [&](const int &fid, const int &vi, const int &vj)->int
		{
			if (fid == -1)
			{
				return -1;
			}

			size_t cnt = 0;
			int id = -1;
			for (size_t i = 0; i < 3; ++i)
			{
				if (F(fid, i) == vi || F(fid, i) == vj)
				{
					cnt++;
				}
				else
				{
					id = F(fid, i);
				}
			}
			return id;
		};

		const auto isInteriorEdge = [&](const int &vi, const int &vj)
		{
			for (const auto &ueid : VuE[vi])
			{
				if ((uE(ueid, 0) == vi && uE(ueid, 1) == vj) || (uE(ueid, 0) == vj && uE(ueid, 1) == vi))
				{
					const auto edges = uE2E[ueid];
					if (edges.size() > 1)
					{
						return true;
					}
				}
			}
			return false;
		};

		const auto initializeDP = [&](const LoopType &loop, MatrixXS &AN, MatrixXS &AR, MatrixXI &I)->bool
		{
			const size_t n = loop.size();
			AN.setConstant(n, n, SMAX);
			AR.setConstant(n, n, SMAX);
			I.setConstant(n, n, -1);
			for (size_t i = 0; i < n - 1; ++i)
			{
				AN(i, i + 1) = 0;
				AR(i, i + 1) = 0;
			}
			for (size_t i = 0; i < n; ++i)
			{
				const size_t _0 = i;
				const size_t _1 = (i + 1) % n;

				const int vi = loop[_0];
				const int vj = loop[_1];

				const auto triangle_vid = findTriangleVertex(e2f(vi, vj), vi, vj);
				I(_0, _1) = (triangle_vid == -1) ? -1 : nv + triangle_vid;
			}
			return true;
		};

		const auto calculateNormal = [&](const Vector3S &v0, const Vector3S &v1, const Vector3S &v2, Vector3S &n)-> bool
		{
			n = (v1 - v0).cross(v2 - v0);
			n.normalize();
			return true;
		};

		const auto calculateArea = [&](const Vector3S &v0, const Vector3S &v1, const Vector3S &v2, Scalar &area)-> bool
		{
			const Vector3S n = (v1 - v0).cross(v2 - v0);
			area = 0.5 * n.norm();
			return true;
		};

		const auto triangleNormal = [&](const int &_0, const int &_1, const LoopType &loop, const MatrixXI &I, Vector3S &n)->bool
		{
			int _s = I(_0, _1);
			if (_s == -1)
			{
				return false;
			}
			int _v0, _v1, _v2;
			if (_s < nv)
			{
				_v0 = loop[_0];
				_v1 = loop[_s];
				_v2 = loop[_1];
			}
			else
			{
				_v0 = loop[_0];
				_v1 = _s % nv;
				_v2 = loop[_1];
			}
			const Vector3S v0 = V.row(_v0);
			const Vector3S v1 = V.row(_v1);
			const Vector3S v2 = V.row(_v2);
#ifdef FILLING_HOLES_DEBUG
			std::cout << "Calculate triangle normal: " << _0 << " " << _1 << " I(_0, _1) = " << I(_0, _1) << " | " << _v0 << " " << _v1 << " " << _v2 << std::endl;
#endif
			return calculateNormal(v0, v1, v2, n);
		};

		const auto calculateDihedralAngle = [&](const Vector3S &n0, const Vector3S &n1, Scalar &angle)
		{
			//make sure all angle measures are positive.
			angle = static_cast<Scalar>(1.0) - std::max(static_cast<Scalar>(-1.0), std::min(static_cast<Scalar>(1.0), n0.dot(n1)));
		};

		const auto calculateTriangleWeight = [&](
			const int &i, const int &m, const int &j,
			const LoopType &loop,
			const MatrixXI &I,
			Scalar &angle, Scalar& area)->bool
		{
			const Vector3S v0 = V.row(loop[i]);
			const Vector3S v1 = V.row(loop[m]);
			const Vector3S v2 = V.row(loop[j]);

			calculateArea(v0, v1, v2, area);

			Vector3S n;
			calculateNormal(v0, v1, v2, n);

			Vector3S n_im, n_mj;
			if(!triangleNormal(i, m, loop, I, n_im))
			{
				n_im = n; 
			}
			if(!triangleNormal(m, j, loop, I, n_mj))
			{
				n_mj = n;
			}

			Scalar angle_im, angle_mj;
			calculateDihedralAngle(n, n_im, angle_im);
			calculateDihedralAngle(n, n_mj, angle_mj);

			angle = std::max(angle, angle_im);
			angle = std::max(angle, angle_mj);

			if (i == 0 && j == loop.size() - 1)
			{
				Vector3S n_ij;
				triangleNormal(j, i, loop, I, n_ij);
				Scalar angle_ij;
				calculateDihedralAngle(n, n_ij, angle_ij);
				angle = std::max(angle, angle_ij);
			}
			return true;
		};

		const auto cacheDP = [&](const LoopType &loop, MatrixXS &ANgle, MatrixXS &ARea, MatrixXI &I)->bool
		{
			const size_t n = loop.size();
			if (n < 3)
			{
				return false;
			}

			if (!initializeDP(loop, ANgle, ARea, I))
			{
				std::cerr << "!initializeDP(loop, ANgle, ARea, I)" << std::endl;
				return false;
			}

			//for all vertexes but last two. 
			for (int k = 2; k < n; ++k)
			{
				//for all polygons belong to [i,i+k],
				//e.g. for all single triangle when i == 2.
				for (int i = 0; i < n - k; ++i)
				{
					int j = i + k;
					Scalar min_angle = SMAX;
					Scalar min_area = SMAX;
					int index = -1;

					for (int m = i + 1; m < j; ++m)
					{
						Scalar angle = 0;
						Scalar area = 0;

						if (isInteriorEdge(loop[i], loop[m]) || isInteriorEdge(loop[m], loop[j]) || isInteriorEdge(loop[j], loop[i]))
						{
							angle = SMAX;
							area = SMAX;
						}
						else
						{
							calculateTriangleWeight(i, m, j, loop, I, angle, area);
						}
#ifdef FILLING_HOLES_DEBUG
						std::cout << "Current process: " << i << " " << m << " " << j << std::endl;
						std::cout << "Triangle weight: " << angle << " " << area << std::endl;
						std::cout << "w" << i << m << ": " << ANgle(i, m) << ", w" << m << j << ": " << ANgle(m, j) << std::endl;
#endif
						angle = std::max(angle, ANgle(i, m));
						angle = std::max(angle, ANgle(m, j));
						area = area + ARea(i, m);
						area = area + ARea(m, j);

						if (angle < min_angle || (angle == min_angle && area < min_area))
						{
							min_angle = angle;
							min_area = area;
							index = m;
						}
					}
					ANgle(i, j) = min_angle;
					ARea(i, j) = min_area;
					I(i, j) = index;
#ifdef FILLING_HOLES_DEBUG
					std::cout << "Result: w" << i << j << ", I(i, j) = " << index << std::endl;
					std::cout << std::endl;
#endif
				}
			}

			return true;
		};

		const auto fillingHole = [&](const LoopType &loop, const MatrixXI &I, FacesCache &cache)
		{
			Vector2I e0;
			e0 << 0, loop.size() - 1;
			std::vector<Vector2I> edges;
			edges.push_back(e0);
			while (!edges.empty())
			{
				Vector2I edge = edges.back();
				edges.pop_back();

				const int i = edge(0);
				const int k = edge(1);
				const int j = I(i, k);
				if (j == -1 || j >= nv)
				{
					continue;
				}
				if (k - i < 2)
				{
					continue;
				}
				cache.push_back(loop[i]);
				cache.push_back(loop[j]);
				cache.push_back(loop[k]);
				{
					Vector2I edge;
					edge << i, j;
					edges.push_back(edge);
				}

				{
					Vector2I edge;
					edge << j, k;
					edges.push_back(edge);
				}
			}
		};

		std::reverse(polyline_of_hole.begin(), polyline_of_hole.end());

#ifdef FILLING_HOLES_DEBUG
		std::cout << "Boundary loop: ";
		for (auto val : polyline_of_hole)
		{
			std::cout << val << " ";
		}
		std::cout << std::endl;
#endif
		if (polyline_of_hole.size() < 3)
		{
			std::cerr << "loop.size() < 3 && Impossible error." << std::endl;
			return false;
		}

		MatrixXS ANgle;
		MatrixXS ARea;
		MatrixXI I;
		if (!cacheDP(polyline_of_hole, ANgle, ARea, I))
		{
			return false;
		}

#ifdef FILLING_HOLES_DEBUG
		std::cout << "ANgle: " << std::endl;
		std::cout << ANgle << std::endl;
		std::cout << std::endl;

		std::cout << "ARea: " << std::endl;
		std::cout << ARea << std::endl;
		std::cout << std::endl;

		std::cout << "I: " << std::endl;
		std::cout << I << std::endl;
#endif
		FacesCache faces_cache;
		fillingHole(polyline_of_hole, I, faces_cache);

		const size_t nnf = faces_cache.size() / 3;
		Fr.resize(nf + nnf, 3);
		Fr.block(0, 0, nf, 3) = F;
		Fr.block(nf, 0, nnf, 3) = Eigen::Map<FrMapType>(faces_cache.data(), nnf, 3);
		return true;
	}
}