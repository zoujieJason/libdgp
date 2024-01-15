#ifndef LIBDGP_KDTREE_H
#define LIBDGP_KDTREE_H

#include "dgp_inline.h"

#include <Eigen/Dense>

#include <cassert>
#include <type_traits>
#include <vector>
#include <variant>

namespace dgp
{
    template<typename PlainObject>
    class KDTree
    {
    private:
        using Scalar        = typename PlainObject::Scalar;
        using RowVectorType = typename Eigen::internal::plain_row_type<PlainObject>::type;

        struct Node
        {
            Scalar split_value;
            int dim;
            int child; 
        };

        struct Leaf
        {
            int start;
            int size;
            Leaf(int _1, int _2) : start(_1), size(_2) {}
        };
        using TreeNode = std::variant<Node, Leaf>;

        PlainObject V_;
        RowVectorType min_; 
        RowVectorType max_;

        int max_depth_; 
        int leaf_capacity_;
        int dim_;
        std::vector<TreeNode> tree_nodes_;

    public:
        template<typename Derived>
        KDTree(const Eigen::MatrixBase<Derived> &V, int leaf_capacity = 16, int max_depth = 64)
            : V_{V}, leaf_capacity_{leaf_capacity}, max_depth_{max_depth}
        {
            dim_ = V_.cols();
            min_ = V_.colwise().minCoeff();
            max_ = V_.colwise().maxCoeff();
            tree_nodes_.reserve(4 * V_.rows() / leaf_capacity_);
            tree_nodes_.resize(1);
            tree_nodes_[0] = Node{};
            create_tree(0, 0, V_.rows(), 1);
        }
        
        template<typename Derived>
        DGP_INLINE void brutal_search(const Eigen::MatrixBase<Derived> &q, Eigen::MatrixBase<Derived> &p, typename Derived::Scalar &min_dist_sq)
        {
            static_assert(std::is_same<typename Derived::Scalar, Scalar>::value, "MIXED SCALAR TYPES");
            EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived, 1, Eigen::Dynamic)  
            assert(q.size() == dim_);

            min_dist_sq = std::numeric_limits<Scalar>::max();
            for(int i = 0; i < V_.rows(); ++i)
            {
                Scalar dist_sq = (q - V_.row(i)).squaredNorm();
                if(dist_sq < min_dist_sq)
                {
                    min_dist_sq = dist_sq;
                    p = V_.row(i);
                }
            }
        }
        
        template<typename Derived>
        DGP_INLINE int recursive_search(const Eigen::MatrixBase<Derived> &q, Eigen::MatrixBase<Derived> &p, typename Derived::Scalar &min_dist_sq)
        {
            static_assert(std::is_same<typename Derived::Scalar, Scalar>::value, "MIXED SCALAR TYPES");
            EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived, 1, Eigen::Dynamic)  
            assert(q.size() == dim_);

            int number_of_reaching_leaves = 0;
            min_dist_sq = std::numeric_limits<Scalar>::max();
            recursive_search(q, 0, min_dist_sq, p, number_of_reaching_leaves);
            return number_of_reaching_leaves;
        }

        template<typename Derived>
        DGP_INLINE int recursive_pruning_search(const Eigen::MatrixBase<Derived> &q, Eigen::MatrixBase<Derived> &p, typename Derived::Scalar &min_dist_sq)
        {
            static_assert(std::is_same<typename Derived::Scalar, Scalar>::value, "MIXED SCALAR TYPES");
            EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived, 1, Eigen::Dynamic)  
            assert(q.size() == dim_);

            int number_of_reaching_leaves = 0;
            min_dist_sq = std::numeric_limits<Scalar>::max();
            recursive_search(min_, max_, q, 0, min_dist_sq, p, number_of_reaching_leaves);
            return number_of_reaching_leaves; 
        }

        template<typename Derived>
        DGP_INLINE int iterative_search(const Eigen::MatrixBase<Derived> &q, Eigen::MatrixBase<Derived> &p, typename Derived::Scalar &min_dist_sq)
        {
            static_assert(std::is_same<typename Derived::Scalar, Scalar>::value, "MIXED SCALAR TYPES");
            EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived, 1, Eigen::Dynamic)  
            assert(q.size() == dim_);

            using QueryNode = std::pair<int, Scalar>;
            std::vector<QueryNode> stack(max_depth_); 
            stack[0].first = 0;
            stack[0].second = 0.0;
            int number_of_reaching_leaves = 0; 

            int probe = 1 ;
            p = V_.row(V_.rows() / 2);
            min_dist_sq = (q - p).squaredNorm(); 

            while(probe)
            {
                auto &stack_curr = stack[probe - 1];
                auto &stack_next = stack[probe];
                auto &tree_node = tree_nodes_[stack_curr.first];

                if(stack_curr.second < min_dist_sq)
                {
                    if(std::holds_alternative<Leaf>(tree_node))
                    {
                        --probe;
                        const Leaf &leaf = std::get<Leaf>(tree_node);
                        for(int i = leaf.start; i < leaf.start + leaf.size; ++i)
                        {
                            Scalar dist = (q - V_.row(i)).squaredNorm();
                            if(dist < min_dist_sq)
                            {
                                min_dist_sq = dist;
                                p = V_.row(i);
                            }
                        }
                        number_of_reaching_leaves++;
                    }
                    else
                    {
                        const Node &node = std::get<Node>(tree_node);
                        Scalar dist = q[node.dim] - node.split_value;
                        if(dist < 0.0)
                        {
                            stack_next.first = node.child;
                            stack_curr.first = node.child + 1;
                        }
                        else 
                        {
                            stack_next.first = node.child + 1;
                            stack_curr.first = node.child; 
                        }
                        stack_next.second = stack_curr.second;
                        stack_curr.second = dist * dist;
                        ++probe; 
                    }
                }
                else 
                {
                    --probe;
                }
            }
            return number_of_reaching_leaves;
        }

    private:
        DGP_INLINE int split(int start, int end, int dim, Scalar split_value)
        {
            int left_probe(start), right_probe(end - 1);
            for(;left_probe<right_probe; ++left_probe, --right_probe)
            {
                while(left_probe < end && V_(left_probe, dim) < split_value)
                {
                    ++left_probe;
                }

                while(right_probe >= start && V_(right_probe, dim) >= split_value)
                {
                    --right_probe;
                }

                if (left_probe > right_probe)
                {
                    break;
                }

                RowVectorType v = V_.row(left_probe);
                V_.row(left_probe) = V_.row(right_probe);
                V_.row(right_probe) = v;
            }
            return V_(left_probe, dim) < split_value ? left_probe + 1 : left_probe;
        }

        DGP_INLINE void create_tree(int nid, int start, int end, int level)
        {
            const RowVectorType &min = V_.block(start, 0, end - start, dim_).colwise().minCoeff();
            const RowVectorType &max = V_.block(start, 0, end - start, dim_).colwise().maxCoeff();

            Node &node = std::get<Node>(tree_nodes_[nid]);
            (max - min).maxCoeff(&node.dim);
            node.split_value = 0.5 * (min[node.dim] + max[node.dim]);

            const auto mid = split(start, end, node.dim, node.split_value);
            node.child = tree_nodes_.size();
            int child = node.child; 
            //  the resize operation may invalidate the reference to the node
            tree_nodes_.resize(tree_nodes_.size() + 2);
            
            if(mid - start <= leaf_capacity_ || level >= max_depth_)
            {
                tree_nodes_[child] = Leaf{start, mid - start};
            }
            else 
            {
                tree_nodes_[child] = Node{};
                create_tree(child, start, mid, level + 1);
            }

            if(end - mid <= leaf_capacity_ || level >= max_depth_)
            {
                tree_nodes_[child + 1] = Leaf{mid, end - mid};
            }
            else 
            {
                tree_nodes_[child + 1] = Node{};
                create_tree(child + 1, mid, end, level + 1);
            }
        }

        template<typename Derived>
        DGP_INLINE typename Derived::Scalar distance(const Eigen::MatrixBase<Derived> &q, const TreeNode &tree_node, typename Derived::Scalar &min_dist_sq, Eigen::MatrixBase<Derived> &p, int &cnt)
        {
            if(std::holds_alternative<Leaf>(tree_node))
            {
                const Leaf &leaf = std::get<Leaf>(tree_node);   
                for(int i = leaf.start; i < leaf.start + leaf.size; ++i)
                {
                    Scalar dist_sq = (q - V_.row(i)).squaredNorm();
                    if(dist_sq < min_dist_sq)
                    {
                        min_dist_sq = dist_sq;
                        p = V_.row(i);
                    }
                }
                cnt++;
                return min_dist_sq;
            }
            else 
            {
                const Node &node = std::get<Node>(tree_node);
                Scalar dist = q[node.dim] - node.split_value;
                return dist * dist;
            }
        }

        template<typename Derived>
        DGP_INLINE void recursive_search(const Eigen::MatrixBase<Derived> &q, int nid, typename Derived::Scalar &min_dist_sq, Eigen::MatrixBase<Derived> &p, int &cnt)
        {
            auto &tree_node = tree_nodes_[nid];
            auto dist_sq = distance(q, tree_node, min_dist_sq, p, cnt);

            if(std::holds_alternative<Node>(tree_node))
            {
                const Node &node = std::get<Node>(tree_node);
                if(q[node.dim] < node.split_value)
                {
                    recursive_search(q, node.child, min_dist_sq, p, cnt);
                    if(dist_sq < min_dist_sq)
                    {
                        recursive_search(q, node.child + 1, min_dist_sq, p, cnt);
                    }
                }
                else 
                {
                    recursive_search(q, node.child + 1, min_dist_sq, p, cnt);
                    if(dist_sq < min_dist_sq)
                    {
                        recursive_search(q, node.child, min_dist_sq, p, cnt);
                    }
                }
            }
        }

        template<typename Derived>
        DGP_INLINE Scalar distance(const RowVectorType &min, const RowVectorType &max, const Eigen::MatrixBase<Derived> &q)
        {
            Scalar dist_sq = 0.0;
            for(int i = 0; i < q.size(); ++i)
            {
                if(q[i] < min[i])
                {
                    Scalar diff = min[i] - q[i];
                    dist_sq += diff * diff;
                }
                else if(q(i) > max(i))
                {
                    Scalar diff = q[i] - max[i];
                    dist_sq += diff * diff;
                }
            }
            return dist_sq;
        }

        template<typename Derived>
        DGP_INLINE void recursive_search(RowVectorType min, RowVectorType max, const Eigen::MatrixBase<Derived> &q, int nid, Scalar &min_dist_sq, Eigen::MatrixBase<Derived> &p, int &cnt)
        {
            auto &tree_node = tree_nodes_[nid];
            if(std::holds_alternative<Leaf>(tree_node))
            {
                const Leaf &leaf = std::get<Leaf>(tree_node);   
                for(int i = leaf.start; i < leaf.start + leaf.size; ++i)
                {
                    Scalar dist_sq = (q - V_.row(i)).squaredNorm();
                    if(dist_sq < min_dist_sq)
                    {
                        min_dist_sq = dist_sq;
                        p = V_.row(i);
                    }
                }
                cnt++;
            }
            else if(std::holds_alternative<Node>(tree_node))
            {
                const Node &node = std::get<Node>(tree_node);
                if(q[node.dim] < node.split_value)
                {
                    {       
                        auto left_max = max;             
                        left_max[node.dim] = node.split_value;
                        recursive_search(min, left_max, q, node.child, min_dist_sq, p, cnt);
                    }

                    {
                        auto right_min = min;
                        right_min[node.dim] = node.split_value;
                        const auto dist_sq = distance(right_min, max, q);
                        if(dist_sq < min_dist_sq)
                        {
                            recursive_search(right_min, max, q, node.child + 1, min_dist_sq, p, cnt);
                        }
                    }
                }
                else 
                {
                    {
                        auto right_min = min;
                        right_min[node.dim] = node.split_value;
                        recursive_search(right_min, max, q, node.child + 1, min_dist_sq, p, cnt);
                    }
                    
                    {
                        auto left_max = max;             
                        left_max[node.dim] = node.split_value;
                        const auto dist_sq = distance(min, left_max, q);
                        if(dist_sq < min_dist_sq)
                        {
                            recursive_search(min, left_max, q, node.child, min_dist_sq, p, cnt);
                        }
                    }
                }
            }
        }
    };

    template<typename Derived>
    KDTree(const Eigen::MatrixBase<Derived> &) -> KDTree<typename Derived::PlainObject>;
}

#endif //LIBDGP_KDTREE_H