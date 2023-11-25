#include "mesh_tree.h"

#include <igl/facet_adjacency_matrix.h>

#include <queue>

namespace dgp
{
    void MeshTree::build(int nv, const Eigen::MatrixXi &F)
    {
        assert(F.cols() == 3 && "MeshTree::build(): F must be #F by 3");
        const auto nf = F.rows(); 
        for(int fid = 0; fid < nf; ++fid)
        {
            nodes_.push_back(std::make_shared<Node>(fid));
        } 

		Eigen::SparseMatrix<int> A;
		igl::facet_adjacency_matrix(F, A);
        std::vector<int> is_visited(nf, false);
        const auto connectedComponent = [&](int seed, int cid)
        {
            std::queue<int> Q;
            Q.push(seed);
            is_visited[seed] = true;
            while(!Q.empty())
            {
                const auto fid = Q.front();
                Q.pop();

                std::weak_ptr<Node> node_parent = nodes_[fid];
                node_parent.lock()->cid_ = cid;
                for(Eigen::SparseMatrix<int>::InnerIterator it(A, fid); it; ++it)
                {
                    const auto adj_fid = it.row();
                    if(!is_visited[adj_fid])
                    {
                        is_visited[adj_fid] = true;
                        Q.push(adj_fid);
                        std::weak_ptr<Node> node_child = nodes_[adj_fid];
                        node_child.lock()->parent_ = node_parent;
                        node_parent.lock()->children_.push_back(node_child);
                    }
                }
            }
        };

        for(int fid = 0, cid = 0; fid < nf; ++fid)
        {
            if(is_visited[fid])
            {
                continue;
            }

            heads_.push_back(nodes_[fid]);
            connectedComponent(fid, cid++);
        }
    }

    void MeshTree::traverse(std::weak_ptr<Node> node, const std::function<void(std::weak_ptr<Node>)> &apply) const
    {
        apply(node);
        for(auto child : node.lock()->children_)
        {
            traverse(child, apply);
        }
    }

    std::vector<std::weak_ptr<MeshTree::Node>>::iterator MeshTree::begin()
    {
        return heads_.begin();
    }

    std::vector<std::weak_ptr<MeshTree::Node>>::iterator MeshTree::end()
    {
        return heads_.end();
    }
}