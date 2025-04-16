#include "snaking.h"
#include "intersection2d.h"
#include "halfedge_opposite.h"
#include "vertex_halfedge_adjacency.h"

#include <igl/vertex_triangle_adjacency.h>
#include <igl/unique_edge_map.h>
#include <igl/orient_halfedges.h>

namespace dgp
{
    template <typename stdit>
    DGP_INLINE void extract_snaxls(
        const Eigen::MatrixXd &V, 
        stdit begin, 
        stdit end,
        Eigen::MatrixXd &iV,
        Eigen::MatrixXi &iE,
        bool loop)
    {
        auto n = std::distance(begin, end); 
        iV = Eigen::MatrixXd::Zero(n, 3); 
        for(stdit it = begin ; it != end; ++it)
        {
            if(it->GetSnaxlType() == Snaxl::Type::Vertex)
            {
                iV.row(std::distance(begin, it)) = V.row(it->AsVertexSnaxl().v);
            }
            else if(it->GetSnaxlType() == Snaxl::Type::Edge)
            {
                auto snaxl = it->AsEdgeSnaxl(); 
                iV.row(std::distance(begin, it)) = (1. - snaxl.t) * V.row(snaxl.v0) + snaxl.t * V.row(snaxl.v1); 
            }
        }
        
        auto iEn = n - (loop ? 0: 1); 
        if(iEn > 0)
        {
            iE = Eigen::MatrixXi::Zero(iEn, 2);
            iE.col(1) = Eigen::VectorXi::LinSpaced(iEn, 1, iEn);
            iE.col(0) = iE.col(1).array() - 1;
            iE.bottomRightCorner(1, 1)(0, 0) = iE.bottomRightCorner(1, 1)(0, 0) % n;
        }
    }

    DGP_INLINE bool snaking_precompute(
        const Eigen::MatrixXi &F,
        SnakingData &snaking_data)
    {
        const auto nf = F.rows(); 
        igl::unique_edge_map(F, snaking_data.E, snaking_data.uE, snaking_data.EMAP, snaking_data.uE2E);
        snaking_data.uE2Fsorted.resize(snaking_data.uE2E.size());
        std::transform(snaking_data.uE2E.begin(), snaking_data.uE2E.end(), snaking_data.uE2Fsorted.begin(), [&nf](auto e)
        {
            std::for_each(e.begin(), e.end(), [&nf](auto &i){ i = i % nf; }); 
            std::sort(e.begin(), e.end());
            return e; 
        });

        halfedge_opposite(F, snaking_data.HeOpp);

        igl::vertex_triangle_adjacency(F.maxCoeff() + 1, F, snaking_data.VF, snaking_data.VFi);

        return vertex_halfedge_adjacency(snaking_data.VF, snaking_data.VFi, F, snaking_data.HeOpp, snaking_data.VHe);
    }

    namespace snaking_utils
    {
        int e2f(int nf, int e) { return e % nf; }; 
        int fi2e(int nf, int f, int i) { return nf * i + f; }; 
        int e2fi(int nf, int e) { return e / nf; };

        bool TestEdgeCase(
            const std::vector<std::vector<int> > &uE2Fsorted, 
            const Eigen::MatrixXi &uE,
            const Eigen::MatrixXi &F,
            const Snaxl::EdgeSnaxl &si,
            const Snaxl &snaxl,
            int &fj)
        {
            fj = -1; 
            if(snaxl.GetSnaxlType() == Snaxl::Type::Edge)
            {
                auto sj = snaxl.AsEdgeSnaxl(); 
                if(si.ue == sj.ue)
                {
                    return false;
                }

                std::vector<int> intersections; 
                std::set_intersection(
                    uE2Fsorted[si.ue].begin(), uE2Fsorted[si.ue].end(),
                    uE2Fsorted[sj.ue].begin(), uE2Fsorted[sj.ue].end(),
                    std::back_inserter(intersections)
                );

                if(intersections.size() == 1)
                {
                    fj = intersections.front();
                }
            }
            else if(snaxl.GetSnaxlType() == Snaxl::Type::Vertex)
            {
                auto sj = snaxl.AsVertexSnaxl();
                if(uE(si.ue, 0) == sj.v || uE(si.ue, 1) == sj.v)
                {
                    return false; 
                }

                for(auto f: uE2Fsorted[si.ue])
                {
                    if(F(f, 0) == sj.v || F(f, 1) == sj.v || F(f, 2) == sj.v)
                    {
                        fj = f;
                        break;
                    }
                }
            }
            return true; 
        }

        bool TestVertexCase(
            const Eigen::MatrixXi &F,
            const Eigen::MatrixXi &uE, 
            const std::vector<std::vector<int> > &uE2E,
            const std::vector<std::tuple<int, int, bool>> &vhe, 
            const Snaxl::VertexSnaxl &si,
            const Snaxl &snaxl,
            int &vj)
        {
            const auto nf = F.rows(); 
            vj = -1; 
            if(snaxl.GetSnaxlType() == Snaxl::Type::Edge)
            {
                auto sj = snaxl.AsEdgeSnaxl(); 
                if(uE(sj.ue, 0) == si.v || uE(sj.ue, 1) == si.v)
                {
                    vj = uE(sj.ue, 0) == si.v ? uE(sj.ue, 1) : uE(sj.ue, 0); 
                }
                else 
                {
                    for(auto e: uE2E[sj.ue])
                    {
                        auto f = snaking_utils::e2f(nf, e); 
                        auto i = snaking_utils::e2fi(nf, e);
                        if(F(f, i) == si.v)
                        {
                            vj = F(f, (i + 1) % 3);
                            return true; 
                        }
                    }
                    return false; 
                }
            }
            else if(snaxl.GetSnaxlType() == Snaxl::Type::Vertex)
            {
                // TODO: should process duplicate vertices
                auto sj = snaxl.AsVertexSnaxl(); 
                for(auto he: vhe)
                {
                    int f, i, j, k;
                    decode_halfedge(he, f, i, j, k);
                    if(sj.v == F(f, (std::get<2>(he) ? j: k)))
                    {
                        vj = sj.v;
                        return true; 
                    }
                }
                return false; 
            }
            return true; 
        }

        bool Projector(
            const Eigen::MatrixXd &V,
            const Eigen::MatrixXi &uE, 
            const std::vector<std::vector<int> > &uE2E,
            const std::vector<int> &ues,
            const Eigen::Vector3d &v0proj,
            int &common_v,
            std::vector<std::pair<int, int> > &vv, 
            std::map<int, Eigen::Vector3d> &vproj)
        {
            if(ues.size() < 3)
            {
                return false; 
            }
            vv.clear();
            vv.reserve(ues.size()); 

            std::unordered_map<int, int> vcount;
            for(auto ue: ues)
            {
                vcount[uE(ue, 0)]++;
                vcount[uE(ue, 1)]++;
            }

            common_v = -1;
            int max_count = -1;
            for (const auto& pair : vcount) 
            {
                if (pair.second > max_count)
                {
                    max_count = pair.second;
                    common_v = pair.first;
                }
            }

            if(max_count != ues.size() || common_v == -1)
            {
                return false; 
            }

            for(auto ue: ues)
            {
                if(uE(ue, 0) == common_v)
                {
                    vv.emplace_back(uE(ue, 1), ue); 
                }
                else if(uE(ue, 1) == common_v)
                {
                    vv.emplace_back(uE(ue, 0), ue); 
                }
                else 
                {
                    return false;
                }
            }
            
            vproj[common_v] = v0proj; 
            Eigen::Vector3d v0 = V.row(common_v); 
            Eigen::Vector3d v1 = V.row(vv[0].first);
            Eigen::Vector3d e0 = v1 - v0;
            Eigen::Vector3d v1proj = v0proj; 
            v1proj(0) = e0.norm() + vproj[common_v](0);
            v1proj(1) = vproj[common_v](1);
            v1proj(2) = vproj[common_v](2);
            vproj[vv[0].first] = v1proj; 
            e0.normalize();

            double angle = 0.;
            for(int i = 1; i < vv.size(); ++i)
            {
                v1 = V.row(vv[i].first);
                Eigen::Vector3d e1 = v1 - v0; 
                double dot = std::max(-1., std::min(1., e0.dot(e1.normalized()))); 
                angle += std::acos(dot);
                double e1len = e1.norm();
                v1proj(0) = vproj[common_v](0) + e1len * std::cos(angle);
                v1proj(1) = vproj[common_v](1) + e1len * std::sin(angle);
                v1proj(2) = vproj[common_v](2);
                vproj[vv[i].first] = v1proj;
                e0 = e1.normalized();
            }
            return true; 
        }
    
        bool SnaxlProj(
            const std::map<int, Eigen::Vector3d> &vproj,
            const Snaxl &snaxl,
            Eigen::Vector3d &pt)
        {
            if(snaxl.GetSnaxlType() == Snaxl::Type::Edge)
            {
                auto edge_snaxl = snaxl.AsEdgeSnaxl();
                auto t = edge_snaxl.t; 
                auto v0 = edge_snaxl.v0;
                auto v1 = edge_snaxl.v1; 
                if(vproj.find(v0) == vproj.end() || vproj.find(v1) == vproj.end())
                {
                    return false; 
                }
                pt = (1. - t) * vproj.at(v0) + t * vproj.at(v1); 
            }
            else if(snaxl.GetSnaxlType() == Snaxl::Type::Vertex)
            {
                auto v = snaxl.AsVertexSnaxl().v;
                if(vproj.find(v) == vproj.end())
                {
                    return false; 
                }
                pt = vproj.at(v);
            }
            return true; 
        }

        void optimizer(
            const Eigen::Vector3d &prev_pt, 
            const Eigen::Vector3d &next_pt, 
            int common_v,
            const std::vector<std::pair<int, int> > &vv, 
            const std::map<int, Eigen::Vector3d> &vproj,
            std::vector<Snaxl> &snaxls,
            double &dist)
        {
            dist = 0.;
            Eigen::Vector3d v0 = vproj.at(common_v); 
            Eigen::Vector3d curr_prev_pt = prev_pt; 
            std::vector<double> ts;
            bool try_avg = false;
            for(int i = 1; i + 1 < vv.size(); ++i)
            {
                // hit = (1. - t1) * v0 + t1 * v1; 
                // hit = (1. - t2) * prev_pt + t2 * next_pt;
                Eigen::Vector3d v1 = vproj.at(vv[i].first); 
                double t1, t2;
                Eigen::Vector3d hit; 
                if(!intersection2d(v0, v1, prev_pt, next_pt, hit, t1, t2))
                {
                    ts.push_back(0.5);
                    try_avg = true; 
                }
                else 
                {
                    ts.push_back(t1);
                }
                dist += (hit - curr_prev_pt).norm();
                curr_prev_pt = hit; 
                snaxls.emplace_back(Snaxl::EdgeSnaxl(vv[i].second, t1, common_v, vv[i].first)); 
            }
            dist += (curr_prev_pt - next_pt).norm();

            if(!try_avg)
            {
                return;
            }

            snaxls.clear();
            dist = std::numeric_limits<double>::max();

            double avg_dist = 0.;
            curr_prev_pt = prev_pt;
            for(int i = 1; i + 1 < vv.size(); ++i)
            {
                auto &t = ts[i - 1];
                Eigen::Vector3d v1 = vproj.at(vv[i].first); 
                Eigen::Vector3d hit = (1. - t) * v0 + t * v1;
                avg_dist += (hit - curr_prev_pt).norm();
                curr_prev_pt = hit; 
                snaxls.emplace_back(Snaxl::EdgeSnaxl(vv[i].second, t, common_v, vv[i].first)); 
            }
            avg_dist += (curr_prev_pt - next_pt).norm();

            double dist_ori = (v0 - prev_pt).norm() + (v0 - next_pt).norm();

            if(avg_dist > dist_ori)
            {
                snaxls.clear();
                dist = std::numeric_limits<double>::max();
                return;
            }
            dist = avg_dist;
        }
    }

    DGP_INLINE void handle_critical_vertices(
        std::list<Snaxl> &snaxls)
    {
        auto curr = snaxls.begin();
        while(curr != snaxls.end())
        {
            if(curr->GetSnaxlType() == Snaxl::Type::Edge)
            {
                auto begin = curr;
                while(curr->GetSnaxlType() == Snaxl::Type::Edge) 
                {
                    curr++; 
                    if(curr == snaxls.end())
                    {
                        break;
                    }
                }

                begin++; 
                std::unordered_map<int, int> vcount;
                auto it = begin; 
                while(it != curr)
                {
                    auto s = it->AsEdgeSnaxl();
                    vcount[s.v0]++;
                    vcount[s.v1]++;
                    const auto max_count = std::distance(begin, it);
                    if(vcount[s.v0] != max_count + 1 && vcount[s.v1] != max_count + 1)
                    {
                        if(max_count < 3) 
                        {
                            vcount.clear(); 
                            begin = it;
                            continue;
                        }

                        auto prev = std::prev(begin)->AsEdgeSnaxl();
                        if(vcount.count(prev.v0) == max_count || vcount.count(prev.v1) == max_count)
                        {
                            vcount.clear(); 
                            begin = it;
                            continue;
                        }

                        int v = -1; 
                        for(auto &pair: vcount)
                        {   
                            if(pair.second == max_count)
                            {
                                if(v != -1)
                                {
                                    vcount.clear(); 
                                    begin = it;
                                    continue;
                                }
                                v = pair.first; 
                            }
                        }

                        if(v == -1)
                        {
                            vcount.clear(); 
                            begin = it;
                            continue; 
                        }

                        it = snaxls.erase(begin, it);
                        it = snaxls.insert(it, dgp::Snaxl::VertexSnaxl(v));
                        std::advance(it, 1);
                        begin = it; 
                        vcount.clear();
                    }
                    else 
                    {
                        ++it; 
                    }
                }
            }
            else
            {
                ++curr;
            }
        }
    }

    DGP_INLINE bool snaking(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const SnakingData &snaking_data,
        std::list<Snaxl> &snaxls)
    {
        const auto nf = F.rows(); 
        const auto &E = snaking_data.E;
        const auto &uE = snaking_data.uE;
        const auto &EMAP = snaking_data.EMAP;
        const auto &uE2E = snaking_data.uE2E; 
        const auto &uE2Fsorted = snaking_data.uE2Fsorted; 
        const auto &VHe = snaking_data.VHe; 
        const auto &HeOpp = snaking_data.HeOpp; 

        const auto countVertex = [](
            const Snaxl &snaxl,
            std::unordered_map<int,int> &vcount)
        {
            if(snaxl.GetSnaxlType() == Snaxl::Type::Edge)
            {
                auto s = snaxl.AsEdgeSnaxl();
                vcount[s.v0]++;
                vcount[s.v1]++;
            }
            else if (snaxl.GetSnaxlType() == Snaxl::Type::Vertex)
            {
                auto s = snaxl.AsVertexSnaxl(); 
                vcount[s.v]++;
            }
        };

        const auto commonCode = [&](
            const Snaxl &prev,
            const Snaxl &next, 
            const std::vector<int> &ues,
            std::vector<Snaxl> &local_snaxls,
            double &dist)->bool
        {
            int common_v;
            std::vector<std::pair<int, int> > vv; 
            std::map<int, Eigen::Vector3d> vproj; 
            if(!snaking_utils::Projector(V, uE, uE2E, ues, Eigen::Vector3d::Zero(), common_v, vv, vproj))
            {
                return false; 
            }

            Eigen::Vector3d prev_pt_proj; 
            if(!snaking_utils::SnaxlProj(vproj, prev, prev_pt_proj))
            {
                return false; 
            }

            Eigen::Vector3d next_pt_proj; 
            if(!snaking_utils::SnaxlProj(vproj, next, next_pt_proj))
            {
                return false; 
            }
            
            snaking_utils::optimizer(prev_pt_proj, next_pt_proj, common_v, vv, vproj, local_snaxls, dist);
            return true; 
        };

        auto curr = ++snaxls.begin();
        while(std::next(curr) != snaxls.end())
        {
            auto prev = std::prev(curr);
            auto next = std::next(curr);

            if(curr->GetSnaxlType() == Snaxl::Type::Edge)
            {
                auto curr_snaxl = curr->AsEdgeSnaxl();

                int fp;
                if(!snaking_utils::TestEdgeCase(uE2Fsorted, uE, F, curr_snaxl, *prev, fp))
                {
                    curr = snaxls.erase(curr); 
                    continue;    
                }

                int fn;
                if(!snaking_utils::TestEdgeCase(uE2Fsorted, uE, F,curr_snaxl, *next, fn))
                {
                    curr = snaxls.erase(curr); 
                    continue; 
                }

                if(fp == -1 || fn == -1)
                {
                    return false; 
                }

                if(uE2E[curr_snaxl.ue].size() != 2)
                {
                    return false; 
                }

                std::vector<int> ues;
                auto e0 = uE2E[curr_snaxl.ue].front(); 
                auto e0f = snaking_utils::e2f(nf, e0); 
                auto e0i = snaking_utils::e2fi(nf, e0);
                auto e0v = F(e0f, (e0i+1)%3); 
                ues.push_back(EMAP(snaking_utils::fi2e(nf, e0f, (e0i+1)%3)));
                
                ues.push_back(curr_snaxl.ue);

                auto e1 = uE2E[curr_snaxl.ue].back(); 
                auto e1f = snaking_utils::e2f(nf, e1); 
                auto e1i = snaking_utils::e2fi(nf, e1);
                auto e1v = F(e1f, (e1i+2)%3);
                if(e0v != e1v)
                {
                    return false; 
                } 
                ues.push_back(EMAP(snaking_utils::fi2e(nf, e1f, (e1i+2)%3)));
                
                double dist; 
                std::vector<Snaxl> local_snaxls; 
                if(!commonCode(*prev, *next, ues, local_snaxls, dist))
                {
                    return false;
                }

                if(local_snaxls.size() == 1)
                {
                    *curr = local_snaxls.front(); 
                    ++curr;
                    continue;
                }
            }
            else if(curr->GetSnaxlType() == Snaxl::Type::Vertex)
            {
                auto curr_snaxl = curr->AsVertexSnaxl();
                const auto &is_circle = VHe[curr_snaxl.v].first;
                const auto &vhe = VHe[curr_snaxl.v].second;

                int vp; 
                if(!snaking_utils::TestVertexCase(F, uE, uE2E, vhe, curr_snaxl, *prev, vp))
                {
                    return false; 
                }

                int vn; 
                if(!snaking_utils::TestVertexCase(F, uE, uE2E, vhe, curr_snaxl, *next, vn))
                {
                    return false; 
                }

                std::vector<int> hepis; 
                std::vector<int> henis;

                std::unordered_map<int, int> prev_vcount; 
                countVertex(*prev, prev_vcount);
                prev_vcount[curr_snaxl.v] = 0;
                std::unordered_map<int, int> next_vcount; 
                countVertex(*next, next_vcount);
                next_vcount[curr_snaxl.v] = 0;
                for(int index = 0; index < vhe.size(); ++index)
                {
                    int f, i, j, k;
                    decode_halfedge(vhe[index], f, i, j, k);
                    auto v = F(f, std::get<2>(vhe[index])?j:k);
                    if(prev_vcount.count(v) > 0)
                    {
                        hepis.push_back(index); 
                    }
                    if(next_vcount.count(v) > 0)
                    {
                        henis.push_back(index); 
                    }
                }

                auto cmp = [](const auto& a, const auto& b) { return a.first > b.first; };
                const auto n = vhe.size();
                std::priority_queue<
                    std::pair<double, std::vector<Snaxl>>,          
                    std::vector<std::pair<double, std::vector<Snaxl>>>, 
                    decltype(cmp)> min_heap(cmp);

                for(auto hepi: hepis)
                {
                    for(auto heni: henis)
                    {
                        if(is_circle)
                        {
                            std::vector<int> b2e; 
                            for(auto begin = hepi, end = heni; begin != end; begin = (begin + 1) % n)
                            {
                                int f, i, j, k;
                                decode_halfedge(vhe[begin], f, i, j, k);
                                b2e.push_back(EMAP(snaking_utils::fi2e(nf, f, std::get<2>(vhe[begin])?k:j)));
                            }
                            int f, i, j, k;
                            decode_halfedge(vhe[heni], f, i, j, k);
                            b2e.push_back(EMAP(snaking_utils::fi2e(nf, f, std::get<2>(vhe[heni])?k:j)));
        
                            std::vector<Snaxl> s0; 
                            double s0_sq_dist; 
                            if(commonCode(*prev, *next, b2e, s0, s0_sq_dist))
                            {
                                min_heap.push(std::make_pair(s0_sq_dist, s0));
                            }

                            std::vector<int> e2b; 
                            for(auto begin = heni, end = hepi; begin != end; begin = (begin + 1) % n)
                            {
                                int f, i, j, k;
                                decode_halfedge(vhe[begin], f, i, j, k);
                                e2b.push_back(EMAP(snaking_utils::fi2e(nf, f, std::get<2>(vhe[begin])?k:j)));
                            }
                            decode_halfedge(vhe[hepi], f, i, j, k);
                            e2b.push_back(EMAP(snaking_utils::fi2e(nf, f, std::get<2>(vhe[hepi])?k:j)));
                            std::reverse(e2b.begin(), e2b.end());
        
                            std::vector<Snaxl> s1; 
                            double s1_sq_dist;
                            if(commonCode(*prev, *next, e2b, s1, s1_sq_dist))
                            {
                                min_heap.push(std::make_pair(s1_sq_dist, s1));
                            }
                        }
                        else 
                        {  
                            std::vector<int> ues;
                            auto begin = std::min(hepi, heni);
                            auto end = std::max(hepi, heni);
                            ues.reserve(end + 1 - begin ); 
                            for(; begin <= end; ++begin)
                            {
                                int f, i, j, k;
                                decode_halfedge(vhe[begin], f, i, j, k);
                                ues.push_back(EMAP(snaking_utils::fi2e(nf, f, std::get<2>(vhe[begin])?k:j)));
                            }
                            if(hepi > heni)
                            {
                                std::reverse(ues.begin(), ues.end());
                            }

                            double dist;
                            std::vector<Snaxl> s; 
                            if(commonCode(*prev, *next, ues, s, dist))
                            {
                                min_heap.push(std::make_pair(dist, s));
                            }
                        }
                    }
                }

                if(!min_heap.empty())
                {
                    auto min_elem = min_heap.top(); 
                    if(!min_elem.second.empty())
                    {
                        curr = snaxls.erase(curr);
                        curr = snaxls.insert(curr, min_elem.second.begin(), min_elem.second.end());
                        std::advance(curr, min_elem.second.size());
                        continue;
                    }
                }
            }
            ++curr;
        }
        return true;
    }
}

#ifdef DGP_STATIC_LIBRARY

#endif
