#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <igl/read_triangle_mesh.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/pathinfo.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/jet.h>

#include <cmath>
#include <queue>
#include <iostream>
#include <sstream>
#include <fstream>

void dual_graph(const Eigen::MatrixXi &F, Eigen::MatrixXi &E)
{
    Eigen::MatrixXi TT, TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);
    const auto nf = F.rows();

    std::vector<int> edges; edges.reserve(2 * F.maxCoeff() * 3); // e = 3 * v - 3 - b
    std::vector<bool> visited(nf, false);
    for(int fid = 0; fid < nf; ++ fid)
    {
        if(visited[fid])
        {
            continue;
        }

        std::queue<int> queue;
        queue.push(fid);
        while (!queue.empty())
        {
            const auto fid = queue.front();
            queue.pop();
            if(visited[fid])
            {
                continue;
            }
            
            for(int i = 0; i < 3; ++i)
            {
                const auto adjacent_fid = TT(fid, i);
                if(adjacent_fid == -1 || visited[fid])
                {
                    continue;
                }
                edges.push_back(fid);
                edges.push_back(adjacent_fid);
                queue.push(adjacent_fid);
            }
            visited[fid] = true;
        }
    }
    E = Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor>>(edges.data(), edges.size() / 2, 2);
}

void calculate_costs(
    const Eigen::MatrixXd &V, 
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXd &BC,
    const Eigen::MatrixXd &N, 
    const Eigen::MatrixXi &E, 
    Eigen::VectorXd &W)
{
    const auto calculateGeodesicDistance = [&](int fid, int adjacent_fid)->double
    {
        int edge_start = -1; 
        int edge_end = -1; 
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                if(F(fid, i) == F(adjacent_fid, j))
                {
                    if(F(fid, (i + 1) % 3) == F(adjacent_fid, (j + 3 - 1) % 3))
                    {
                        edge_start = F(fid, i); 
                        edge_end = F(fid, (i + 1) % 3);
                    }
                    else if(F(fid, (i + 3 - 1) % 3) == F(adjacent_fid, (j + 1) % 3))
                    {
                        edge_end = F(fid, i);
                        edge_start = F(fid, (i + 3 - 1) % 3);
                    }
                    else 
                    {
                        assert(false && "The two faces are not adjacent!");
                    }
                }
            }
        }
        const Eigen::RowVector3d v0 = V.row(edge_start); 
        const Eigen::RowVector3d v1 = V.row(edge_end);
        const Eigen::RowVector3d e = v1 - v0;
        const double edge_length = e.norm();
        double distance = 0.0;
        {
            const Eigen::RowVector3d bc = BC.row(fid);
            const Eigen::RowVector3d v = bc - v0;
            distance += v.cross(e).norm() / edge_length;
        }
        {
            const Eigen::RowVector3d bc = BC.row(adjacent_fid);
            const Eigen::RowVector3d v = bc - v0;
            distance += v.cross(e).norm() / edge_length;
        }
        return distance;
    };

    const double pi = 3.14159265358979323846;
    const double rate = 1000.0; 
    const double lambdas = 100.0;

    const auto ne = E.rows();
    W.setZero(ne);
    for(int eid = 0; eid < ne; ++eid)
    {
        const auto fid = E(eid, 0);
        const auto adjacent_fid = E(eid, 1);

        const double dot = std::max(-1.0, std::min(1.0, N.row(fid).dot(N.row(adjacent_fid))));
        const double dihedral = std::acos(dot) + 1e-8; 
        const double geodesic = calculateGeodesicDistance(fid, adjacent_fid);

        double cost = rate * lambdas * (-std::log(dihedral / pi) * geodesic);
        const Eigen::RowVector3d vc = (BC.row(adjacent_fid) - BC.row(fid)).normalized();
        vc.dot(N.row(fid)) < 0 ? cost *= (1.0 + std::abs(dot)) : cost *= 1.0;
        W(eid) = cost;
    }
}

void add_graph(
    const Eigen::MatrixXi &E,
    const Eigen::MatrixXd &BC,
    const Eigen::MatrixXd &C,
    const Eigen::VectorXd &W, 
    igl::opengl::glfw::Viewer &viewer)
{
    viewer.data().clear_edges(); 
    viewer.data().clear_points(); 
    viewer.data().clear_labels();

    const auto ne = E.rows();
    for(int eid = 0; eid < ne; ++eid)
    {
        viewer.data().add_edges(
            BC.row(E(eid, 0)),
            BC.row(E(eid, 1)),
            C.row(eid)
        );

        std::stringstream ss;
        const Eigen::RowVector3d bc = 0.5 * (BC.row(E(eid, 0)) + BC.row(E(eid, 1)));
        ss << W(eid) << std::endl;
        viewer.data().add_label(bc, ss.str());
    }
    viewer.data().add_points(BC, Eigen::RowVector3d(0.0, 1.0, 0.0));
}

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V; 
    Eigen::MatrixXi F;
    if(!igl::read_triangle_mesh(argv[1], V, F))
    {
        return 0;
    }
    std::string dir,base,ext,name;
    igl::pathinfo(argv[1],dir,base,ext,name);
    if(ext=="stl")
    {
        Eigen::MatrixXd temp_V;
        Eigen::VectorXi SVI,SVJ;
        igl::remove_duplicate_vertices(temp_V,0,V,SVI,SVJ);
        std::for_each(F.data(),F.data()+F.size(),[&SVJ](int & f){f=SVJ(f);});   
    }

    Eigen::MatrixXd BC;
    igl::barycenter(V, F, BC);

    Eigen::MatrixXd N;
    igl::per_face_normals(V, F, N);

    Eigen::MatrixXi E; 
    dual_graph(F, E);

    Eigen::VectorXd W;
    calculate_costs(V, F, BC, N, E, W);

    Eigen::MatrixXd C; 
    igl::jet(W, true, C);


    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);    

    viewer.data().set_mesh(V, F);

    bool show_dual_graph = false;

    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
    {
        std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
        if (key == 'W')
        {
            viewer.data().show_custom_labels = !viewer.data().show_custom_labels;
            return true;
        }

        if(key == 'G')
        {
            show_dual_graph = !show_dual_graph;
            if(show_dual_graph)
            {
                add_graph(E, BC, C, W, viewer);
            }
            else
            {
                viewer.data().clear_edges();
                viewer.data().clear_points();
                viewer.data().clear_labels();
            }
            return true;
        }

        return false;
    };

    double point_size = 5.f;
    double line_width = 0.5f * point_size;
    double label_size = 0.5f * point_size;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &v) -> bool
    {
        viewer.data().point_size = point_size;
        viewer.data().line_width = line_width;
        viewer.data().label_size = label_size;
        return false;
    };

    menu.callback_draw_custom_window = [&]()
    {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiCond_FirstUseEver);
        ImGui::Begin(
            "Dual-Graph Parameters", nullptr,
            ImGuiWindowFlags_NoSavedSettings
        );

        ImGui::PushItemWidth(-80);
        ImGui::DragScalar("Point Size", ImGuiDataType_Double, &point_size, 0.5, 0, 0, "%.4f");
        ImGui::PopItemWidth();

        ImGui::PushItemWidth(-80);
        ImGui::DragScalar("Line Width", ImGuiDataType_Double, &line_width, 0.5, 0, 0, "%.4f");
        ImGui::PopItemWidth();

        ImGui::PushItemWidth(-80);
        ImGui::DragScalar("Label Size", ImGuiDataType_Double, &label_size, 0.5, 0, 0, "%.4f");
        ImGui::PopItemWidth();

        ImGui::InputText("Name", name);

        ImGui::End();
    };

    viewer.launch();
    return 0; 
}
