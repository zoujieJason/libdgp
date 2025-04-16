
// #include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
// #include <igl/opengl/glfw/imgui/ImGuiMenu.h>
// #include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
// #include <igl/opengl/glfw/imgui/SelectionWidget.h>

// #include <igl/remove_duplicate_vertices.h>

// #include <dgp/snaking.h>

// const Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);
// const Eigen::RowVector3d red(igl::MAYA_RED(0), igl::MAYA_RED(1), igl::MAYA_RED(2));

// int main(int argc, char **argv)
// {
//     igl::opengl::glfw::Viewer viewer;
//     igl::opengl::glfw::imgui::ImGuiPlugin plugin;
//     viewer.plugins.push_back(&plugin);
//     igl::opengl::glfw::imgui::ImGuiMenu menu;
//     plugin.widgets.push_back(&menu);

//     const Eigen::MatrixXd V = (Eigen::MatrixXd(6, 3)
//     << 0.0, 0.0, 0.0,
//       1.0, -0.5, 0.0,
//       1.0, 1.0, 0.0, 
//       2.0, 0.0, 0.0,
//       2.5, 1.0, 0.0,
//       1.0, 2.0, 0.0).finished();
//     const Eigen::MatrixXi F = (Eigen::MatrixXi(5, 3)
//     << 0, 1, 2, 1, 3, 2, 2, 3, 4, 2, 4, 5, 5, 0, 2).finished();

//     // const Eigen::MatrixXd V = (Eigen::MatrixXd(5, 3)
//     // << 0.0, 0.0, 0.0,
//     //   1.0, -0.5, 0.0,
//     //   1.0, 1.0, 0.0, 
//     //   2.0, 0.0, 0.0,
//     //   2.5, 1.0, 0.0).finished();
//     // const Eigen::MatrixXi F = (Eigen::MatrixXi(3, 3)
//     // << 0, 1, 2, 1, 3, 2, 2, 3, 4).finished();

//     dgp::SnakingData snaking_data; 
//     dgp::snaking_precompute(F, snaking_data); 
//     std::list<dgp::Snaxl> snaxls;
//     // snaxls.emplace_back(dgp::Snaxl::VertexSnaxl(0));
//     // snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(0), 0.1, 1, 2));
//     // snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(0), 0.2, 1, 2));
//     // snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(1), 0.1, 3, 2));
//     // snaxls.emplace_back(dgp::Snaxl::VertexSnaxl(4));
//     // snaxls.emplace_back(dgp::Snaxl::VertexSnaxl(0));
//     // snaxls.emplace_back(dgp::Snaxl::VertexSnaxl(2));
//     // snaxls.emplace_back(dgp::Snaxl::VertexSnaxl(4));

//     std::cout << snaking_data.uE  << "\n";

//     std::cout << snaking_data.EMAP.transpose() << "\n";

//     snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(14), 0.1, 5, 0));
//     snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(4), 0.8, 2, 0));
//     snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(0), 0.1, 1, 2));
//     snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(1), 0.1, 3, 2));
//     snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(2), 0.5, 3, 4));

//     // snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(6), 0.5, 0, 1));
//     // snaxls.emplace_back(dgp::Snaxl::VertexSnaxl(2));
//     // snaxls.emplace_back(dgp::Snaxl::EdgeSnaxl(snaking_data.EMAP(2), 0.5, 3, 4));

//     Eigen::MatrixXd iV;
//     Eigen::MatrixXi iE; 
//     dgp::extract_snaxls(V, snaxls.begin(), snaxls.end(), iV, iE); 

//     viewer.data().set_mesh(V, F);
//     viewer.data().set_edges(iV, iE, sea_green);
//     viewer.data().set_points(iV, red);

//     viewer.callback_key_pressed = [&](
//       igl::opengl::glfw::Viewer& viewer,
//       unsigned char key,
//       int modifier) -> bool
//   {
//       switch (key)
//       {
//           case '0':
//           {
//               viewer.data().clear();
//               viewer.data().set_mesh(V,F);
//               viewer.core().align_camera_center(V, F);

//               dgp::snaking(V, F, snaking_data, snaxls);
//               dgp::handle_critical_vertices(snaxls);
//               dgp::extract_snaxls(V, snaxls.begin(), snaxls.end(), iV, iE); 
//               viewer.data().set_edges(iV, iE, sea_green);
//               viewer.data().set_points(iV, red);
//               return true;
//           }
//           default: return false; 
//       }
//       return true; 
//   };
//     viewer.launch();
//     return 0;
// }

// #include <dgp/halfedge_opposite.h>
// #include <dgp/vertex_halfedge_adjacency.h>

// #include <igl/oriented_facets.h>
// #include <igl/orient_halfedges.h>
// #include <igl/vertex_triangle_adjacency.h>

// int main(int argc, char **argv)
// {
//     const Eigen::MatrixXi F = (Eigen::MatrixXi(3, 3)
//     << 0, 1, 3, 1, 2, 3, 3, 2, 0).finished();

//     Eigen::MatrixXi E, oE;
//     igl::orient_halfedges(F, E, oE);
//     std::cout << E << std::endl;
    
//     Eigen::MatrixXi HeOpp;
//     dgp::halfedge_opposite(E, HeOpp);
//     std::cout << "HeOpp: \n";
//     std::cout << HeOpp << "\n";
//     std::cout << "\n";

//     std::vector<std::vector<int> > VF, VFi; 
//     igl::vertex_triangle_adjacency(4, F, VF, VFi);

//     std::vector<std::vector<std::tuple<int, int, bool>>> VHe;
//     if(!dgp::vertex_halfedge_adjacency(VF, VFi, F, HeOpp, VHe)) 
//     {
//         return 1;
//     }

//     for(auto vhe: VHe)
//     {
//         for(auto he: vhe)
//         {
//             std::cout << std::get<0>(he) << " " << std::get<1>(he)  << " " <<  std::get<2>(he) << " | ";
//         }
//         std::cout << std::endl;
//     }

//     return 0;
// }

#include <dgp/read_labels.h>
#include <dgp/remove_unref.h>
#include <dgp/snaking.h>

#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/SelectionWidget.h>

#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/AABB.h>
#include <igl/screen_space_selection.h>
#include <igl/boundary_loop.h>

int iLabel = 0;
Eigen::MatrixXd V; 
Eigen::MatrixXi F; 
Eigen::VectorXi I, J, FMAP;
igl::AABB<Eigen::MatrixXd, 3> tree;
Eigen::VectorXd W; 
Eigen::Array<double,Eigen::Dynamic,1> and_visible;
Eigen::VectorXd S; 
std::vector<size_t> selected_indices;
std::vector<int> boundary; 
std::list<dgp::Snaxl> snaxls; 
dgp::SnakingData snaking_data; 

Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);
Eigen::RowVector3d red(igl::MAYA_RED(0), igl::MAYA_RED(1), igl::MAYA_RED(2));
Eigen::RowVector3d grey(igl::MAYA_GREY(0), igl::MAYA_GREY(1), igl::MAYA_GREY(2));
Eigen::RowVector3d green(igl::MAYA_GREEN(0), igl::MAYA_GREEN(1), igl::MAYA_GREEN(2));
bool key_swith = true; 

int main(int argc, char **argv)
{
    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);
    igl::opengl::glfw::imgui::SelectionWidget widget;
    plugin.widgets.push_back(&widget);
    widget.mode = igl::opengl::glfw::imgui::SelectionWidget::Mode::OFF; 

    Eigen::MatrixXd inputV; 
    Eigen::MatrixXi inputF;
    std::string meshname(argv[1]);
    if(!igl::read_triangle_mesh(meshname,inputV,inputF)) 
    {
      std::cout << "Failed to load mesh " << meshname << std::endl;
      return 0;
    }
    std::string dirname, basename, extension, filename; 
    igl::pathinfo(meshname, dirname, basename, extension, filename);
    if(extension == "stl")
    {
        Eigen::MatrixXd tempV = inputV; 
        Eigen::VectorXi SVI, SVJ; 
        igl::remove_duplicate_vertices(tempV,0,inputV,SVI,SVJ);
        std::for_each(inputF.data(),inputF.data()+inputF.size(),[&SVJ](int & f){f=SVJ(f);});
    }
  
    std::string labelname(argv[2]);
    std::vector<int> vLabel;
    dgp::read_labels(labelname, vLabel);
    if(vLabel.size() != inputF.rows()) 
    {
      std::cout << "Failed to load label " << labelname << std::endl;
      return 0;
    }
    std::vector<std::string> vsLabel; 
    std::set<int> mLabel(vLabel.begin(), vLabel.end());
    std::for_each(mLabel.begin(), mLabel.end(), [&](int label) { if(label != 0) vsLabel.push_back(std::to_string(label)); });

    Eigen::VectorXi input_faces_label = Eigen::Map<Eigen::VectorXi>(vLabel.data(), vLabel.size());

    menu.callback_draw_custom_window = [&]()
    {
        // base
        float window_current_height = 10.f; 
        {
            float window_size_y = 80.f;
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), window_current_height), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, window_size_y), ImGuiCond_FirstUseEver);
            ImGui::Begin("Base", nullptr, ImGuiWindowFlags_NoSavedSettings);
            ImGui::Combo("label", &iLabel, vsLabel);
            if(ImGui::Button("clear"))
            {
                selected_indices.clear();
            }
            ImGui::End();
            window_current_height += window_size_y; 
        }

        // precompute
        {
            float window_size_y = 80.f;
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), window_current_height), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(210, window_size_y), ImGuiCond_FirstUseEver);
            ImGui::Begin("precompute", nullptr, ImGuiWindowFlags_NoSavedSettings);
            if(ImGui::Button("precompute"))
            {
                Eigen::VectorXi Ki = input_faces_label.unaryExpr([&](auto fl) { return fl == std::stoi(vsLabel[iLabel]) ? 1 : 0; });
				dgp::remove_unref(inputV, inputF, Ki, V, F, I, J, FMAP);
                and_visible.setZero(V.rows());
                W.setZero(V.rows());
                S.setZero(V.rows());
                tree.init(V, F);
                dgp::snaking_precompute(F, snaking_data); 
            }
            if(ImGui::Button("build snaxls"))
            {
                std::vector<bool> selected(V.rows(), false); 
                std::for_each(selected_indices.begin(), selected_indices.end(), [&selected](auto v){ selected[v] = true; }); 
                igl::boundary_loop(F, boundary);
                snaxls.clear(); 
                for(int i = 0; i < boundary.size(); ++i)
                {
                    if(selected[boundary[i]])
                    {
                        for(int j = i; j < boundary.size(); ++j)
                        {
                            if(!selected[boundary[j]])
                            {
                                break;
                            }
                            snaxls.emplace_back(dgp::Snaxl::VertexSnaxl(boundary[j]));
                        }
                        break;
                    }
                }
            }
            ImGui::End();
            window_current_height += window_size_y; 
        }
    };

    widget.callback = [&]()
    {
        igl::screen_space_selection(
            V,F,tree,
            viewer.core().view,
            viewer.core().proj,
            viewer.core().viewport,
            widget.L,
            W,
            and_visible);

        S = W;
        S.array() *= and_visible;
        selected_indices.clear();
        for(Eigen::Index i = 0; i < S.size(); ++i)
        {
            if(S(i)>0.)
            {
                selected_indices.push_back(i); 
            }
        }
        W.setZero(V.rows());
        and_visible.setZero(V.rows());
    };

    viewer.data().set_mesh(inputV,inputF);
    viewer.data().show_lines = false;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

    viewer.callback_key_pressed = [&](
        igl::opengl::glfw::Viewer& viewer,
        unsigned char key,
        int modifier) -> bool
    {
        switch (key)
        {
            case '0':
            {
                const Eigen::MatrixXd CM = (Eigen::MatrixXd(2,3)<<
                0.3,0.3,0.5,                 
                255.0/255.0,228.0/255.0,58.0/255.0).finished();

                viewer.data().clear();
                viewer.data().set_mesh(V,F);
                viewer.core().align_camera_center(V, F);
                viewer.data().set_data(S,0,1,igl::COLOR_MAP_TYPE_PLASMA,2);
                viewer.data().set_colormap(CM);
                viewer.data().set_points(V(selected_indices, Eigen::all), red);
                return true;
            }
            case '1':
            {
                Eigen::MatrixXd colors(F.rows(), 3);    
                for(int f = 0; f < F.rows(); ++f)
                {
                    //colors.row(f) = static_cast<bool>(input_faces_label(FMAP(f)))? red : grey;
                    colors.row(f) = grey;
                }
                                       
                Eigen::MatrixXd iV;
                Eigen::MatrixXi iE; 
                dgp::extract_snaxls(V, snaxls.begin(), snaxls.end(), iV, iE); 

                viewer.data().clear();
                viewer.data().set_mesh(V,F);
                viewer.core().align_camera_center(V, F);
                viewer.data().set_colors(colors);
                viewer.data().set_edges(iV, iE, sea_green);
                viewer.data().set_points(iV, red);
                return true;
            }
            case '2':
            {          
                dgp::snaking(V, F, snaking_data, snaxls);
                dgp::handle_critical_vertices(snaxls);
                return true;
            }
            default: return false; 
        }
        return true; 
    };

    viewer.launch();
    return 0;
}