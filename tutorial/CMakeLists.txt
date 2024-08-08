option(LIBDGP_TUTORIALS_GLFW_TEST "Compile libdgp tutorials glfw_test" OFF)
option(LIBDGP_TUTORIALS_BARYCENTRIC_EMBEDDING "Compile libdgp tutorials barycentric_embedding" OFF)
option(LIBDGP_TUTORIALS_BIJECTIVE_PARAMETERIZATION "Compile libdgp tutorials bijective_parameterization" OFF)
option(LIBDGP_TUTORIALS_ARAP "Compile libdgp tutorials arap" ON)
option(LIBDGP_TUTORIALS_LSCM "Compile libdgp tutorials lscm" ON)
option(LIBDGP_TUTORIALS_MESH_TREE "Compile libdgp tutorials mesh_tree" ON)
option(LIBDGP_TUTORIALS_KDTREE "Compile libdgp tutorials KDTree" ON)
option(LIBDGP_TUTORIALS_ANG_DIST "Compile libdgp tutorials ang_dist" ON)

if(NOT TARGET igl::glfw)
    message(WARNING "Tutorial executables depend on OpenGL and glfw.")
endif()

if(LIBDGP_TUTORIALS_GLFW_TEST)
    dgp_add_tutorial(glfw_test igl::glfw)
endif()     

if(LIBDGP_TUTORIALS_BARYCENTRIC_EMBEDDING)
    dgp_add_tutorial(barycentric_embedding igl::core igl::glfw dgp)
endif()

if(LIBDGP_TUTORIALS_BIJECTIVE_PARAMETERIZATION)
    dgp_add_tutorial(bijective_parameterization igl::core igl::glfw dgp)
endif()

if(LIBDGP_TUTORIALS_ARAP)
    dgp_add_tutorial(arap igl::core igl::glfw dgp)
endif()

if(LIBDGP_TUTORIALS_LSCM)
    dgp_add_tutorial(lscm igl::core igl::glfw dgp)
endif()

if(LIBDGP_TUTORIALS_MESH_TREE)
    dgp_add_tutorial(mesh_tree igl::core igl::glfw dgp)
endif()

if(LIBDGP_TUTORIALS_KDTREE)
    dgp_add_tutorial(kdtree igl::core igl::glfw dgp)
endif()

if(LIBDGP_TUTORIALS_ANG_DIST)
    dgp_add_tutorial(ang_dist igl::core igl::glfw igl::imgui dgp)
endif()