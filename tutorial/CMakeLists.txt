option(LIBDGP_TUTORIALS_GLFW_TEST "Compile libdgp tutorials glfw_test" ON)
option(LIBDGP_TUTORIALS_BARYCENTRIC_EMBEDDING "Compile libdgp tutorials barycentric_embedding" ON)
option(LIBDGP_TUTORIALS_BIJECTIVE_PARAMETERIZATION "Compile libdgp tutorials barycentric_embedding" ON)

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
