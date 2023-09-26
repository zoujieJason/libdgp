function(dgp_add_tutorial name)
    # Only create tutorial if dependencies have been enabled in the current build
    foreach(module IN ITEMS ${ARGN})
        if(NOT TARGET ${module})
            return()
        endif()
    endforeach()

    message(STATUS "Creating libdgp tutorial: ${name}")
    add_executable(${name} ${CMAKE_CURRENT_SOURCE_DIR}/${name}/main.cpp)
    target_link_libraries(${name} PRIVATE ${ARGN})

    set_target_properties(${name} PROPERTIES FOLDER Libdgp_Tutorials)
endfunction()
