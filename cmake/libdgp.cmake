if(LIBDGP_USE_STATIC_LIBRARY)
    set(DGP_SCOPE PUBLIC)
    add_library(dgp STATIC)
    set_target_properties(dpg PROPERTIES OUTPUT_NAME dpg)
else()
    set(DGP_SCOPE INTERFACE)
    add_library(dgp INTERFACE)
endif()

include(GNUInstallDirs)
target_include_directories(dgp ${DGP_SCOPE}
        $<BUILD_INTERFACE:${libdgp_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )

file(GLOB INC_FILES "${libdgp_SOURCE_DIR}/include/dgp/*.h")
file(GLOB SRC_FILES "${libdgp_SOURCE_DIR}/include/dgp/*.cpp")
file(GLOB INC_FILES "${libdgp_SOURCE_DIR}/include/HLBFGS/*.h")
file(GLOB SRC_FILES "${libdgp_SOURCE_DIR}/include/HLBFGS/*.cpp")
if(${CMAKE_VERSION} VERSION_LESS 3.19 AND NOT LIBDGP_USE_STATIC_LIBRARY)
    # Old approach: defines a dummy custom target for the IDE
    add_custom_target(dgp_ SOURCES ${INC_FILES} ${SRC_FILES})
else()
    target_sources(dgp PRIVATE ${INC_FILES} ${SRC_FILES})
endif()

target_link_libraries(dgp ${DGP_SCOPE} Eigen3::Eigen)
target_compile_definitions(dgp INTERFACE -DLIBDGP_DATA_PATH="${libdgp_SOURCE_DIR}/data")
