cmake -B build
cd ./build
make 
./target

cmake -G Ninja

cmake --build build
./build/target

find_package
CMAKE_MODULE_PATH -> Find<Target>.cmake -> not found 
-> CMAKE_PREFIX_PATH | <Target>_DIR -> <Target>Config.cmake | <Target>_config.cmake

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
*.cmake -> function | FetchContent
function(function_name name) 
endfunction()
${name}->function_name 
${ARGN}->function(function_name, 1, 2, 3, ...)
foreach(argc IN ITEMS ${ARGN})
    message(${argc})
endforeach()

FetchContent->FETCHCONTENT_BASE_DIR:PATH=*/build/_deps
-> build _deps/libigl/CMakeLists.txt -> eigen | glfw | glad | ... -> ${project_name}/CMakeLists.txt

get_directory_property(LIBIGL_PARENT_DIR PARENT_DIRECTORY) -> set some variables on or off

make edit_cache

$<BUILD_INTERFACE:${libdgp_SOURCE_DIR}/include> -> if BUILD_INTERFACE { ... }
$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> -> if INSTALL_INTERFACE { ... }

cmake debug on vscode:
shift+cmd+p -> cmake debug(ctrl+f5)