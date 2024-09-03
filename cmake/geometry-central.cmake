include(FetchContent)

FetchContent_Declare(
  geometry-central
  GIT_REPOSITORY https://github.com/nmwsharp/geometry-central.git
  GIT_TAG master
) 

FetchContent_MakeAvailable(geometry-central)