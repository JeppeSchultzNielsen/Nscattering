cmake_minimum_required(VERSION 3.21)
project(aluscattering)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(match)
include_directories(setup)

find_package(AUSALIB REQUIRED)

add_executable(Analyzer Analyse.cpp)
target_link_libraries(Analyzer PRIVATE AUSA::AUSA)
add_executable(Vinkler Vinkelhistogrammer.cpp)
target_link_libraries(Vinkler PRIVATE AUSA::AUSA)
add_executable(Run run.cpp Vinkelhistogrammer.cpp Analyse.cpp Thickness.cpp FindCurrent.cpp include/runner.h include/runner2.h)
target_link_libraries(Run PRIVATE AUSA::AUSA)
add_executable(Thickness Thickness.cpp)
target_link_libraries(Thickness PRIVATE AUSA::AUSA)