cmake_minimum_required(VERSION 3.3)

project(demo_cfftw LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_BUILD_TYPE Release)

find_package(LAPACK)

find_package(VTK 8.2.0 REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

add_executable(demo_cfftw main.cpp)
target_link_libraries(demo_cfftw PRIVATE ${VTK_LIBRARIES} -lfftw3)
target_include_directories(demo_cfftw PRIVATE ../castor-fftw/include/)