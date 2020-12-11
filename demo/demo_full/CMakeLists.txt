cmake_minimum_required(VERSION 3.3)

# create a project
project(demo_cfftw LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_BUILD_TYPE Release)

# find VTK
find_package(VTK 8.2.0 REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

# add a new executable 'demo_cfftw' from the file 'main.cpp'
add_executable(demo_cfftw main.cpp)

# link the VTK library and the FFTW3 for double precision
target_link_libraries(demo_cfftw PRIVATE ${VTK_LIBRARIES} -lfftw3)

# add the 'non-standard' include directories, more particularly
# the one for the headers of 'castor-fftw'
target_include_directories(demo_cfftw PRIVATE ../../../castor-fftw/include/)
