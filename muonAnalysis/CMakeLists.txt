# CMakeists.txt for muon diffusion analysis root files

#----------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(muonAnalysis)

#----------------------------------------------------------
# Location of root fies
list(APPEND CMAKE_PREFIX_PATH $ROOTSYS/etc/cmake/FindROOT.cmake})
find_package(ROOT REQUIRED)
set(ROOT_INCLUDE_PUBLIC /gscratch/damic/sw/include)

message(${ROOT_INCLUDE_PUBLIC})
#----------------------------------------------------------
# Include root and project directories
include_directories(${ROOT_INCLUDE_PUBLIC})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------
# Locate source and header files for project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------
# Create an executable and link to root libraries
add_executable(muonanalysis main.cxx ${sources} ${headers})
target_link_libraries(muonanalysis ${ROOT_LIBRARIES})
