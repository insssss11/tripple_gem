# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(tgem)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(tgem tgem.cpp)

target_link_libraries(tgem Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file
  tgemcell.geo tgemcell.sif tgemcell/dielectrics.dat)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
