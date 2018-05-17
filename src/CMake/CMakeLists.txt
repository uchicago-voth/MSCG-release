########################################
# CMake build system
# Created by Christoph Junghans
########################################
cmake_minimum_required(VERSION 3.1)
cmake_minimum_required(VERSION 3.1)

project(mscg LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

include(GNUInstallDirs)
find_package(GSL REQUIRED)
find_package(LAPACK REQUIRED)

set(MSCG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

set(SOVERSION 0)
file(GLOB MSCG_LIB_SOURCES ${MSCG_SOURCE_DIR}/*.cpp)
foreach(_APP newfm rangefinder combinefm)
  file(GLOB MSCG_${_APP}_SOURCES ${MSCG_SOURCE_DIR}/${_APP}.cpp)
  list(REMOVE_ITEM MSCG_LIB_SOURCES ${MSCG_${_APP}_SOURCES})
  add_executable(${_APP} ${MSCG_${_APP}_SOURCES})
  target_link_libraries(${_APP} mscg)
  set_target_properties(${_APP} PROPERTIES OUTPUT_NAME ${_APP}_no_gro.x)
  install(TARGETS ${_APP} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()

add_library(mscg ${MSCG_LIB_SOURCES})
set_target_properties(mscg PROPERTIES SOVERSION ${SOVERSION})
target_compile_options(mscg PRIVATE -DDIMENSION=3 -D_exclude_gromacs=1)
target_include_directories(mscg PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(mscg ${GSL_LIBRARIES} ${LAPACK_LIBRARIES})
install(TARGETS mscg LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

file(GLOB MSCG_HEADERS ${MSCG_SOURCE_DIR}/*.h)
install(FILES ${MSCG_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mscg)
