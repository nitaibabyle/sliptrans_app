cmake_minimum_required(VERSION 3.5.0)

set(PROJ_NAME "slip_app")


find_package(pinocchio REQUIRED)
find_package(Eigen3 REQUIRED)
message(Eigen: ${EIGEN3_INCLUDE_DIR})


aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} SLIP_TRANS_SRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/IKQP/ IKQP_SRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/IKQP/QP_Pierre QP_SRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/IKQP/QP_Pierre/f2c SUB_QP_F2C_SRC)

set(ZJT_trans_SRC
  ${SLIP_TRANS_SRC}
  ${IKQP_SRC}
  ${QP_SRC}
  ${SUB_QP_F2C_SRC}
  )

  include_directories(
include
  ${EIGEN3_INCLUDE_DIR}
)


# set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_executable(slip_app ${ZJT_trans_SRC})

# add_library(slip_app STATIC "sliptrans.cpp" ${ZJT_trans_SRC})


target_include_directories(slip_app PUBLIC
  ${pinocchio_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/IKQP>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/IKQP/QP_Pierre>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/IKQP/QP_Pierre/f2c>)

target_link_libraries(slip_app
  ${pinocchio_LIBRARIES}
)
# target_link_libraries(main_app Lib_HLQwalkingmain)
target_compile_options(slip_app PUBLIC)
target_compile_features(slip_app PUBLIC cxx_std_20)
