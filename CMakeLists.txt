cmake_minimum_required(VERSION 3.0)
project(MPCLIB)
find_package(OsqpEigen REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(include)

add_executable(DoubleIntegratorExample examples/doubleIntegratorExample.cpp src/MpcController.cpp)
add_executable(AckermannExample examples/ackermannExample.cpp src/MpcController.cpp)

target_link_libraries(DoubleIntegratorExample OsqpEigen::OsqpEigen Eigen3::Eigen)
target_link_libraries(AckermannExample OsqpEigen::OsqpEigen Eigen3::Eigen)