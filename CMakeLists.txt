cmake_minimum_required(VERSION 3.0)
project(MPCLIB)
find_package(OsqpEigen REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(include)

add_executable(TestMpcLib src/testMpcClass.cpp src/MpcController.cpp)

target_link_libraries(TestMpcLib OsqpEigen::OsqpEigen Eigen3::Eigen)