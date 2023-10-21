ENVIRONMENT SETUP:

osqp:
- following these instructions: https://osqp.org/docs/get_started/sources.html
- in appropriate directory: git clone https://github.com/osqp/osqp
- make sure to run the "cmake --build . --target install" to get library onto path
  - will need sudo permissions

osqp-eigen:
- clone from here: https://github.com/robotology/osqp-eigen/tree/master
- looks like osqp-eigen is not up-to-date with osqp as of v0.8.1. I needed to make the following changes to the source:
  - PROFILING -> OSQP_ENABLE_PROFILING in:
    - include/OsqpEigen/Constants.hpp
    - src/Settings.cpp
- setup and install following these directions: https://github.com/robotology/osqp-eigen/tree/master
  - need to build and install after making above changes to source
- will need sudo permissions for "make install"
- will need to set environment variable:
  - cd ~
  - nano .bashrc
  - add line: export OsqpEigen_DIR=<INSERT_CORRECT_PATH>/osqp-eigen/build

eigen:
- in an appropriate directory: git clone https://gitlab.com/libeigen/eigen.git
- cd eigen
- sudo cp Eigen usr/local/include
  - alternatively make sure all those header files are on the include path

TO BUILD:
- mkdir build
- cd build
- cmake ..
- cmake --build .