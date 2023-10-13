ENVIRONMENT SETUP:

osqp-eigen:
- setup and install following these directions: https://github.com/robotology/osqp-eigen/tree/master
- I put everything in a "Libraries" folder, but the "home" directory may be best
- will need sudo permissions for "make install"
- will need to set environment variable:
  - cd ~
  - nano .bashrc
  - add line: export OsqpEigen_DIR=~/Libraries/osqp-eigen/build
    - sub in appropriate directory

eigen:
- in an appropriate directory (probably "home"): git clone https://gitlab.com/libeigen/eigen.git
- cd eigen
- sudo cp Eigen usr/local/include
  - alternatively make sure all those header files are on the include path

osqp:
- following these instructions: https://osqp.org/docs/get_started/sources.html
- in an appropriate directory (probably "home"): git clone https://github.com/osqp/osqp
- make sure to run the "cmake --build . --target install" to get library onto path
  - will need sudo permissions

TO BUILD:
- mkdir build
- cd build
- cmake ..
- cmake --build .