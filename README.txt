ENVIRONMENT SETUP:

osqp:
- in an appropriate directory: 
- git clone https://github.com/osqp/osqp
- cd osqp
- mkdir build
- cd build
- cmake -G "Unix Makefiles" ..
- cmake --build .
- sudo cmake --build . --target install

osqp-eigen:
- in an appropriate directory:
- git clone https://github.com/robotology/osqp-eigen.git
- cd osqp-eigen
- mkdir build
- cd build
- cmake ..
- make
- sudo make install
- nano ~/.bashrc
  - add the line: "export OsqpEigen_DIR=<INSERT_CORRECT_PATH>/osqp-eigen/build"

eigen:
- in an appropriate directory: 
- git clone https://gitlab.com/libeigen/eigen.git
- cd eigen
- sudo cp -r Eigen usr/local/include

TO BUILD:

- mkdir build
- cd build
- cmake ..
- cmake --build .