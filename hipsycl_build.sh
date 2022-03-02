mkdir build
cd build
cmake -D CMAKE_CXX_COMPILER=/opt/hipSYCL/bin/syclcc-clang ..
make
cd ..
