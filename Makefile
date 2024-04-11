CXX = /mnt/hdd_baracuda/opt/llvm-project/build/bin/clang++ # icpx # g++
CXXFLAGS = -Wall -O3 -std=c++20 -fopenmp -march=native
INCLUDES = -I/projectnb/qfe/nmatsumo/opt/eigen/

NVCC = nvcc
NVCCFLAGS = -arch=sm_70 -O3 -lcusolver -std=c++17
INCLUDES_CUDA =


DIR = ./


# # all: solve.o solve.o eps.o tt.o

all: wolff.o
# all: tt.o eps.o t_vev.o psipsi.o eig.o

wolff.o: wolff_hex.cc header.hpp
	$(CXX) $< $(CXXFLAGS) $(INCLUDES) -o $(DIR)$@
