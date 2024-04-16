CXX = g++ # /usr/local/bin/clang++ # icpx # g++
CXXFLAGS = -Wall -O3 -std=c++17 -fopenmp # -g
INCLUDES = -I/usr/lib/gcc/x86_64-linux-gnu/11/include/ # -I/usr/lib/gcc/x86_64-linux-gnu/11/
LDFLAGS = -L/lib/gcc/x86_64-linux-gnu/11/

DIR = ./

# # all: solve.o solve.o eps.o tt.o

all: wolff.o
# all: tt.o eps.o t_vev.o psipsi.o eig.o

wolff.o: wolff_hex.cc header.hpp
	$(CXX) $< $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -o $(DIR)$@
