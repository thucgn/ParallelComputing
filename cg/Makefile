CXX=g++
CXXFLAGS=-std=c++11 -fopenmp -O3 

all: main1 main2 main4 main8

main1: conjugate_gradient.cpp sparse_matrix_csr.h util.h
	$(CXX) $(CXXFLAGS) -DTHREADS=1 -o $@ $^

main2: conjugate_gradient.cpp sparse_matrix_csr.h util.h
	$(CXX) $(CXXFLAGS) -DTHREADS=2 -o $@ $^

main4: conjugate_gradient.cpp sparse_matrix_csr.h util.h
	$(CXX) $(CXXFLAGS) -DTHREADS=4 -o $@ $^

main8: conjugate_gradient.cpp sparse_matrix_csr.h util.h
	$(CXX) $(CXXFLAGS) -DTHREADS=8 -o $@ $^

run:
	./main8 matrix x x0
	./main4 matrix x x0
	./main2 matrix x x0
	./main1 matrix x x0

clean:
	rm -f ./main*
