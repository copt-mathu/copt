default: all

include ../Makefile.in

I = -I../include

all: array vector

array: bin/array
	
bin/array: obj/array.o
	$(CXX) obj/array.o -L$(DIR_LIB) -o $@ $(LIB)
obj/array.o: array.cpp ../include/Base/Array.hpp ../include/Base/ArrayImpl.hpp
	$(CXX) $(CXXFLAGS) $(I) -c $< -o $@

vector: bin/vector
	
bin/vector: obj/vector.o
	$(CXX) obj/vector.o -L$(DIR_LIB) -o $@ $(LIB)
obj/vector.o: vector.cpp ../include/Base/Vector.hpp ../include/Base/VectorImpl.hpp
	$(CXX) $(CXXFLAGS) $(I) -c $< -o $@

matrix: bin/matrix
	
bin/matrix: obj/matrix.o
	$(CXX) obj/matrix.o -L$(DIR_LIB) -o $@ $(LIB)
obj/matrix.o: matrix.cpp ../include/Base/Matrix.hpp ../include/Base/MatrixImpl.hpp ../include/Base/ArrayImpl.hpp ../include/Base/Array.hpp
	$(CXX) $(CXXFLAGS) $(I) -c $< -o $@

solver: bin/solver
	
bin/solver: obj/solver.o
	$(CXX) obj/solver.o -L$(DIR_LIB) -o $@ $(LIB)
obj/solver.o: solver.cpp ../include/Base/Solver.hpp
	$(CXX) $(CXXFLAGS) $(I) -c $< -o $@
