# makefile for optimization framework
# written by Ruimin Wang
# powered by 'MathU'
# copyright@MathU

include ./Makefile.in

all: matrix simplex vecmat spmat umfpack umfpackwrapper omp ls nls eigen lapack proximal mpitest lasso eigensolver iotest lpproblem lu 

vector: bin/vector
bin/vector: obj/vector.o
	$(CXX) obj/vector.o -L$(DIR_LIB) -o $@ $(LIB)
obj/vector.o: examples/vector.cpp 
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

	
matrix: bin/matrix
bin/matrix: obj/matrix.o
	$(CXX) obj/matrix.o -L$(DIR_LIB) -o $@ $(LIB)
obj/matrix.o: test/matrix.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

simplex: bin/simplex
bin/simplex: obj/simplex.o
	$(CXX) obj/simplex.o -L$(DIR_LIB) -o $@ $(LIB)
obj/simplex.o: test/simplex.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

vecmat: bin/vecmat
bin/vecmat: obj/vecmat.o
	$(CXX) obj/vecmat.o -L$(DIR_LIB) -o $@ $(LIB)
obj/vecmat.o: test/matrix_vector.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

spmat: bin/spmat
bin/spmat: obj/spmat.o
	$(CXX) obj/spmat.o -L$(DIR_LIB) -o $@ $(LIB)
obj/spmat.o:test/spmat.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

umfpack: bin/umfpack
bin/umfpack: obj/umfpack.o
	$(CXX) obj/umfpack.o -L$(DIR_LIB) -o $@ $(LIB)
obj/umfpack.o:test/umfpack.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

umfpackwrapper: bin/umfpackwrapper
bin/umfpackwrapper: obj/umfpackwrapper.o
	$(CXX) obj/umfpackwrapper.o -L$(DIR_LIB) -o $@ $(LIB)
obj/umfpackwrapper.o: test/umfpack_wrapper.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

omp: bin/omp
bin/omp: obj/omp.o
	$(CXX) obj/omp.o -L$(DIR_LIB) -o $@ $(LIB)
obj/omp.o: test/omp.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

ls: bin/ls
bin/ls: obj/ls.o
	$(CXX) obj/ls.o -L$(DIR_LIB) -o $@ $(LIB)
obj/ls.o: test/ls.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

nls: bin/nls
bin/nls: obj/nls.o
	$(CXX) obj/nls.o -L$(DIR_LIB) -o $@ $(LIB)
obj/nls.o: test/nls.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

eigen: bin/eigen
bin/eigen: obj/eigen.o
	$(CXX) obj/eigen.o -L$(DIR_LIB) -o $@ $(LIB)
obj/eigen.o: test/eigen_umfpack.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

lapack: bin/lapack
bin/lapack: obj/lapack.o
	$(CXX) obj/lapack.o -L$(DIR_LIB) -o $@ $(LIB)
obj/lapack.o: test/lapack_test.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

proximal: bin/proximal
bin/proximal: obj/proximal.o
	$(CXX) obj/proximal.o -L$(DIR_LIB) -o $@ $(LIB)
obj/proximal.o: test/proximal.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

mpitest: bin/mpitest
bin/mpitest: obj/mpitest.o
	mpic++ obj/mpitest.o -L$(DIR_LIB) -o $@ $(LIB)
obj/mpitest.o: test/openmpi_test.cpp
	mpic++ $(CXXFLAGS) $(DIR_INC) -c $< -o $@

lasso: bin/lasso
bin/lasso: obj/lasso.o
	$(CXX) obj/lasso.o -L$(DIR_LIB) -o $@ $(LIB)
obj/lasso.o: test/lasso.cpp include/Algorithms/Lasso.hpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

eigensolver: bin/eigensolver
bin/eigensolver: obj/eigensolver.o
	$(CXX) obj/eigensolver.o -L$(DIR_LIB) -o $@ $(LIB)
obj/eigensolver.o: test/eigensolver.cpp include/LinearAlgo/Eigenvalue.hpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

iotest: bin/iotest
bin/iotest: obj/iotest.o
	$(CXX) obj/iotest.o -L$(DIR_LIB) -o $@ $(LIB)
obj/iotest.o: test/io_test.cpp include/IOs/MtxFile.hpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

lpproblem: bin/lpproblem
bin/lpproblem: obj/lpproblem.o
	$(CXX) obj/lpproblem.o -L$(DIR_LIB) -o $@ $(LIB)
obj/lpproblem.o: test/lp_problem.cpp include/Problems/LPProblem.hpp include/Problems/BPProblem.hpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

lu: bin/lu
bin/lu: obj/lu.o
	$(CXX) obj/lu.o -L$(DIR_LIB) -o $@ $(LIB)
obj/lu.o: test/lu.cpp include/LinearAlgo/LUSolver.hpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

size: bin/size
bin/size: obj/size.o
	$(CXX) obj/size.o -L$(DIR_LIB) -o $@ $(LIB)
obj/size.o: test/size.cpp include/Base/Size.hpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

on: bin/on
bin/on: obj/on.o
	$(CXX) obj/on.o -L$(DIR_LIB) -o $@ $(LIB)
obj/on.o: test/operation_norm.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

svd: bin/svd
bin/svd: obj/svd.o
	$(CXX) obj/svd.o -L$(DIR_LIB) -o $@ $(LIB)
obj/svd.o: test/svd.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

lowrank: bin/lowrank
bin/lowrank: obj/lowrank.o
	$(CXX) obj/lowrank.o -L$(DIR_LIB) -o $@ $(LIB)
obj/lowrank.o: test/lowrank.cpp include/Algorithms/LowRank.hpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

pca: bin/pca
bin/pca: obj/pca.o
	$(CXX) obj/pca.o -L$(DIR_LIB) -o $@ $(LIB)
obj/pca.o: test/pca.cpp include/Algorithms/PCA.hpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@

help: $(TEST_BIN) 
	
$(TEST_BIN): $(TEST_OBJ)
	$(CXX) $(TEST_OBJ) -L$(DIR_LIB) -o $@ -lcblas -lblas
$(TEST_OBJ):$(TEST_SRC)
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c $< -o $@


.PHONY:clean
clean:
	rm -f $(DIR_OBJ)/*.o
