#mpicxx -O3 -openmp -lmpi  -Iinclude parallelMVmul.cpp  -o parallelMVmul
mpicxx -O3 -openmp -lmpi  -Iinclude/ -Iinclude/tblinalg -Iinclude/quantum_correlators dos_test.cpp  -o dos_test
