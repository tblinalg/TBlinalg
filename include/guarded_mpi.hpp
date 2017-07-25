#ifndef GUARDED_MPI_HPP
#define GUARDED_MPI_HPP


#ifndef MPI_VERSION
#define GUARDED_MPI(x) {}
#endif

#ifdef MPI_VERSION
#define GUARDED_MPI(x)  x
#endif


#endif


