#ifndef TYPES_DEFINITIONS
#define TYPES_DEFINITIONS

#include <complex>
#include <limits>

namespace kpm
{

typedef double real;
//typedef float real;
typedef std::complex<real> complex;
typedef int integer;
}
const kpm::complex I(0.,1.);



namespace tblinalg
{
#define	MPI_TBLINALG_COMPLEX MPI_DOUBLE_COMPLEX

typedef double real;
//typedef float real;
typedef std::complex<real> complex;
typedef int integer;
real RealZero=std::numeric_limits<real>::epsilon();

}
//const kpm::complex I(0.,1.);

#endif
