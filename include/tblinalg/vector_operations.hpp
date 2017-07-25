#ifndef TBLINALG_VECTOR_OPERATIONS_HPP
#define TBLINALG_VECTOR_OPERATIONS_HPP

#include <mpi.h>
#include <vector>
#include <algorithm> //for std::sort
#include <cassert>

namespace tblinalg
{
	void Dot(	const size_t root,
				const std::vector<complex>& x,
				const std::vector<complex>& y,
				complex& result)
	{
		const size_t dim = x.size() ;
		complex partial_result= 0 ;
		
		for( std::size_t i = 0 ; i< dim ; i++ )
			partial_result+= std::conj(x[i])*y[i] ;

		result=0.0;		
		MPI_Reduce( &partial_result, &result, 1, 
					MPI_TBLINALG_COMPLEX,
					MPI_SUM,
					root,
					MPI_COMM_WORLD);

	}


	void Dot(	const std::vector<complex>& x,
				const std::vector<complex>& y,
				complex& result)
	{
		const size_t dim = x.size() ;
		complex partial_result= 0 ;
		
		for( std::size_t i = 0 ; i< dim ; i++ )
			partial_result+= std::conj(x[i])*y[i] ;

		result=0.0;		
		MPI_Allreduce( &partial_result, &result, 1, 
						MPI_TBLINALG_COMPLEX,
						MPI_SUM,
						MPI_COMM_WORLD);

	};
	


	void Norm(	const size_t root,
				const std::vector<complex>& x,
				real& result)
	{
		const size_t dim = x.size() ;
		real partial_result= 0 ;
		
		for( std::size_t i = 0 ; i< dim ; i++ )
			partial_result+= std::norm(x[i]);

		MPI_Reduce( &partial_result, &result, 1, 
						MPI_TBLINALG_COMPLEX,
						MPI_SUM,
						root,
						MPI_COMM_WORLD);
		result = sqrt( result ) ;
	};


	void Norm(	const std::vector<complex>& x,
				real& result)
	{
		const size_t dim = x.size() ;
		real partial_result= 0 ;
		
		for( std::size_t i = 0 ; i< dim ; i++ )
			partial_result+= std::norm(x[i]);

		MPI_Allreduce( &partial_result, &result, 1, 
						MPI_TBLINALG_COMPLEX,
						MPI_SUM,
						MPI_COMM_WORLD);
		result = sqrt( result ) ;
	};

	void Normalize( std::vector<complex>& x )
	{
		const size_t dim = x.size() ;
		real  nrm = 0 ;  
		Norm( x, nrm ); 
		
		const real invNorm = 1/nrm ;
		
		for( std::size_t i = 0 ; i< dim ; i++ )
			 x[i]*=invNorm;

	};
	
	
	void FillRandomFill(std::vector<complex>& x)
	{
		const size_t dim = x.size() ;
		for( std::size_t i = 0 ; i< dim ; i++ )
		{
			const real phi = 2.0*M_PI*( (real)rand()/(real)RAND_MAX - 0.5 ) ;
			x[i] = std::exp( I * phi);
		}
	};


	void Fill(std::vector<complex>& x, complex C)
	{
		const size_t dim = x.size() ;
		for( std::size_t i = 0 ; i< dim ; i++ )
			x[i] = C;
	};



	void SyncRandomFill(	const size_t globDim,
					const size_t rowOrig,
					std::vector<complex>& x)
	{
		const size_t dim = x.size() ;
		for( std::size_t i = 0 ; i< globDim ; i++ )
		{
			const integer idx = i - rowOrig ;
			const real phi = 2.0*M_PI*( (real)rand()/(real)RAND_MAX - 0.5 ) ;

			if ( idx >= 0  && idx < dim )
				x[idx] = std::exp( I * phi);
		}
	};	
}
#endif
