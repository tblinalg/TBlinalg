//C libraries
#include <cstdlib>
#include <mpi.h>
//C++99 libraries
#include <iostream>
//Custom Libraries
#include "tblinalg/structural.hpp"
#include "tblinalg/indexes.hpp"
#include "tblinalg/string_util.h"
#include "tblinalg/tb_operator.hpp"
#include "tblinalg/vector_operations.hpp"
#include "quantum_correlators/kpm.hpp"


int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	int root=0,rank, worldSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	
	if( rank==root)
		std::cout<<"Running parallelMVmul using "<<worldSize<<" processes"<<std::endl;

	if ( argc != 1+ 1 )
	{
		if( rank==root)
			std::cerr<<"This function require the argument LABEL, please submit it."<<std::endl;
		std::exit(-1);
	}
	std::string label=argv[1];
	std::string structFilename=label+".Structural";

	if( rank==root)
		std::cout<<"Reading Structural file: "<<structFilename<<std::endl;

	tblinalg::Structural strPar;
	strPar.loadtxt(structFilename);
	if( rank==root)
		strPar.printParams();
	//Store structural parameters in index idx
	tblinalg::Indexes idx(5);
	idx.Dim(0) = strPar.dim[0] ;
	idx.Dim(1) = strPar.dim[1] ;
	idx.Dim(2) = strPar.dim[2] ;
	idx.Dim(3) = strPar.norb   ;
	idx.Dim(4) = strPar.nspin  ;
	size_t totNumOrbs=idx.Dim(0)*idx.Dim(1)*idx.Dim(2)*idx.Dim(3)*idx.Dim(4);
	size_t spatialDim=2;

	//Store block parameters in index bidx
	tblinalg::Indexes bidx(3);
	std::string DDfilename=label+".DD";
	if ( !tblinalg::fileExists(DDfilename) )
	{
		std::cerr<<"File "<<DDfilename<< " not found, return -1 "<<std::endl;
		std::exit(-1);
	}
	std::ifstream domFile( DDfilename.c_str(),std::ifstream::binary );
	domFile>>bidx.Dim(0)>>bidx.Dim(1)>>bidx.Dim(2);
	bidx.ScatterToAll(rank);
	const size_t gridDim = bidx.Dim(0)*bidx.Dim(1)*bidx.Dim(2) ;

	if( rank==root)
	std::cout<<"From the DD the number of domains is: "
			 <<bidx.Dim(0)<<" "<<bidx.Dim(1)
			 <<" "<<bidx.Dim(2)<<std::endl;

	//Store thread parameters in index bidx
	tblinalg::Indexes tidx(3);
	for( size_t d=0; d<3 ;d++)
		tidx.Dim(d) = idx.Dim(d)/bidx.Dim(d);
	const size_t blockDim = tidx.Dim(0)*tidx.Dim(1)*tidx.Dim(2) ;

	if( rank==root)
		std::cout<<"BlockDim: "<<blockDim<<std::endl;


	
	if ( gridDim != worldSize )
	{
		if(rank==0)
			std::cerr<<"This version of the program, only support worldsiz==gridDim"<<std::endl;
		MPI_Finalize( );	
		return 0;
	}


	tblinalg::TBOperator Op(spatialDim);
	Op.SetRowRange(label); 

	//Reads the on-site matrix and its neighborhood, corresponding to this block
	tblinalg::real
	Emin = -4 ,
	Emax =  4 ,
	alpha= 0.9 ,
	Escal = 2.0*alpha/(Emax-Emin),
	Esshit= 0.5*(Emax+Emin);

	for(int dx=-1;dx<=1; dx++)	
	for(int dy=-1;dy<=1; dy++)
	{
		std::string OpName="Ham";
		std::string OpOutput( "operators/"+label+"."+OpName+"["+tblinalg::ToString(dx)+","+tblinalg::ToString(dy)+"]" );
		Op.ReadOpFromFile(dx,dy, OpOutput );
		Op.Rescale(dx,dy,Esshit,Escal);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	Op.SetMpiWorkspace(bidx);

	const size_t 
	effDim 	=Op.RowDim(),
	globalDim 	=totNumOrbs ;
	
	srand(  544545 );
	
	std::vector< tblinalg::complex>
	Psi(effDim,0), jm1(effDim,0) , jm2(effDim,0) ,
	*pjm1=&jm1, *pjm2=&jm2 ;

	if ( rank == 0 )
		Psi[0]=1;

	//tblinalg::SyncRandomFill(globalDim,rowOrig, Psi );
	//tblinalg::Normalize( Psi );

	std::vector<tblinalg::real> output(globalDim);
	jm1=Psi;
	size_t Mom = 1024 ;
	std::vector<tblinalg::complex> muPar( Mom ,0  ); 
	pjm1=&jm1; pjm2=&jm2;
	for(size_t m=0;m< Mom ; m++ )
	{
		tblinalg::complex mutmp=0;
		cheb_evolve( m , Op, pjm1, pjm2 );
		tblinalg::Dot(root, Psi, *pjm1 , mutmp);
		muPar[m]+=mutmp;
	}

	if (root == rank )
	{
		std::ofstream testfile( "dos.dat" ) ;
		for( double x=-0.9;x< 0.9 ; x=x+0.001 )
		{
			double dos = 0.5*muPar[0].real() ;
			for(size_t m=1;m< Mom ; m++ )
				dos+= muPar[m].real()*cos( m* acos( x  ) ) ;
		
			dos=dos/sqrt(1 - x*x );
			testfile<< x/Escal <<" "<< dos <<std::endl;
		}
	}	


	MPI_Finalize( );
	return 0;

}
