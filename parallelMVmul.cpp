
/** Simple program to compute the spin_conductivity using KPM + MPI + openmp**/
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>
#include "types_definitions.hpp"
#include "tb_operator.hpp"
#include "string_util.h"
#include "mpi_ofile.hpp"

#include "kpm.hpp"


#include "structural.hpp"

#include <vector_operations.hpp>


int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int rank, worldSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);


	
	if( rank==0)
		std::cout<<"Running parallelMVmul using "<<worldSize<<" processes"<<std::endl;

	
	std::string label=argv[1];
	std::string structFilename=label+".Structural";



	if( rank==0)
		std::cout<<"Reading Structural file: "<<structFilename<<std::endl;
	Structural strPar;
	strPar.loadtxt(structFilename);
	if( rank==0)
		strPar.printParams();


	int n0=strPar.dim[0];
	int n1=strPar.dim[1];
	int norb=strPar.norb;
	int nspin=strPar.nspin;
	kpm::integer totNumOrbs= n0*n1*norb*nspin;
	int spatialDim=2;



	int blockNum[3]={8,8,1};
	int gridDim = blockNum[0]*blockNum[1]*blockNum[2];
	int blockRank = rank ;
	int blockIdx0 = rank/blockNum[1] ;
	int blockIdx1 = rank%blockNum[1] ;
	
	int  threadNum[3] = {n0/blockNum[0] , n1/blockNum[1] , 1/blockNum[2] };
	int  blockDim     = threadNum[2]*threadNum[1]*threadNum[0];
	long threadRank   = ( blockIdx0*threadNum[1] +  blockIdx1 )*blockDim;

	
	kpm::integer
	OrbsPerBlock = norb*nspin;
	const int effDim=blockDim*OrbsPerBlock;

	long rowMin=blockRank*OrbsPerBlock;
	long rowMax=blockRank*OrbsPerBlock + effDim;

	
	if ( gridDim != worldSize)
	{
		if(rank==0)
			std::cerr<<"This version of the program, only support worldsiz==gridDim"<<std::endl;
		MPI_Finalize( );	
		return 0;
	}


	TBOperator Ham(spatialDim);
	Ham.ReadOpFromFile(rowMin,rowMax,"test.Ham");


	kpm::real a=1.0, b=0.0;
	std::vector<kpm::complex> x(effDim,0),y(effDim,0);
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Comm cartComm;
	int PBC[3]={1,1,1};
	MPI_Cart_create(MPI_COMM_WORLD, spatialDim , blockNum, PBC,0,&cartComm);
	Ham.Multiply(cartComm,blockIdx0,blockIdx1,a,x,b,y);




	TBOperator Hamiltonian(2);

	const size_t root=0;;
	const size_t globalDim =totNumOrbs ;
	//Creates a coordinated random vector
	srand(  544545 );


	std::vector<kpm::complex> 
	Psi(effDim,0), jm1(effDim,0) , jm2(effDim,0) ,
	*pjm1=&jm1, *pjm2=&jm2 ;

	tblinalg::SyncRandomFill(globalDim,rowMin, Psi );
	//Normalize the vector
	tblinalg::Normalize( Psi );

	//Pass the vector to the first Chebyshev vector
	jm1=Psi;


	kpm::real 
	Emin = -4 ,
	Emax =  4 ,
	alpha= 0.9 ;
	kpm::integer Mom = 128 ;

	std::vector<tblinalg::complex> mu( Mom ); 

	pjm1=&jm1; pjm2=&jm2;
	for(int m=0;m< Mom ; m++ )
	{
		tblinalg::complex mutmp=0;
		cheb_evolve(cartComm,blockIdx0,blockIdx1, m , Hamiltonian, pjm1, pjm2 );
		tblinalg::Dot(root, Psi, *pjm1 , mutmp);
		mu[m]+=mutmp;
	}


/*
	std::string  chebOutName=label+"ChebMomTrGreen.mom1D"; 
	std::ofstream mom1D(chebOutName.c_str());
	mom1D<<"test"<<" "<<Mom<<" "<<Emin
		 <<" "   <<Emax<<" "<<alpha<<" "<<norb<<" "
		 <<nspin <<" "<<1.0<<std::endl;
*/


/*
	

	x[6]= 1.;
	Ham.Multiply(a,x,b,y,blockIdx0,blockIdx1);
//	for(int row= 0; row < totNumOrbs ; row ++)
//	if ( row >= threadRank && row < threadRank+blockDim)
//		std::cout<<row<<std::endl;
	
	//if( rank==0)
	//{
	//	std::cout<<"The number of domains is 32 32 1 (FIXED BY HAND)"<<std::endl;
	//	std::cout<<" each process will deal with " <<(TotDomNum+worldSize-1)/worldSize<<" domains"<<std::endl;
	//	std::cout<<"In the following implementation the matrix should enter partitioned, therefore abort if nproc!=ndom"<<std::endl;
	//}
	
*/
	MPI_Finalize( );
	return 0;

}
/***************************REEScalando el hamiltoniano******************/
	
