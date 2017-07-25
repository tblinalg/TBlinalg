

//C libraries
#include <cstdlib>
#include <cassert>
//C++99 libraries
#include <iostream>
#include <vector>
//Custom Libraries
#include "types_definitions.hpp"
#include "tblinalg/indexes.hpp"
#include "tblinalg/structural.hpp"
#include "tblinalg/tb_operator_writer.hpp"

//Exclusive libraries
#include "lattice_geometry.h"


int main(int argc, char *argv[])
{
	//Input parameters
	const size_t spatialDim = 2;
	const size_t colBufferSize= 3;
	const size_t NumOfIndexes=5;
	tblinalg::Indexes idx(NumOfIndexes);
	tblinalg::real Lat[spatialDim][spatialDim];
	for (int i=0;i<spatialDim ;i++)
	for (int j=0;j<spatialDim ;j++)
		Lat[i][j]= A[i][j];
	const size_t norb=2;
	const size_t nspin=1;
	

	if( argc != 3+1)
	{
		std::cout<<"Please submit: label, dim0, dim1 "<<std::endl;
		return 0;
	}
	std::string label(argv[1]);
	idx.Dim(0) = atoi(argv[2]);
	idx.Dim(1) = atoi(argv[3]);
	idx.Dim(2) = 1;
	idx.Dim(3) = norb;
	idx.Dim(4) = nspin;

	// Compute the operator dimension
	size_t OpDim=1; 
	for(size_t i=0;i< NumOfIndexes;i++ )
		OpDim*= idx.Dim(i) ;


	//Save the structural file name
	std::string outputFilename = label+".Structural"; 
	tblinalg::Structural strucParams;


	for(size_t i=0; i< spatialDim; i++)
	{
		strucParams.dim[i]=idx.Dim(i);
		for(size_t j=0; j< spatialDim; j++)
			strucParams.lat[i][j]=Lat[i][j];
	}
	strucParams.norb=norb;
	strucParams.nspin=nspin;
	strucParams.printParams();
	strucParams.savetxt(outputFilename);


	const Complex t1      = -1;
	tblinalg::TBOperatorWriter Ham(OpDim ,colBufferSize);


	tblinalg::Indexes colIdx= idx;
	for(idx(0)=0;idx(0)<idx.Dim(0); idx(0)++)	
	for(idx(1)=0;idx(1)<idx.Dim(1); idx(1)++)	
	for(idx(2)=0;idx(2)<idx.Dim(2); idx(2)++)	
	for(idx(3)=0;idx(3)<idx.Dim(3); idx(3)++)	
	for(idx(4)=0;idx(4)<idx.Dim(4); idx(4)++)	
	{	
		const tblinalg::integer 
		io = idx(3) ,
		jo=1-io;

		const tblinalg::integer 
		is = idx(4) ,
		js=1-is;

		//Compute the sign associate with each index
		const tblinalg::real	
		ss= 1 - 2*is,	//is=(0)1 implies Sz=(1)-1
		oo= 1 - 2*io ;

		//define the indexes for the nearest neighbors
		const tblinalg::integer
		DI0nn[3]={0,oo,0 },
		DI1nn[3]={0,0 ,oo};
		//define the indexes of the next nearest neighbors
		const tblinalg::integer
		DI0nnn[6]={-1,+1, 0, 0 ,-1 ,+1 },
		DI1nnn[6]={ 0, 0,-1,+1 ,+1 ,-1 };

		const tblinalg::integer
		row= idx.ReduceToOne();

		const tblinalg::real 
		ri[3]={ idx(0)*A[0][0] + idx(1)*A[1][0]+ io*Delta[0],
				idx(0)*A[0][1] + idx(1)*A[1][1]+ io*Delta[1], 
				 0 
				};			    

	//----------------------------NEAREST NEIGHBORS/----------------------------/
		for(int n=0;n<3;n++)
		{
			
			tblinalg::real 
			rj[2]=	{ 
					(idx(0)+DI0nn[n])*A[0][0] + (idx(1)+DI1nn[n])*A[1][0]+jo*Delta[0],
					(idx(0)+DI0nn[n])*A[0][1] + (idx(1)+DI1nn[n])*A[1][1]+jo*Delta[1] 
					};

				colIdx(0) = (idx(0) + DI0nn[n] + idx.Dim(0) )%idx.Dim(0);
				colIdx(1) = (idx(1) + DI1nn[n] + idx.Dim(1) )%idx.Dim(1);
				colIdx(2) = 0 ;
				colIdx(3) = jo;
				colIdx(4) = is;

				const size_t
				col= colIdx.ReduceToOne();

				tblinalg::complex val = DI1nn[n];
				Ham.AddEntry(row,col,val);
		}
	}		
	

	Ham.WriteIntoFile( "operators/"+label+".Ham" );
	Ham.saveTxtCoo( "operators/"+label+".Ham");

return 0;}



	
