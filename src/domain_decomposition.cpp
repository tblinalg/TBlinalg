
//C libraries
#include <cstdlib>
#include <cassert>
//C++99 libraries
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
//Custom Libraries
#include "mpi/process_grid.hpp"
#include "tblinalg/structural.hpp"
#include "tblinalg/indexes.hpp"
#include "tblinalg/domain_decomposition.hpp"
#include "tblinalg/tb_operator.hpp"
//#include "tblinalg/string_util.h"


 
 


int main(int argc, char *argv[])
{
	if( argc != 6+1)
	{
		std::cout<<"Please submit a .cfg file, label , Op, NumDom0, NumDom1, NumDom2  "<<std::endl;
		return 0;
	}
	std::string cfgFilename=argv[1];

	//Structural parameters
	
	tblinalg::Structural strParams;
	strParams.GetFromCFG(cfgFilename);
	const size_t 
	OpDim=	strParams.dim[0]*strParams.dim[1]*strParams.dim[2]*
			strParams.nspin*strParams.norb;
	tblinalg::Indexes idx(5);
	idx.Dim(0) = strParams.dim[0] ;
	idx.Dim(1) = strParams.dim[1] ;
	idx.Dim(2) = strParams.dim[2] ;
	idx.Dim(3) = strParams.norb   ;
	idx.Dim(4) = strParams.nspin  ;	
	strParams.printParams();


	//MPI_STRUCTURe
	MyMPI::ProcessGrid proc;
	proc.GetFromCFG(cfgFilename);
	proc.Print();
	size_t procGridDim =	proc.GriDim();
	tblinalg::Indexes pidx(3);
	pidx.Dim(0)= proc.BlockNumber(0) ;
	pidx.Dim(1)= proc.BlockNumber(1);
	pidx.Dim(2)= proc.BlockNumber(2);


	//DOMAIN DECOMPOSITION 
	tblinalg::Indexes bidx(3);
	std::string label = argv[2];
	std::string OpName= argv[3];
	bidx.Dim(0)= atoi ( argv[4]) ;
	bidx.Dim(1)= atoi ( argv[5]);
	bidx.Dim(2)= atoi ( argv[6]);
	size_t GridDim =	bidx.Dim(0)*bidx.Dim(1)*bidx.Dim(2) ;
	std::map<size_t,size_t> idxMap; 	//Map from the standard indexation
	

	//<<<<<<<<<<<<<<<<<<<<DOMAIN DECOMPOSITION>>>>>>>>>>>>>>>>>>>>>>//
	tblinalg::MapFromLatticeToDD(pidx,bidx,idx, idxMap );

 	tblinalg::Indexes tidx(3);
	for(size_t i=0;i<3;i++) 
		tidx.Dim(i)= (idx.Dim(i) + pidx.Dim(i)-1)/pidx.Dim(i) ; 

	std::ofstream ddoutput((label+".DD").c_str(),std::ofstream::binary );
	ddoutput<<pidx.Dim(0)<<" "<<pidx.Dim(1)<<" "<<pidx.Dim(2)<<" ";
	GridDim =	pidx.Dim(0)*pidx.Dim(1)*pidx.Dim(2) ;
	
	// 	Once all the sites are properly indexed, we can create the
	// hamiltonian in the usual fashion, but usin map instead of row
	// or col
	////////////////////////////////////////////////////////////////////////
	std::string opFilename( ("operators/"+label+"."+OpName ).c_str() );
	std::cout<<"Reading Operator"<<opFilename<<std::endl;

	tblinalg::TBOperator Op; 
	Op.ReadOpFromFile(opFilename);
	std::cout<<"Succeed"<<std::endl;



	//Get the number of elements in each block
	size_t numOrbs=0;
	ddoutput<<numOrbs<<" ";
	for( size_t rank= 0 ; rank< GridDim ; rank ++ ) 
	{
		for( size_t row = 0 ; row < Op.TotNumOrbs() ; row++)
		{
			const tblinalg::MatRow matrow=Op.matGrid.Block().GetMatRow(row);
			const size_t elemsPerRow = matrow.Dim() ;
			bidx.ScatterToAll(rank);
			idx.ScatterToAll(row);

			bool isOrbInBlock=true; 
			for(size_t d=0;d<3;d++)
			{
				tidx(d) = idx(d)/ tidx.Dim(d);
				isOrbInBlock*= (tidx(d) == bidx(d) ); 
			}
			if ( isOrbInBlock )
				numOrbs+=1;
		}
		ddoutput<<numOrbs<<" ";
	}
	ddoutput.close();


	
	tblinalg::Indexes rowIdx=idx;
	tblinalg::Indexes colIdx=idx;
	tblinalg::Indexes trowIdx=tidx;
	tblinalg::Indexes tcolIdx=tidx;

	std::cout<<"Thread Dim "<<bidx.Dim(0)<<" "<<bidx.Dim(1)<<std::endl;
	int dr[3]={ 0, 0, 0 } ;
	for( dr[0]=-1;dr[0]<=1; dr[0]++)	
	for( dr[1]=-1;dr[1]<=1; dr[1]++)
	//if ( dr[1] == -1 && dr[0] == 0)
	{
		std::string OpOutput( "operators/"+label+"."+OpName+"["+tblinalg::ToString(dr[0])+","+tblinalg::ToString(dr[1])+"]" );
		tblinalg::TBOperator  GridOp = Op;
		
		for( size_t row = Op.RowOrigin() ; row < Op.RowDim() ; row++)
		{
			rowIdx.ScatterToAll(row);
			
			//Create a row to store the elements 
			const tblinalg::MatRow matrow=Op.matGrid.Block().GetMatRow(row);
			const size_t elemsPerRow = matrow.Dim() ;

			tblinalg::MatRow outputMatRow;	
			for( size_t elidx=0; elidx< elemsPerRow; elidx++)
			{	

				const tblinalg::MatTuple t = matrow.GetMatTuple(elidx);
				colIdx.ScatterToAll(t.col);
				
				bool isMyTuple=true;
				for(size_t d=0;d<3;d++)
				{
					trowIdx(d) = rowIdx(d)/ trowIdx.Dim(d); 
					tcolIdx(d) = colIdx(d)/ tcolIdx.Dim(d); 
					int CheckEdge =( trowIdx(d) + dr[d] + bidx.Dim(d) )%bidx.Dim(d); 
					trowIdx(d) = CheckEdge ;
					isMyTuple *= (tcolIdx(d) ==trowIdx(d) );
				}
				
					if (isMyTuple)
				{
					if ( row ==1	 ) 
					std::cout	<<"DX="<<dr[0]<<" DY="<<dr[1]<<std::endl
								<<"Checking row->col: ("<<row<<"->"<<t.col<<")"<<std::endl
								<<"Mapped Version row->col: ("<<idxMap[row]<<"->"<<idxMap[t.col]<<")"<<std::endl
								<<"with rowBlock indices "<<trowIdx(0)<<" "<<trowIdx(1)<<" | "
								<<"with colBlock indices "<<tcolIdx(0)<<" "<<tcolIdx(1)							
								<<std::endl<<std::endl;
						outputMatRow.AddMatTuple( idxMap[t.col], t.val );
				}	
			}
			GridOp.matGrid.Block().SetMatRow(idxMap[row],outputMatRow);
		}

		GridOp.WriteIntoFile(OpOutput);
		GridOp.saveTxtCoo( OpOutput);
	}
		



return 0;}


