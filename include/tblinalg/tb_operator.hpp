/// @brief This method adds two integers.
/// @param a First integer to add.
/// @param b Second integer to add.
/// @return The sum of both parameters.

#ifndef TBLINALG_TBOPERATOR_HPP
#define TBLINALG_TBOPERATOR_HPP

//C libraries
#include <cstdlib>
#include <cassert>
#include <mpi.h>
//C++99 libraries
#include <iostream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
//Custom Libraries
#include "types_definitions.hpp"
#include "tblinalg/structural.hpp"
#include "tblinalg/string_util.h"
#include "tblinalg/input_output.hpp"
#include "tblinalg/mat_tuple.hpp"
#include "tblinalg/mat_row.hpp"
#include "tblinalg/mat_grid.hpp"
#include "tblinalg/indexes.hpp"


namespace tblinalg{	class TBOperator{
	
	//Public variables
	public:
	MatGrid matGrid;
	
	
//<<<<<<<<<<<<<<<<<<<<<<<< CLASS CONSTRUCTORS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//

	//Constructor for the reader phase
	public:
	TBOperator( ) :
	numNonZero_(0),bIdx(3),matGrid(0), 
	rowDim_(0), rowOrig_(0), rowEnd_(0)
	{}

	//Constructor for the writer phase
	TBOperator(const size_t _spatial_dim ) :
	numNonZero_(0),bIdx(3),matGrid(_spatial_dim), 
	rowDim_(0), rowOrig_(0), rowEnd_(0)
	{}

	//Constructor for the filling phase
	TBOperator(const size_t _OPDim, const size_t _colBufferSize) :
	numNonZero_(0),bIdx(3),matGrid(0) , 
	rowDim_(0), rowOrig_(0), rowEnd_(0)
	{
		colBufferSize_ = _colBufferSize ;
		matGrid.Block().SetDim(_OPDim);
		SetRowOrigin(0);
		SetRowDim(_OPDim);	
	}
//<<<<<<<<<<<<<<<<<<<<<<<< GETTERS CONSTRUCTORS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
	size_t BufferSize() const  {  return colBufferSize_;  }

	size_t Dim() const { return matGrid.Block().Dim(); }
	
	size_t TotNumOrbs() const {  return totNumOrbs_; }

	size_t RowDim() const { return rowDim_; }  

	size_t RowOrigin() const { return rowOrig_; }

	size_t RowEnd() const { return rowEnd_; }

//<<<<<<<<<<<<<<<<<<<<<<< SETTERS CONSTRUCTORS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
	void SetTotNumOrbs(const size_t _totNumOrbs)
	{ 
		totNumOrbs_=_totNumOrbs;
	}

	void SetRowOrigin(const size_t _rowOrig)   
	{
		rowOrig_=_rowOrig;
		if ( rowOrig_ >= TotNumOrbs() )
			rowOrig_= TotNumOrbs() ;

		rowEnd_=rowOrig_ + rowDim_;
		if (rowEnd_ >= TotNumOrbs() )
		{
			rowEnd_= TotNumOrbs() ;
			std::cout<<"A mistake was made before!"<<std::endl;
		}
	}


	void SetRowEnd(const size_t _rowEnd)
	{
		rowEnd_=_rowEnd;
		if (rowEnd_ >= TotNumOrbs() )
		rowEnd_=TotNumOrbs();
	}

	void SetRowDim(const size_t _rowDim)
	{
		rowDim_=_rowDim;
		rowEnd_= RowOrigin()+ _rowDim;
		if (rowEnd_ >= TotNumOrbs() )
			rowEnd_ = TotNumOrbs() ;
	}

//<<<<<<<<<<<<<<<<<<<<<< COMPLEX PUBLIC METHODS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>//

	void SetRowRange( const std::string ddFilename)  ;
						
						
	void AddEntry(const size_t row,const size_t col,const complex val );	

	void Rescale(	const size_t dx, const size_t dy , 
					const tblinalg::real shift,	const tblinalg::real scal ); 
	
	void Multiply(	const integer dx,
					const integer dy,
					const real a, 
					const std::vector<complex>& X, 
					const real b, 
					std::vector<complex>& Y 
					);
	
	void Multiply(	const real a, 
					const std::vector<complex>& X, 
					const real b,
					std::vector<complex>& Y 
					) ;

	size_t NumNonZero( ); 
	
	size_t NumNonZeroInRow(size_t row ) ;

	void ReadOpFromFile(const int dx, const int dy, const std::string opFilename);	

	void WriteIntoFile(std::string label ) ;
	
	void saveTxtCoo(std::string label ); 

	void PrintTBOp(const integer dx,const integer dy);

//<<<<<<<<<<<<<<<<<<<<<< RECYCLE METHODS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>//

	void ReadOpFromFile(  const std::string opFilename)
	{
		ReadOpFromFile(-1,-1, opFilename);
	};

	void PrintTBOp() 
	{
		PrintTBOp(-1,-1);
	}
	
	void Rescale( const tblinalg::real shift,	const tblinalg::real scal )
	{	
		Rescale(-1,-1,  shift, scal ); 
	}; 

//<<<<<<<<<<<<<<<<<<<<<<MPI METHODS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
	void SetMpiWorkspace( Indexes _bIdx)
	{ 
		bIdx			= _bIdx;
		int sp=2;
		int PBC[3]={1,1,1};
		int blockNum[3]={bIdx.Dim(0),bIdx.Dim(1),bIdx.Dim(2)};

		MPI_Cart_create(MPI_COMM_WORLD, sp , blockNum, PBC,0,&cartComm);

		destX			= std::vector<complex>(RowDim());
		sourceX			= std::vector<complex>(RowDim());
		sourceColOrig	= RowOrigin();
		destColOrig		= RowOrigin();
	};


	private: 
	size_t spatialDim_, rowDim_, rowOrig_, rowEnd_, colBufferSize_, numNonZero_;
	size_t	totNumOrbs_;
	
	//MPI VARIABLES
	size_t sourceColOrig;
	size_t destColOrig;
	std::vector<complex> destX;
	std::vector<complex> sourceX;
	MPI_Comm cartComm;
	Indexes bIdx;
	
};};

#include "tboperator/add_entry.cpp"
#include "tboperator/rescale.cpp"
#include "tboperator/multiply.cpp"
#include "tboperator/print_tb_op.cpp"
#include "tboperator/save_txt_coo.cpp"
#include "tboperator/num_non_zero.cpp"
#include "tboperator/num_non_zero_in_row.cpp"
#include "tboperator/read_op_from_file.cpp"
#include "tboperator/write_into_file.cpp"
#include "tboperator/set_row_range.cpp"


#endif


