#ifndef TBLINALG_MAT_GRID_HPP
#define TBLINALG_MAT_GRID_HPP

//C libraries
#include <cassert>
//C++99 libraries
#include <vector>
//Custom Libraries
#include "mat_block.hpp"
#include "types_definitions.hpp"
#include "structural.hpp"


struct BlockParameters
{
	int blockIdx[3];
	int blockDim[3];
	
};


namespace tblinalg
{

class MatGrid
{
	public:
	MatGrid()
	{
		std::vector< MatBlock > matGrid_(1);
	}

	MatGrid(const size_t _spDim):numNei_(pow(3,_spDim))
	{
		matGrid_= std::vector< MatBlock >( numNei_ );
	}


	MatBlock Block(const size_t mIdx=0 ) const
	{
		assert( mIdx< numNei_ );
		return 	matGrid_[mIdx];		
	}


	MatBlock& Block(const size_t mIdx=0 )
	{
		assert( mIdx< numNei_ );
		return 	matGrid_[mIdx];		
	}


	MatBlock operator()(const integer i0,const integer i1 ) const
	{
		kpm::integer mIdx=MooreIDX(i0, i1 );
		return 	matGrid_[mIdx];		
	}


	MatBlock& operator()(const integer i0,const integer i1 )
	{
		kpm::integer mIdx=MooreIDX(i0, i1 );
		return 	matGrid_[mIdx];		
	}
	
	size_t MooreIDX(const integer i0,const integer i1 ) const
	{
		return   4 + i0*3 + i1 ; 						
	}

		size_t numNei_;
		std::vector< MatBlock > matGrid_;
		std::vector<size_t> blockLength;	
		std::vector<size_t> blockOrg;	

};
};
#endif
