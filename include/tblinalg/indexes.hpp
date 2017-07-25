#ifndef TBLINALG_INDEXES_HPP
#define TBLINALG_INDEXES_HPP


//C libraries
#include <cassert>
//C++99 libraries
#include <vector>
//Custom Libraries
#include "types_definitions.hpp"


namespace tblinalg{

class Indexes
{
	public:

	//Class constructor 
	Indexes(const size_t _numIndexes):
	numIndexes_(_numIndexes)
	{
		idxArray_= std::vector<integer>( _numIndexes );
		dimArray_= std::vector<size_t>( _numIndexes );
	}



	integer operator()( const  size_t pos) const
	{
		assert( pos < NumOfIndexes());
		return idxArray_[pos];
	}  

	integer& operator()( const  size_t pos)
	{
		assert( pos < NumOfIndexes());
		return idxArray_[pos];
	}  


	size_t Dim(const size_t pos) const 
	{
		assert( pos < NumOfIndexes());
		return dimArray_[pos];
	}

	size_t& Dim(const size_t pos) 
	{
		assert( pos < NumOfIndexes());
		return dimArray_[pos];
	}


	size_t NumOfIndexes() const
	{
		return numIndexes_ ;
	}
	
	integer ReduceToOne() const
	{
		integer singleIdx=idxArray_[0];
		
		for ( size_t pos=1 ; pos < NumOfIndexes(); pos++ ) 
			singleIdx = singleIdx*dimArray_[pos] + idxArray_[pos] ;

		return singleIdx;  
	}

	integer ScatterToAll(size_t idx ) 
	{
		for ( size_t pos=0 ; pos < NumOfIndexes() ; pos++ )
		{	const size_t inv_pos=NumOfIndexes()-1-pos;
			idxArray_[inv_pos] =( idx )%dimArray_[inv_pos];
			idx=idx/dimArray_[inv_pos];
		}
	}

	private:
		size_t numIndexes_;
		std::vector<integer> idxArray_;
		std::vector<size_t> dimArray_;
};
};

#endif
