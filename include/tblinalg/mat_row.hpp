#ifndef TBLINALG_MAT_ROW_HPP
#define TBLINALG_MAT_ROW_HPP

//C libraries
#include <cassert>
//C++99 libraries
#include <vector>
#include <algorithm> //for std::sort
//Custom Libraries
#include "mat_tuple.hpp"
#include "types_definitions.hpp"

namespace tblinalg
{
class MatRow
{
	
	public:
	MatRow(){}

	MatRow(const size_t _dim)
	{
		matRow_=std::vector< MatTuple >(_dim);
	}	
	
	void AddMatTuple(  MatTuple  _matTuple )
	{
		matRow_.push_back( _matTuple );
	}

	void AddMatTuple(const size_t _col, complex _val )
	{
		matRow_.push_back( MatTuple(_col , _val ) );
	}

	MatRow SortRow( )
	{
		std::sort( matRow_.begin(), matRow_.end(), sortTuplebyCol());
		return *this;
	}
	

	size_t Dim() const 
	{
		return matRow_.size();
	}

	MatTuple GetMatTuple(const size_t idx ) const 
	{
		assert( idx< Dim() );
		return matRow_[idx];
	}
	

	MatTuple& GetMatTuple(const size_t idx )
	{
		assert( idx< Dim() );
		return matRow_[idx];
	}

	
	std::vector< MatTuple> matRow_;
};

};
#endif
