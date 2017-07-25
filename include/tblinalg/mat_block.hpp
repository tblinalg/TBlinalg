#ifndef TBLINALG_MAT_BLOCK_HPP
#define TBLINALG_MAT_BLOCK_HPP

//C libraries
#include <cassert>
//C++99 libraries
#include <vector>
//Custom Libraries
#include "mat_row.hpp"
#include "types_definitions.hpp"

namespace tblinalg
{
class MatBlock
{
	public:
	MatBlock(){	}

	MatBlock(const size_t _dim)
	{
		matBlock_=std::vector< MatRow >(_dim);
	}	


	size_t Dim( ) const { return matBlock_.size(); }


	void SetDim( const size_t _dim )
	{
		matBlock_=std::vector< MatRow >(_dim);
	}
	
	void SetMatBlock( std::vector< MatRow > _matBlock )
	{
		matBlock_=_matBlock;
	}
	
	void SetMatRow(const size_t row, const  MatRow  _matRow )
	{
		assert( row <  Dim() );
		matBlock_[row]=_matRow;
	}

	MatRow& GetMatRow(const size_t row )
	{
		assert( row <  Dim() );
		return matBlock_[row];
	}

	MatRow GetMatRow(const size_t row ) const 
	{
		assert( row <  Dim() );
		return matBlock_[row];
	}

	std::vector< MatRow > matBlock_;


};
};
#endif
