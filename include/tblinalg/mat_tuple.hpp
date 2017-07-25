#ifndef TBLINALG_MAT_TUPLE_HPP
#define TBLINALG_MAT_TUPLE_HPP

#include "types_definitions.hpp"

namespace tblinalg{

struct MatTuple
{
	public:
	size_t col;
	complex val;

	//<<<<<<<<<<<<<<<CONSTRUCTOR>>>>>>>>>>>>>>>>>>>//
	MatTuple():
	 col(0), val(0) {}

	MatTuple(const  size_t _col,const  complex _val):
	 col(_col), val(_val) {}

};


struct sortTuplebyCol 
{ 
    bool operator()(MatTuple const &a, MatTuple const &b)
    const{ 
        return a.col < b.col ;
		}
};

};

#endif
