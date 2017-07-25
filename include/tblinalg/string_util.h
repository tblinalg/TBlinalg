#ifndef TBLINALG_STRING_UTIL
#define TBLINALG_STRING_UTIL

#include <string>
#include <sstream>
#include "types_definitions.hpp"

namespace tblinalg
{
	std::string ToString(const int val)
	{
		std::stringstream ss;
		ss<<val;
		return ss.str();
	};

};





#endif
