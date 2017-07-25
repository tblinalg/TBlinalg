#ifndef STRING_UTIL
#define STRING_UTIL

#include <string>
#include <sstream>

namespace std
{

	std::string to_string(const int x)
	{
		std::string Result;		// string which will contain the result
		std::ostringstream convert;	// stream used for the conversion
		convert << x;			// insert the textual representation of 'Number' in the characters in the stream
		Result = convert.str();		// set 'Result' to the contents of the stream
		return Result;
	}

}







#endif
