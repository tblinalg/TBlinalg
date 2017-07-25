

#ifndef INPUT_OUTPUT_HPP
#define INPUT_OUTPUT_HPP

#include <sys/stat.h>

inline bool fileExists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
};


#endif


