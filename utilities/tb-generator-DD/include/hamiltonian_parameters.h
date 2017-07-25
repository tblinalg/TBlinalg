
#ifndef HAMILTONIAN_PARAMETERS_H
#define HAMILTONIAN_PARAMETERS_H

#include "lattice_geometry.h"
#include <fstream>
#include <string>
#include <iostream>
//Pristine
static const Complex E0[2]   ={ 0.00 , 0.00 };
static const Complex t1      = 1.0;
static const Complex t2      =0. ;
static const Complex VI[2]   ={ 0.01 ,  0.01 };
static const Complex VR      =  0.0 ;
static const Complex VPIA[2] = { 0., 0.};

struct HamParams
{
	Complex E0[2],t1,t2 ;
	Complex VI[2] , VR , VPIA[2];

	void load_from_file(const std::string filename)
	{
		std::cout<<"Reading hamiltonian from input file: "<<filename<<std::endl;
		std::ifstream input( filename.c_str() );
		std::cout<<"Reading E0=[E0A,E0B]= ";
		input>>E0[0]>>E0[1];
		std::cout<<"[ "<<E0[0]<<","<<E0[1]<<" ]"<<std::endl;
		std::cout<<"Reading t1= ";
		input>>t1;
		std::cout<<t1<<std::endl;
		std::cout<<"Reading t2= ";
		input>>t2;
		std::cout<<t2<<std::endl;
		std::cout<<"Reading VI=[VIA,VIB]= ";
		input>>VI[0]>>VI[1];
		std::cout<<"[ "<<VI[0]<<","<<VI[1]<<" ]"<<std::endl;
		std::cout<<"Reading VR= ";
		input>>VR;
		std::cout<<VR<<std::endl;
		std::cout<<"Reading PIA=[PIAA,PIAB]= ";
		input>>VPIA[0]>>VPIA[1];
		std::cout<<"[ "<<VPIA[0]<<","<<VPIA[1]<<" ]"<<std::endl;
		input.close();
	}

};


#endif
