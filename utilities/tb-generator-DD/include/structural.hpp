

#ifndef STRUCTURAL_HPP
#define STRUCTURAL_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>	

class Structural
{
	public:
	int dim[3];
	double lat[3][3];
	int norb;
	int nspin; 
		
	void savetxt(std::string filename)
	{
		std::ofstream os(filename.c_str(),std::ofstream::binary );	
		for(int i=0; i< 3; i++)
		{
			os<<dim[i]<<" ";
			for(int j=0; j< 3; j++)
				os<<lat[i][j]<<" ";
		}
		os<<norb<<" ";
		os<<nspin<<" ";
		os.close();
	};

	void loadtxt(std::string filename)
	{
		std::ifstream is(filename.c_str(),std::ifstream::binary );	
		for(int i=0; i< 3; i++)
		{
			is>>dim[i];
			for(int j=0; j< 3; j++)
				is>>lat[i][j];
		}
		is>>norb;
		is>>nspin;
		is.close();
	}
	
	void printParams()
	{
	std::cout<<"The structural parameters are:"<<std::endl
			 <<"lattice vectors :"<<std::endl
			 <<"Lat1: ( "<<lat[0][0]<<" , "<<lat[0][1]<<" , "<<lat[0][2]<<" )"<<std::endl
			 <<"Lat2: ( "<<lat[1][0]<<" , "<<lat[1][1]<<" , "<<lat[1][2]<<" )"<<std::endl
			 <<"Lat3: ( "<<lat[2][0]<<" , "<<lat[2][1]<<" , "<<lat[2][2]<<" )"<<std::endl
			 <<"Supercell with dimensions : "
			 <<dim[0]<<" x "<<dim[1]<<" x "<<dim[2]<<std::endl
			 <<"Number of orbital per unit cell is : "
			 <<norb<<std::endl
			 <<"Number os spin per orbital spin:"
			 <<nspin
			 <<std::endl;
	}
	
};	  


#endif
