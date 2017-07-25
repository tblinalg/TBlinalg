

#ifndef MPI_PROCESS_GRID_HPP
#define MPI_PROCESS_GRIDL_HPP



//C libraries
#include <cstdlib>
#include <cassert>
//C++99 libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
//Custom Libraries
#include "types_definitions.hpp"
#include "input_output.hpp"



namespace MyMPI
{


struct ProcessGrid
{
	public:
	
	size_t BC[3];

	ProcessGrid():	
	spatialDim_(3)
	{
		BC[0]=1; BC[1]=1; BC[2]=1;
	}
	
	ProcessGrid(const size_t _spatialDim):	
	spatialDim_(_spatialDim)
	{
		BC[0]=1; BC[1]=1; BC[2]=1;
	}
	


	size_t GriDim() const
	{
		return gridDim_;
	}

	size_t ProcessID() const
	{
		return processID_;
	}

	size_t BlockNumber(const size_t i) const
	{
		assert( i< spatialDim_ );
		return blockNum_[i];
	}

	void SetGridDim(const size_t _gridDim)
	{
		gridDim_=_gridDim;
	}

	void SetBlockNumber(const size_t dim0,const size_t dim1, const size_t dim2)
	{
		blockNum_[0]=dim0;
		blockNum_[1]=dim1;
		blockNum_[2]=dim2;
		
	}

	void SetProcessID(size_t _processID)
	{
		processID_=_processID;
	}

	
	void Print()
	{
		std::cout<<std::endl<<"The grid of processes is distributed as: "<<blockNum_[0]<<" "<<blockNum_[1]<<" "<<blockNum_[2]<<std::endl;
		
	};
	


	void savetxt(std::string filename)
	{
		std::ofstream os(filename.c_str(),std::ofstream::binary );	
		os<<blockNum_[0]<<" "<<blockNum_[1]<<" "<<blockNum_[2];
		os.close();
	}

	void loadtxt(std::string filename)
	{
		if ( !fileExists(filename) )
		{
			std::cout<<"the .proc file: "<< filename<<" not found. ABORTING "<<std::endl;
			std::exit(-1);
		}
		std::ifstream is(filename.c_str(),std::ifstream::binary );	
		for(size_t i=0; i< spatialDim_; i++)

		is>>blockNum_[0]>>blockNum_[1]>>blockNum_[2];
		size_t gridim=blockNum_[0]*blockNum_[1]*blockNum_[2];
		assert( gridim>0); 
		SetGridDim( gridim );

		is.close();
	}

	void GetFromCFG(std::string filename)
	{
		if ( !fileExists(filename) )
		{
			std::cout<<"the .cfg file: "<< filename<<" not found. ABORTING "<<std::endl;
			std::exit(-1);
		}
		std::ifstream is(filename.c_str(),std::ifstream::binary );	
		
		
		const std::string header("PROCESS_GRID");
		std::string line;
		bool header_found=false;
		while( getline(is, line ) ) 
		{
			if (line.find(header, 0) != std::string::npos)
			{
				int bdim0,bdim1,bdim2;
				is>>bdim0>>bdim1>>bdim2;
				SetBlockNumber(bdim0,bdim1,bdim2);
				header_found=true;
			}
		}	
		is.close();

		if( !header_found )
		{
			std::cout << "The header " << header<<" was not found, therefore aborting"<< std::endl;
			std::exit(-1);
		}
	}

	private:
	size_t processID_;
	size_t blockNum_[3];
	size_t gridDim_;
	size_t spatialDim_;
};
};

#endif
