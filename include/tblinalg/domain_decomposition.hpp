

#ifndef TBLINALG_DOMAIN_DECOMPOSITION_HPP
#define TBLINALG_DOMAIN_DECOMPOSITION_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>	
#include "tblinalg/input_output.hpp"
#include "tblinalg/indexes.hpp"
#include <map>

namespace tblinalg
{

	void MapFromLatticeToDD(Indexes idx,Indexes bIdx,Indexes tIdx, std::map<size_t,size_t>& mymap )
	{
		//	The first step of the domain decomposition is to look within all
		// the sites beloning to the lattice, which one belong to a given
		// domain and index it accordingly
		size_t rowIdx=0;
		for(bIdx(0)=0; bIdx(0)< bIdx.Dim(0) ; bIdx(0)++)
		for(bIdx(1)=0; bIdx(1)< bIdx.Dim(1) ; bIdx(1)++)
		for(bIdx(2)=0; bIdx(2)< bIdx.Dim(2) ; bIdx(2)++)
		{	

			//	It is assume that this loop will go through all the sites in the
			// system, no matter is is disconected from the lattice.
			//	Therefore, the first step is to separate the row indices 
			// in domains
			for(idx(0)=0;idx(0)<idx.Dim(0) ; idx(0)++)	
			for(idx(1)=0;idx(1)<idx.Dim(1) ; idx(1)++)
			for(idx(2)=0;idx(2)<idx.Dim(2) ; idx(2)++)
			for(idx(3)=0;idx(3)<idx.Dim(3) ; idx(3)++)
			for(idx(4)=0;idx(4)<idx.Dim(4) ; idx(4)++)
			{	
				const size_t
				row= idx.ReduceToOne();
				if(
					(idx(0)>= bIdx(0)*tIdx.Dim(0) && idx(0)< (bIdx(0)+1)*tIdx.Dim(0) )
					&&
					(idx(1)>= bIdx(1)*tIdx.Dim(1) && idx(1)< (bIdx(1)+1)*tIdx.Dim(1) )
					&&
					(idx(2)>= bIdx(2)*tIdx.Dim(2) && idx(2)< (bIdx(2)+1)*tIdx.Dim(2) )
					)
				{
					mymap[row]= rowIdx ;
					rowIdx+=1;
//					std::cout<<"map"<<row<<" "<<rowIdx<<std::endl;
				}
			}
		}
	};

};
#endif
