
#include <cstdlib>
#include <sys/time.h>
#include <vector>
#include <algorithm>    // std::sort
#include <fstream>
#include <string>
#include <complex>
#include <iostream>
#include <iomanip>
#include "string_util.h"

#include "puddle_disorder.hpp"
#include "lattice_geometry.h"
#include "mat_operators.hpp"
#include "hamiltonian_parameters.h"
#include "tb_operator.hpp"

#include <climits>

#include "structural.hpp"

#include <map>	

static const Real tol_zero= std::numeric_limits<Real>::epsilon() ;
static const Real epsilon= std::numeric_limits<Real>::epsilon() ;

int IndexesToIndex( const int i0,const int  n0,
					const int i1,const int  n1,
					const int i2,const int  n2,
					const int i3,const int  n3 
					)
{
	return ( ( (i0+n0)%n0 *n1 + (i1+n1)%n1 )*n2 + (i2+n2)%n2 )*n3 + (i3+n3)%n3;  
//	return (i0+n0)%n0 + (i1+n1)%n1*n0   + i2*n0*n1 + i3*n0*n1*n2;  //old ormat

};
// the index convection is ( js*norb+ jo )		

int main(int argc, char *argv[])
{
	const int spDim=3;
	

	if( argc != 3+1)
	{
		std::cout<<"Please submit: label, dim0, dim1 "<<std::endl;
		return 0;
	}
	std::string label(argv[1]);

	const int dim[3]={atoi(argv[2]), atoi(argv[3]),1};
	const int norb =2;
	const int nspin=1;	
	const int OpDim= dim[0]*dim[1]*norb*nspin;



	const Complex t1      = -1;
	
	
	TBOperator Ham(OpDim , 20 );
	TBOperator DDHam(OpDim , 20 );

	
	//<<<<<<<<<<<<<<<<<<<<DOMAIN DECOMPOSITION>>>>>>>>>>>>>>>>>>>>>>//
	//Number of domains one want to break the cell
	std::map<int,int> idxMap; 	//Map from the standard indexation
										//to the domain indexation
	int DomNum[3] = {32,32,1 };
	int TotNumDom =	DomNum[0]*DomNum[1]*DomNum[2] ;
	int DomSize[3]; 
	for(int i=0;i<3;i++) 
		DomSize[i]=(dim[i]+DomNum[i]-1)/DomNum[i];
	std::cout<<"The domain size is: "<<DomSize[0]<<" "<<DomSize[1]<<" "<<DomSize[2]<<std::endl;
	//	The first step of the domain decomposition is to look within all
	// the sites beloning to the lattice, which one belong to a given
	// domain and index it accordingly
	int domRow=0;
	for(int domIdx0=0; domIdx0< DomNum[0] ; domIdx0++)
	for(int domIdx1=0; domIdx1< DomNum[1] ; domIdx1++)
	{	
		const int domIdx= domIdx0 * DomNum[1] + domIdx1 ;

		//	It is assume that this loop will go through all the sites in the
		// system, no matter is is disconected from the lattice.
		//	Therefore, the first step is to separate the row indices 
		// in domains
		for(int i0=0;i0<dim[0] ; i0++)	
		for(int i1=0;i1<dim[1] ; i1++)
		for(int io=0;io<norb ; io++)
		for(int is=0;is<nspin; is++)
		{	
			const int
			row= IndexesToIndex(i0,dim[0],i1,dim[1],io,norb,is,nspin);
			if(
				(i0>= domIdx0*DomSize[0] && i0< (domIdx0+1)*DomSize[0] )
				&&
				(i1>= domIdx1*DomSize[1] && i1< (domIdx1+1)*DomSize[1] )
				)
			{
				idxMap[row]= domRow ;
				domRow+=1;
			}
		}
	}
	// 	Once all the sites are properly indexed, we can create the
	// hamiltonian in the usual fashion, but usin map instead of row
	// or col
////////////////////////////////////////////////////////////////////////

	for(int i0=0;i0<dim[0] ; i0++)	
	for(int i1=0;i1<dim[1] ; i1++)
	for(int io=0;io<norb ; io++)
	for(int is=0;is<nspin; is++)
	{	
		//Compute the sign associate with each index
		const Real	
		ss= 1 - 2*is,	//is=(0)1 implies Sz=(1)-1
		oo= 1 - 2*io ;
//		//compute the final orbital/spin
		const int
		jo=1-io,
		js=1-is;
		//define the indexes for the nearest neighbors
		int DI0nn[3]={0,oo,0 };
		int DI1nn[3]={0,0 ,oo};
		//define the indexes of the next nearest neighbors
		int DI0nnn[6]={-1,+1, 0, 0 ,-1 ,+1 };
		int DI1nnn[6]={ 0, 0,-1,+1 ,+1 ,-1 };
		//define the inital index
		const int
		row= IndexesToIndex(i0,dim[0],i1,dim[1],io,norb,is,nspin);

	//----------------------------NEAREST NEIGHBORS/----------------------------/
		for(int n=0;n<3;n++)
		{	
			const Real 
			rj[2]={ (i0+DI0nn[n])*A[0][0] + (i1+DI1nn[n])*A[1][0]+jo*Delta[0],
					(i0+DI0nn[n])*A[0][1] + (i1+DI1nn[n])*A[1][1]+jo*Delta[1] 
				  };
					const int
					col= IndexesToIndex(i0+DI0nn[n],dim[0],i1+DI1nn[n],dim[1],jo,norb,is,nspin);
					Complex val = t1;
					Ham.AddEntry(row,col,val);
					DDHam.AddEntry(idxMap[row],idxMap[col],val);
		}	
	}		

	Ham.saveTxtCoo( label+".Ham" );
	DDHam.saveTxtCoo( label+".DD.Ham" );


return 0;}



	
