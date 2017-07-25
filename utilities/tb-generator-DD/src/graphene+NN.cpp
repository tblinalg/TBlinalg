
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
	

	if( argc != 7+1)
	{
		std::cout<<"Please submit: label, dim0, dim1 Up xp Rp Seed"<<std::endl;
		return 0;
	}
	std::string _label(argv[1]);
	const int n0 = atoi(argv[2]);
	const int n1 = atoi(argv[3]);
	const Real Up = atof(argv[4]);
	const Real xp = atof(argv[5]);
	const Real Rp = atof(argv[6]);
	const int Seed = atoi(argv[7]);
	const int norb =2;
	const int nspin=1;	
	const int OpDim= n0*n1*norb*nspin;
	const int coord_max= 20;




	std::string label(_label);
	int  dim[spDim]={n0,n1,1};
	Real Lat[spDim][spDim];
	
	for (int i=0;i<spDim ;i++)
	for (int j=0;j<spDim ;j++)
		Lat[i][j]= A[i][j];

	std::cout<<"Creating a tight-binding hamiltonian"
			 <<"for the system: "<<label<<std::endl<<std::endl
			 <<"The lattice vectors are:"<<std::endl
			 <<"Lat1: ( "<<Lat[0][0]<<" , "<<Lat[0][1]<<" , "<<Lat[0][2]<<" )"<<std::endl
			 <<"Lat2: ( "<<Lat[1][0]<<" , "<<Lat[1][1]<<" , "<<Lat[1][2]<<" )"<<std::endl
			 <<"Lat3: ( "<<Lat[2][0]<<" , "<<Lat[2][1]<<" , "<<Lat[2][2]<<" )"<<std::endl
			 <<"Using a supercell with dimensions : "
			 <<dim[0]<<" x "<<dim[1]<<" x "<<dim[2]
			 <<"The number of orbital per unit cell is : "
			 <<norb
			 <<"and the calculation is perfomed considering spin:"
			 <<nspin
			 <<std::endl;
	

	std::string outputFilename = label+".Structural"; 
	
	Structural strucParams;
	for(int i=0; i< spDim; i++)
	{
		strucParams.dim[i]=dim[i];
		for(int j=0; j< spDim; j++)
			strucParams.lat[i][j]=Lat[i][j];
	}
	strucParams.norb=norb;
	strucParams.nspin=nspin;
	
	strucParams.printParams();
	strucParams.savetxt(outputFilename);
	
	Structural strucParamsRead;
	strucParamsRead.loadtxt(outputFilename);
	strucParamsRead.printParams();


	HamParams ham_params;
	ham_params.load_from_file("params.dat");
	//Create the reciprocal lattice vector of the Monkhorst pack
	//<<<<<<<<<<<<<<<<<<<<END INITIAL BLOCK>>>>>>>>>>>>>>>>>>>>>>>>>>>//	
	

	//<<<<<<<<<<<<<<<<<<<<READ HAMILTONIAN PARAMETERS>>>>>>>>>>>>>>>>>>>>>>>>>//
	const Complex E0[2]   ={ ham_params.E0[0] , ham_params.E0[1] };
	const Complex t1      = ham_params.t1;
	const Complex t2      = ham_params.t2 ;
	const Complex VI[2]   ={ ham_params.VI[1], ham_params.VI[0] };
	const Complex VR      =  ham_params.VR ;
	const Complex VPIA[2] = {ham_params.VPIA[0],ham_params.VPIA[1]};
	//<<<<<<<<<<<<<<<<<<<<ENDREAD HAMILTONIAN PARAMETERS>>>>>>>>>>>>>>>>>>>>>>//


	TBOperator Ham(OpDim ,coord_max);
	TBOperator VelX(OpDim,coord_max);
	TBOperator VelY(OpDim,coord_max);
	TBOperator VelXSz(OpDim,coord_max);
	TBOperator VelYSz(OpDim,coord_max);


	for(int i0=0;i0<dim[0] ; i0++)	
	for(int i1=0;i1<dim[1] ; i1++)
	for(int io=0;io<norb ; io++)
	for(int is=0;is<nspin; is++)
	{	
		//Compute the sign associate with each index
		const Real	
		ss= 1 - 2*is,	//is=(0)1 implies Sz=(1)-1
		oo= 1 - 2*io ;
		//compute the final orbital/spin
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
		row= IndexesToIndex(i0,n0,i1,n1,io,norb,is,nspin);
		const Real 
		ri[3]={ i0*A[0][0] + i1*A[1][0]+ io*Delta[0],
				i0*A[0][1] + i1*A[1][1]+ io*Delta[1], 
				 0 
				};			    
	//----------------------------ONSITE ENERGY/----------------------------/
			{
				const int
				col= IndexesToIndex(i0,n0,i1,n1,io,norb,is, nspin);
				Complex val = E0[io];
				Ham.AddEntry(row,col,val);
			}
	//----------------------------NEAREST NEIGHBORS/----------------------------/
			for(int n=0;n<3;n++)
			{	
				const Real 
				rj[2]={ (i0+DI0nn[n])*A[0][0] + (i1+DI1nn[n])*A[1][0]+jo*Delta[0],
						(i0+DI0nn[n])*A[0][1] + (i1+DI1nn[n])*A[1][1]+jo*Delta[1] 
						};
				const int
				col= IndexesToIndex(i0+DI0nn[n],n0,i1+DI1nn[n],n1,jo,norb,is,nspin);
				Complex val = t1;
				Ham.AddEntry(row,col,val);
				VelX.AddEntry(row,col,I*(rj[0]-ri[0])*val);
				VelY.AddEntry(row,col,I*(rj[1]-ri[1])*val);
				VelXSz.AddEntry(row,col,ss*I*(rj[0]-ri[0])*val);
				VelYSz.AddEntry(row,col,ss*I*(rj[1]-ri[1])*val);				
				}	
		}		
	
	//Add diagonal disorder


	if( xp > epsilon && Up > epsilon)
	{
		std::vector<Real> diagElements(OpDim,0.);
		srand(Seed);
		Real puddHeight= Up;
		Real puddConcen= xp;
		Real puddRange = Rp;
	
		PuddleDisorder(	n0, n1, norb, nspin, 
				puddHeight, puddConcen, puddRange,
				diagElements,Seed*45);	
		double onsiteAv=0;
		for(int i0=0;i0<n0   ; i0++)	
		for(int i1=0;i1<n1   ; i1++)
		for(int io=0;io<norb ; io++)
		for(int is=0;is<nspin; is++)
		{	
			const int
			k0= IndexesToIndex(i0,n0,i1,n1,io,norb,is, nspin);
			Complex val = diagElements[k0];
			Ham.AddEntry(k0,k0,val);
			onsiteAv+=std::norm(val);
		}
		std::cout<<"The average onsite is "<<onsiteAv/OpDim<<std::endl;
	}
	

	std::cout<<"Checking the hermiticity of the matrix"<<std::endl;
	bool is_hermitian=true;
	for(int i0=0;i0<n0   ; i0++ )
	for(int i1=0;i1<n1   ; i1++ )
	for(int io=0;io<NORB ; io++ )
	for(int is=0;is<SPIN ; is++ )
	{
		const int 
		row= IndexesToIndex(i0,n0,i1,n1,io,norb,is,nspin),
		coordnum0=Ham.Op_triplet[row].size();

		for(int n=0; n<coordnum0 ; n++)
		{
			const int 
			col= Ham.Op_triplet[row][n].col;
			const Complex 
			val= Ham.Op_triplet[row][n].val;

			const int 
			coordnum1=Ham.Op_triplet[col].size();
			for(int m=0; m<coordnum1 ; m++)
				if ( Ham.Op_triplet[col][m].col == row )
					if ( std::norm( std::conj(Ham.Op_triplet[col][m].val) -val) > tol_zero )
						std::cout<<row<<" "<<col<<" "<<val<<" | "<<col<<" "<<Ham.Op_triplet[col][m].col<<" "<<  Ham.Op_triplet[col][m].val<<std::endl;
			 

		}
		
	}

	Ham.saveTxtCoo( label+".Ham" );

	Ham.WriteIntoFile(label+".Ham");
	VelX.WriteIntoFile(label+".Vx");
	VelY.WriteIntoFile(label+".Vy");
	VelXSz.WriteIntoFile(label+".VxSz");
	VelYSz.WriteIntoFile(label+".VySz");

return 0;}



	
