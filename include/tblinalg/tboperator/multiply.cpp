#define ROOT -1


void 
tblinalg::TBOperator::Multiply(	const integer dx,
								const integer dy,
								const real a, 
								const std::vector<complex>& X, 
								const real b, 
								std::vector<complex>& Y )
{
	const int rank = bIdx.ReduceToOne();

		if (rank == ROOT )
		{
			std::cout<<"Input X "<<std::endl;
			for(int r=0; r< RowDim() ; r++)
			if (X[ r  ].real()!=0 )
				std::cout<<r<<" "<<X[ r  ]<<std::endl;
		}

	for( size_t rowIdx=0; rowIdx < RowDim() ; rowIdx++ )
	{
		
		Y[rowIdx]= Y[rowIdx]*b;

		const size_t numElems=matGrid(dx,dy).GetMatRow(rowIdx).Dim();
		const MatRow matrow = matGrid(dx,dy).GetMatRow(rowIdx);	
		for( size_t elIdx=0; elIdx < numElems ; elIdx++ )
		{
			const MatTuple t	= matrow.GetMatTuple(elIdx);
			const integer colIdx= t.col-destColOrig;
			const complex val   = a*t.val ;

			if( rank == ROOT)
//			if ( colIdx>= 0 && colIdx< RowDim() )
//				if ( colIdx>= 0  && std::norm(X[colIdx])!=0 )
					std::cout//<<"Y("<<rowIdx+RowOrigin()<<") = "
							 //<<b<<"* Y("<<rowIdx+RowOrigin()<<") + "<<a
							 //<<"* H("<< rowIdx+RowOrigin()<<","<<t.col<<")"
							 //<<" X("<<t.col<<")"<<std::endl
							 <<"RowDim"<<RowDim()<<" "
							 <<"a="<<a<<" b="<<b<<" colOrigin= "<<destColOrig<<std::endl
							  <<" aH("<< rowIdx+RowOrigin()<<","<<t.col<<")="<< val.real()
							 <<" X("<<colIdx<<")="<<X[colIdx].real()
							 <<" bY("<<rowIdx+RowOrigin()<<")="<<Y[rowIdx].real()<<std::endl;		

			
			if ( colIdx>= 0 && colIdx< RowDim() )
				Y[rowIdx]+= X[colIdx]*val;

		}
	}
}

void 
tblinalg::TBOperator::Multiply(	const real a, 
								const std::vector<complex>& X, 
								const real b,
								std::vector<complex>& Y 
								) 
{
	const int rank = bIdx.ReduceToOne();
	//On-site multiplication
	if(rank == ROOT )
		std::cout<<" I am 0, I send to myself and receive from myself"<<std::endl;
	
	//Save the X vector to be send as source
	sourceX = X;
	sourceColOrig=RowOrigin();
	destColOrig  =sourceColOrig; //not remove!
	Multiply(0,0,a, X,b, Y );

	//if(rank == ROOT )
	//for( size_t rowIdx=0; rowIdx < RowDim() ; rowIdx++ )
	//{
	//	std::cout<<"Matrix Vector X: "<<X[rowIdx]<<" "<<X[rowIdx]*b<<" "<<b<<" "<<a<<std::endl;
	//	std::cout<<"Matrix Vector Y: "<<Y[rowIdx]<<" "<<Y[rowIdx]*b<<" "<<b<<" "<<a<<std::endl;
	//}
	
	//Off-site multiplication
	for(integer dx=-1; dx <= 1 ; dx++)
	for(integer dy=-1; dy <= 1 ; dy++)
	if ( dx!=0 || dy !=0 ) 
	{
		int rectag=1,sendtag=1, source, dest, sourceBuffSize=RowDim(), destBuffSize;
		int sourceCoord[3] ={bIdx(0)+dx,bIdx(1)+dy,0};
		int destCoord[3]   ={bIdx(0)-dx,bIdx(1)-dy,0};

		MPI_Cart_rank(cartComm,sourceCoord,&source);	
		MPI_Cart_rank(cartComm,destCoord,&dest);

		MPI_Cart_coords(cartComm, source,2,sourceCoord);	
		MPI_Cart_coords(cartComm, dest,2,destCoord);

		MPI_Sendrecv(	&sourceBuffSize ,1 , MPI_INT,dest,sendtag,
						&destBuffSize   ,1 , MPI_INT,source,rectag,
						cartComm,MPI_STATUS_IGNORE );
		destX = std::vector<complex>(destBuffSize);
						
		MPI_Sendrecv(	&sourceColOrig ,1 , MPI_INT,dest,sendtag,
						&destColOrig   ,1 , MPI_INT,source,rectag,
						cartComm,MPI_STATUS_IGNORE );
		MPI_Sendrecv(	&sourceX[0] ,sourceBuffSize, MPI_DOUBLE_COMPLEX,dest,sendtag,
						&destX[0]   ,destBuffSize, MPI_DOUBLE_COMPLEX,source,rectag, 
						cartComm,MPI_STATUS_IGNORE );

		if( rank == ROOT)
		{
			std::cout<<"I am at block: ( "<<bIdx(0)<<" , "<<bIdx(1)<<" )"<<" with "<<dx<<" "<<dy<<std::endl;
			std::cout<<"I send  to dest: "<<dest<<" with coord  ( "<<destCoord[0]<<" , "<<destCoord[1]<<" )"	<<std::endl;
			std::cout<<"The vector Xsource: ";
			for( size_t rowIdx=0; rowIdx < RowDim() ; rowIdx++ )
			if ( std::norm(sourceX[rowIdx] )!=0 )
				std::cout<<"( "<<rowIdx+RowOrigin()<<" , "<<sourceX[rowIdx].real()<<" ) ";
			std::cout<<std::endl;
	
			std::cout<<"and receive from source: "<<source<<" with coord  ( "<<sourceCoord[0]<<" , "<<sourceCoord[1]<<" )"	<<std::endl;
			std::cout<<"The vector Xdest: ";
			for( size_t rowIdx=0; rowIdx < RowDim() ; rowIdx++ )
			if ( std::norm(destX[rowIdx] )!=0 )
				std::cout<<"( "<<rowIdx+RowOrigin()<<" , "<<destX[rowIdx].real()<<" ) ";
			std::cout<<std::endl;

		}



		Multiply(dx,dy,a, destX,1, Y );
	}

};	
