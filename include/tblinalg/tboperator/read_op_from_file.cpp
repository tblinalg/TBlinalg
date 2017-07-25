void tblinalg::TBOperator::ReadOpFromFile(const int dx, const int dy, const std::string opFilename)
{
	//Cheack if the file existas and open the file
	std::string
	localFilename= opFilename+".OP";	
	
	if ( ! fileExists(localFilename) )
	{
		std::cerr<<"File: "<<localFilename<<" not found, return -1;";
		std::exit(-1);
	}
	std::ifstream OpFile  (localFilename.c_str(), std::ofstream::binary);
	
	//Read the dimension of the matrix, and the total number 
	size_t  OpDim, OpNumEntries;
	OpFile>>OpDim>>OpNumEntries;
	SetTotNumOrbs(OpDim);
	//Set initialize the current (x,y) element of the grid
	//to its corresponding site
	if( RowDim()== 0 )
	{
		SetRowDim(OpDim);
		std::cout<<"Warning, OpDim() not setted. Using: "<<OpDim<<std::endl;
	}
	
	matGrid(dx,dy).SetDim( RowDim() );
	
	for( size_t row = 0; row < TotNumOrbs() ; ++row )
	{
		bool isMyRow= (row >= RowOrigin() && row < RowEnd() ) ;
		//Get the number of elements per row
		size_t ElemsPerRow;
		OpFile>>ElemsPerRow;
		
		//Create a row to store the elements 
		MatRow matrow;
		
		//Iterate over the elements
		for( size_t elem=0; elem < ElemsPerRow ; ++elem )
		{
			size_t col; real Reval, Imval;
			OpFile>>col>>Reval>>Imval;

			if( isMyRow )
			matrow.AddMatTuple(col,complex(Reval,Imval) );
		}
		
		if( isMyRow)
		matGrid(dx,dy).SetMatRow(row-RowOrigin(),matrow.SortRow() );	
	}
	OpFile.close();
};	
