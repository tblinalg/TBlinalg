void tblinalg::TBOperator::SetRowRange(const std::string ddFilename)  
{
	int rank, worldSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	//Cheack if the file existas and open the file
	std::string
	localFilename= ddFilename+".DD";	
	if ( ! fileExists(localFilename) )
	{
		std::cerr<<"File: "<<localFilename<<" not found, return -1;";
		std::exit(-1);
	}
	std::ifstream OpFile  (localFilename.c_str(), std::ofstream::binary);
	//Read the dimension of the matrix, and the total number 
	size_t  Dom[3];
	OpFile>>Dom[0]>>Dom[1]>>Dom[2];

	size_t rowOrig,rowEnd;
	OpFile>>rowOrig;
	for ( size_t r=0; r<worldSize ;r++ )
	{
		OpFile>>rowEnd;
		if (r==rank)
		{
			SetRowOrigin(rowOrig);
			SetRowEnd(rowEnd);
			SetRowDim(rowEnd-rowOrig);		
			OpFile.close();
			return ;
		}
	rowOrig=rowEnd;
	}

};	
