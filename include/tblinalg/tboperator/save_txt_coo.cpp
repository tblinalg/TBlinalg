

void tblinalg::TBOperator::saveTxtCoo(std::string label )
{
	const size_t dim = Dim();

	std::string opFilename = label + ".COO" ;
	std::ofstream opFile( opFilename.c_str()  );
		

	for(size_t row=0;row < dim; row++)
	{
		MatRow unsortedRow;
		const MatRow matrow = matGrid.Block().GetMatRow(row);
		const size_t rowSize = matrow.Dim() ;
		
		for(size_t elIdx=0;elIdx < rowSize ; elIdx++)
		{
			const MatTuple tuple = matrow.GetMatTuple(elIdx);
			const size_t col = tuple.col;
			const complex val = tuple.val;
			if( std::norm( val ) > RealZero )
			unsortedRow.AddMatTuple(col,val);
		}

	unsortedRow.SortRow();
	for(size_t elIdx=0;elIdx < unsortedRow.Dim() ; elIdx++)
	{
		size_t  col =unsortedRow.GetMatTuple(elIdx).col;
		complex val =unsortedRow.GetMatTuple(elIdx).val;
		opFile	<<row<<" "
		<<col<<" "
		<<val.real()<<" "
		<<val.imag()<<std::endl;
	}
	}	
opFile.close();
};
