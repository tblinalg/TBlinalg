
void tblinalg::TBOperator::WriteIntoFile(std::string label )
{
	const size_t dim = Dim();
	std::string opFilename = label + ".OP" ;
	std::ofstream opFile( opFilename.c_str(),  std::ofstream::binary );

	
		opFile	<<Dim()<<" "
				<<NumNonZero( ) <<" ";

	for(size_t row=0; row < dim ; row++)
	{

		const MatRow matrow = matGrid.Block().GetMatRow(row);
		const size_t numElems = matrow.Dim() ;
		
		opFile<<NumNonZeroInRow( row )<<" ";
		MatRow unsortedRow;
		
		for(size_t elIdx=0;elIdx < numElems ; elIdx++)
		{
			const MatTuple t = matrow.GetMatTuple(elIdx);
			if( std::norm( t.val ) > RealZero )
				unsortedRow.AddMatTuple(t.col,t.val);
		}
		
		unsortedRow.SortRow();
		for(size_t elIdx=0;elIdx < unsortedRow.Dim() ; elIdx++)
		{
			const MatTuple t = unsortedRow.GetMatTuple(elIdx) ;
			opFile	<<t.col<<" "
			<<t.val.real()<<" "
			<<t.val.imag()<<" ";
		}
	}			

	opFile.close();

};
