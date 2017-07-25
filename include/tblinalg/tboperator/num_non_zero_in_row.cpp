size_t tblinalg::TBOperator::NumNonZeroInRow(size_t row )
{
	size_t nnz=0;
	const MatRow matrow=matGrid.Block().GetMatRow(row);
	const size_t numElems = matrow.Dim() ;

	for( size_t elIdx=0; elIdx < numElems ; elIdx++)
	if( std::norm( matrow.GetMatTuple(elIdx).val ) > RealZero )
		nnz+=1;
	return nnz;
	return 5;
}
