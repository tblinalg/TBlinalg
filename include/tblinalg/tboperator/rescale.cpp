void tblinalg::TBOperator::Rescale(	const size_t dx, const size_t dy , 
									const tblinalg::real shift,	
									const tblinalg::real scal )
{
	for( size_t rowIdx= 0; rowIdx <RowDim() ; rowIdx++ )
	{
		const size_t row = rowIdx + RowOrigin() ;
		const size_t numElems = matGrid(dx,dy).GetMatRow(rowIdx).Dim();
		for( size_t elIdx=0; elIdx < numElems; elIdx++ )
		{
			MatTuple
			t = matGrid(dx,dy).GetMatRow(rowIdx).GetMatTuple(elIdx);
			if ( row  == t.col )
					t.val -= shift ;
				t.val*= scal ;			
				//Return to the tuple
				matGrid(dx,dy).GetMatRow(rowIdx).GetMatTuple(elIdx).val =t.val ;
		}
	}
};
