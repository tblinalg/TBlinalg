size_t tblinalg::TBOperator::NumNonZero( ) 
	{
		const size_t dim = Dim();

		if( numNonZero_ == 0)
		for(size_t row=0; row < dim ; row++)
		{
			const MatRow matrow  = matGrid.Block().GetMatRow(row);
			const size_t numElems=matrow.Dim();
			
			for(size_t elIdx=0;	elIdx < numElems ;  elIdx++)
			if( std::norm( matrow.GetMatTuple(elIdx).val ) > RealZero )
					numNonZero_+=1;
		}
		return numNonZero_;


	}
