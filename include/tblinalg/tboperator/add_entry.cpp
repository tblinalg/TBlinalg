void tblinalg::TBOperator::AddEntry(const size_t row,const size_t col,const complex val )
	{
		MatRow matRow = matGrid.Block().GetMatRow(row);

		size_t rowSize= matRow.Dim() ;
		bool col_found = false;

		if ( rowSize == 0)
		{ 
			matRow.AddMatTuple(col,val);
			col_found = true;
		}
		else for(size_t elIdx=0; elIdx < rowSize ;elIdx++)
		{
			MatTuple elem = matRow.GetMatTuple(elIdx) ;
			if( elem.col == col )
			{
				elem.val+= val ;
				elIdx=rowSize;
				col_found = true;
			}
		}
		if ( !col_found )
			matRow.AddMatTuple(col,val);
		
		if( matRow.Dim() > BufferSize() )
		{
		
			std::cerr	<<std::endl<<"Error:"
						<<"Attempt to add element in row: "<< row<<std::endl
						<<"with dimension "<<matRow.Dim()<<std::endl		
						<<"exceed the buffer size: "<< BufferSize()<<std::endl
						<<"Choose another "
						<<"buffer size or verify the hamiltonian"
						<<std::endl;
			std::cerr<< "The row has the following elements"<<std::endl;
			for(size_t elIdx=0; elIdx < matRow.Dim() ;elIdx++)
			{
				MatTuple elem = matRow.GetMatTuple(elIdx) ;
				std::cerr<<row<<" "<<elem.col<<" "<<elem.val<<std::endl;
			}
			std::cerr<<std::endl;
			
			std::exit(-1);
		}

		
	matGrid.Block().GetMatRow(row)	= 	matRow ;
	};
