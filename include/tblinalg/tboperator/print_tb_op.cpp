void tblinalg::TBOperator::PrintTBOp(const integer dx,const integer dy) 
{
	for(size_t rowIdx= 0; rowIdx < RowDim() ; rowIdx++)
	for(size_t elem=0;elem < matGrid(dx,dy).GetMatRow(rowIdx).Dim() ; elem++)
	{
		size_t row= RowOrigin() + rowIdx ;
		integer col=matGrid(dx,dy).GetMatRow(rowIdx).GetMatTuple(elem).col;
		complex val=matGrid(dx,dy).GetMatRow(rowIdx).GetMatTuple(elem).val;
		std::cout<<row<<" "<<col<<" "<<val.real()<<" "<<val.imag()<<std::endl;
	}		
};
