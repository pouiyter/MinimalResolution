//operations on each rows
template<typename R>
void matrix<R>::row_operation(std::function<vectors<matrix_index,R>(vectors<matrix_index,R> const&,int)> row_rule, matrix<R> *result){
	std::function<vectors<matrix_index,R>(int)> rows = [this, &row_rule](int i){
		return row_rule(find(i),i); };
	result->construct(rank, rows);
}

//filter each row vector
template<typename R>
void matrix<R>::filter(std::function<bool(matrix_index)> rule, matrix<R> *result){
	std::function<vectors<matrix_index,R>(int)> mk = [this,&rule](int i){ 
		return moduleOper->filter(rule, find(i)); };
	result->construct(rank, mk);
}

//filter each row vector
template<typename R>
void matrix<R>::filter(std::function<bool(matrix_index,const R&)> rule, matrix<R> *result){
	std::function<vectors<matrix_index,R>(int)> mk = [this,&rule](int i){ 
		return moduleOper->filter(rule, find(i)); };
	result->construct(rank, mk);
}
