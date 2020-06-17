//set a matrix to the unit matrix
template<typename R>
void matrix<R>::set2unit(int rk){
	std::function<vectors<matrix_index,R>(int)> rows = [this](int i){
		return this->moduleOper->singleton(i); };
	construct(rk, rows);
}

//set a matrix to zero
template<typename R>
void matrix<R>::set2zero(int rk){
	std::function<vectors<matrix_index,R>(int)> rows = [this](int i){ 
		return this->moduleOper->zero(); };
	construct(rk,rows);
}

//check if two matrices are equal
template<typename R>
bool matrix<R>::equal(matrix<R>& Y){
	if(rank != Y.rank)
		return false;
	for(unsigned i=0; i<rank; ++i){
		auto diff = this->moduleOper->add(this->moduleOper->minus(find(i)), Y.find(i));
		if(!this->moduleOper->isZero(diff))
			return false;
	}
	return true;
}
