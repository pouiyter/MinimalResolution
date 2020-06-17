//right multiplication on a row vector v
template<typename R>
vectors<matrix_index,R> matrix<R>::maps_to(const vectors<matrix_index,R> &v){
	//compute the scalor multiplication of the entries of v with corresponding rows of the matrix
	std::function<vectors<matrix_index,R>(int)> summands = [this, &v] (int i){
		return moduleOper->scalor_mult(v.dataArray[i].coeficient, find(v.dataArray[i].ind)); };
	//sum up
	return moduleOper->sum(summands,0,v.size());
}

//compose this with rows, with rank of rows equal to res_rank
template<typename R>
void matrix<R>::compose(std::function<vectors<matrix_index,R>(int)> rows, int res_rank, matrix<R> *result){
	std::function<vectors<matrix_index,R>(int)> res_rows = [this, rows, res_rank](int i){ 
		std::cout << "\r" << i << "/" << res_rank  << std::flush;
		return maps_to(rows(i)); };
	result->construct(res_rank, res_rows);
}

//compose with Y
template<typename R>
void matrix<R>::compose(const matrix<R> *Y, matrix<R> *result){
	std::function<vectors<matrix_index,R>(int)> rows = [Y](int i){ 
		return Y->find(i); };
	return compose(rows, Y->rank, result);
}

//compute the maps with informations on the generators
template<typename index, typename ring, typename module1,typename module2>
module2 maps_to(vectors<index,ring> const &x, std::function<module2(const ring&, const module1&)> scalor_mult, AbGroupOp<module2> *modop2, unsigned start, unsigned stop, std::function<module1(index)> find){
	if(stop<=start) return modop2->zero();
	if(stop-start==1) return scalor_mult(x.dataArray[start].coeficient, find(x.dataArray[start].ind));
	unsigned mid = (start+stop) / 2; 
	return modop2->add(maps_to(x,scalor_mult,modop2,start,mid,find), maps_to(x,scalor_mult,modop2,mid,stop,find));
}

//compute the maps with informations on the generators
template<typename index, typename ring, typename module1,typename module2>
module2 maps_to(vectors<index,ring> const &x, std::function<module2(const ring&, const module1&)> scalor_mult, AbGroupOp<module2> *modop2, std::function<module1(index)> find){
	return maps_to(x,scalor_mult,modop2,0,x.size(), find);
}

//compose with another matrix with different type of coeficients
template<typename R>
template<typename coef_type,typename module2>
module2 matrix<R>::maps_to(vectors<matrix_index,coef_type> const &x, std::function<module2(const coef_type&, const vectors<matrix_index,R>&)> scalor_mult, AbGroupOp<module2> *modop2){
	std::function<vectors<matrix_index,R>(matrix_index)> finder = [this](matrix_index i){
		return this->find(i); };
	return ::maps_to(x,scalor_mult,modop2,finder);
}

//compute the maps with informations on the generators
template<typename index, typename ring, typename module1,typename module2>
module2 maps_to_p(vectors<index,ring> const &x, std::function<module2(const ring&, const module1&)> scalor_mult, AbGroupOp<module2> *modop2, unsigned start, unsigned stop, std::function<module1(index)> find){
	if(stop<=start) return modop2->zero();
	if(stop-start==1) return scalor_mult(x.dataArray[start].coeficient, find(x.dataArray[start].ind));
	unsigned mid = (start+stop) / 2; 
	vectors<index,ring> s1, s2;
#pragma omp parallel sections
	{
#pragma omp section
		{ s1 = maps_to(x,scalor_mult,modop2,start,mid,find); }
#pragma omp section
		{ s2 = maps_to(x,scalor_mult,modop2,mid,stop,find); }
	}
	return modop2->add(std::move(s1), std::move(s2));
}

//compute the maps with informations on the generators
template<typename index, typename ring, typename module1,typename module2>
module2 maps_to_p(vectors<index,ring> const &x, std::function<module2(const ring&, const module1&)> scalor_mult, AbGroupOp<module2> *modop2, std::function<module1(index)> find){
	return maps_to_p(x,scalor_mult,modop2,0,x.size(), find);
}

//compose with another matrix with different type of coeficients
template<typename R>
template<typename coef_type,typename module2>
module2 matrix<R>::maps_to_p(vectors<matrix_index,coef_type> const &x, std::function<module2(const coef_type&, const vectors<matrix_index,R>&)> scalor_mult, AbGroupOp<module2> *modop2){
	std::function<vectors<matrix_index,R>(matrix_index)> finder = [this](matrix_index i){
		return this->find(i); };
		return ::maps_to_p(x,scalor_mult,modop2,finder);
}
