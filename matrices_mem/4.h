//unify a vector to make a given entry 1
template<typename ring>
vectors<matrix_index,ring> matrix_mem<ring>::unify(matrix_index n, vectors<matrix_index,ring> const &x){
	ring po = this->moduleOper->ringOper->inverse(this->moduleOper->component(n,x));
	return this->moduleOper->scalor_mult(po,x);
}

//ellimination a given entry
template<typename ring>
void matrix_mem<ring>::take_away(matrix_index n, vectors<matrix_index,ring> &x, vectors<matrix_index,ring> const &y){
	ring cp = this->moduleOper->ringOper->minus(this->moduleOper->component(n,x));
	auto smd = this->moduleOper->scalor_mult(cp,y);
	x = this->moduleOper->add(std::move(x),std::move(smd));
}

// row transfermations, subtract the given row to cancel the given column
template<typename ring>
void matrix_mem<ring>::row_reduction(matrix_index row, matrix_index col){
	//unify the given row
	auto uv = unify(col,data[row]);

	#pragma omp parallel for schedule(dynamic)
	for(unsigned i=0; i<data.size(); ++i){
		//unify the given row
		if(i==row) data[i] = uv;
		//transform the other rows with vanhishing given column
		else take_away(col, data[i], uv);
	}
}

//Gaussian ellimination, with given rows and colums
template<typename ring>
void matrix_mem<ring>::gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols){
	int i = 0;
	for(auto rc: row_cols){
		std::cout << "\r" << ++i << std::flush;
		row_reduction(rc.first,rc.second);
	}
}
