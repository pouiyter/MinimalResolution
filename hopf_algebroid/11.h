//make a multiplication table
template<typename ring, typename algebroid>
void Hopf_Algebroid<ring,algebroid>::make_multiplication_table(algebroid const &x, int deg_x, matrix<ring> *result) {
	result->clear();
	int n_entries = ranksBelowDeg(maxDeg-deg_x);
	for(int i=0; i<n_entries; ++i) {
		//get the i-th element
		auto h = vector2algebroid(moduleOper->singleton(i));
		//compute the multiplication and transform back to a vector in a cofree comodule
		auto y = algebroidRingOper->multiply(x, h);
		result->insert(i, algebroid2vector(y,0));
	}
}

//using the multiplication table to compute the multiplicative structure
template<typename algebroid, typename degree_type>
template<typename ring>
vectors<matrix_index,ring> cofree_comodule<algebroid,degree_type>::multiply_using_table(unsigned i, matrix<ring> const &mult_table){
	unsigned p = findPos(i);
	int k =  i - position_of_gens[p];
	auto v = mult_table.find(k);
	
	std::function<matrix_index(matrix_index)> shifting = [p,this] (matrix_index j) {
		return j + this->position_of_gens[p]; };
	return mult_table.moduleOper->re_index(shifting, v);
}
