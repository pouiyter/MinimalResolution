//return the leading term index for a quotient map. inj_index lists the leading indices of an injective map. Return the indices left out in the quotient space and the leading term correspondence for the quotient map
template<typename R>
std::pair<std::vector<matrix_index>,std::map<matrix_index,matrix_index>> matrix<R>::quot_index(std::vector<matrix_index> const  &inj_index, unsigned rank){
	std::vector<matrix_index> inverse_ind;
	std::map<matrix_index,matrix_index> quot_ind;
	
	//mark those indices missed by the subspace to be quotiented out
	std::vector<bool> missed(rank,true);
	for(auto i:inj_index) 
		missed[i]=false;
	
	int cur_ind = 0;
	for(int i=0; i<(int)rank; ++i)
		if(missed[i]){
			//if i is not the leading term of the subsapce, then it should appear in the quotient
			inverse_ind.push_back(i);
			
			//mark the index of i in the quotient space
			quot_ind.emplace(i,cur_ind);
			++cur_ind;
		}
	return std::make_pair(std::move(inverse_ind), std::move(quot_ind));
}

//make a quotient matrix. This should be allready in an echelon form. This is essentially tansform (I,A) into the transpose of (-A,I)
template<typename R>
void matrix<R>::make_quotient(std::vector<matrix_index> const &inverse_ind, std::map<matrix_index,matrix_index> const &quot_ind, std::vector<matrix_index> const &inj_index, matrix<R> *result){
	result->clear();
	
	//find the quotient map from a matrix in echelon form. The indices to be quotiented out has to be transformed into the other ones
	for(unsigned i=0; i<inj_index.size(); ++i) {
		auto newTerm = moduleOper->minus(find(i));
		newTerm = moduleOper->add(newTerm, moduleOper->singleton(inj_index[i]));
		
		std::function<matrix_index(matrix_index)>  rule = [&quot_ind] (matrix_index n){ 
			return quot_ind.at(n); };
		result -> insert(inj_index.at(i), moduleOper -> re_index(rule, newTerm));
	}
	//on the remaining indices are identities 
	for(unsigned i=0; i<inverse_ind.size(); i++)
		result -> insert(inverse_ind[i], moduleOper -> singleton(i));
}
