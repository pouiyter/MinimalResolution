//compute the quotient comodule
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::quotient(const CoModule<algebroid,degree_type> *X, matrix<ring> *quot, std::vector<matrix_index> inverse_ind, comodule_generic<algebroid,degree_type> &result){
	//get the rank of the quotient comodule
	result.base_module.rank = inverse_ind.size();

	//compute the degrees of the generators of the quotient
	result.base_module.degree.resize(inverse_ind.size());
	for(unsigned i=0; i<inverse_ind.size(); ++i)
		result.base_module.degree[i] = X->degree(inverse_ind[i]);
	std::function<vectors<matrix_index,algebroid>(const algebroid&, const vectors<matrix_index,ring>&)> right_mult = [this] (const algebroid& A, const vectors<matrix_index,ring> &V){
		return right_scalor_mult(A,V); };
	
	//compose the quotien matrix and the coaction matrix. Note that we need to first apply the right unit to the quotient matrix
	std::function<vectors<matrix_index,algebroid>(int)> coactor = [X,quot,&inverse_ind,right_mult,this] (int i){
		vectors<matrix_index,algebroid> coa = X->coaction(inverse_ind[i]);
		AbGroupOp<vectors<matrix_index,algebroid>> *adop = algebroidModuleOper;
		std::cout << "\r" << i << "/" << inverse_ind.size() << std::flush;
		return quot->maps_to(coa,right_mult, adop); 
	};
	std::cout << "computing the coactions...\n" << std::flush;
	result.coaction_matrix->construct(inverse_ind.size(), coactor);
}

#include"matrices_mem.h"
//compute the quotient comodule
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::quotient_p(const CoModule<algebroid,degree_type> *X, matrix<ring> *quot0, std::vector<matrix_index> inverse_ind, comodule_generic<algebroid,degree_type> &result){
	//load the quot matrix
	matrix_mem<ring> quot1;
	matrix<ring> *quot = &quot1;
	quot->construct(quot0);
	
	//get the rank of the quotient comodule
	result.base_module.rank = inverse_ind.size();
	
	//compute the degrees of the generators of the quotient
	result.base_module.degree.resize(inverse_ind.size());
	for(unsigned i=0; i<inverse_ind.size(); ++i)
		result.base_module.degree[i] = X->degree(inverse_ind[i]);
	std::function<vectors<matrix_index,algebroid>(const algebroid&, const vectors<matrix_index,ring>&)> right_mult = [this] (const algebroid& A, const vectors<matrix_index,ring> &V){
		return right_scalor_mult(A,V); };
		
		//compose the quotien matrix and the coaction matrix. Note that we need to first apply the right unit to the quotient matrix
		std::function<vectors<matrix_index,algebroid>(int)> coactor = [X,quot,&inverse_ind,right_mult,this] (int i){
			vectors<matrix_index,algebroid> coa = X->coaction(inverse_ind[i]);
			AbGroupOp<vectors<matrix_index,algebroid>> *adop = algebroidModuleOper;
			std::cout << "\r" << i << "/" << inverse_ind.size() << std::flush;
			return quot->maps_to_p(coa,right_mult, adop); 
		};
		std::cout << "computing the coactions...\n" << std::flush;
		result.coaction_matrix->construct_parallel(inverse_ind.size(), coactor);
}

//scalor multiplication using the right unit
template<typename ring, typename algebroid>
vectors<matrix_index, algebroid> Hopf_Algebroid<ring,algebroid>::right_scalor_mult (algebroid const &A, vectors<matrix_index, ring> const &V) {
	if(algebroidRingOper->isZero(A))
		return algebroidModuleOper->zero();
	std::function<algebroid(ring const&)> rule = [this,&A] (ring const &r) { return algebroidRingOper->multiply(A,etaR(r)); };
	return moduleOper->termwise_operation(rule, V);
}
