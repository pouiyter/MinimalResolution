//construct a trivial comodule
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::set_to_trivial(comodule_generic<algebroid,degree_type> &X, degree_type deg) {
	X.base_module.rank = 1;
	X.base_module.degree.resize(1);
	X.base_module.degree[0] = deg;
	
	std::function<vectors<matrix_index,algebroid>(int)> rw = [this] (int i){
		return algebroidModuleOper->singleton(0,algebroidRingOper->unit(1)); };
	X.coaction_matrix->construct(1,rw);
}

//initializing the coaction of a cofree comodule
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::init_cofree_data(std::function<degree_type(matrix_index)> cofree_degree){
	//compute the coaction
	cofree_comodule<algebroid,degree_type>::cofree_coaction = [this] (matrix_index n){ 
		return delta(n); };

	//initialize the other data
	cofree_comodule<algebroid,degree_type>::modOpers = algebroidModuleOper;
	cofree_comodule<algebroid,degree_type>::cofree_degree = cofree_degree;
}

//the stable members of cofree comodules
template<typename algebroid, typename degree_type>
std::function<vectors<matrix_index, algebroid>(matrix_index)> cofree_comodule<algebroid,degree_type>::cofree_coaction;
template<typename algebroid, typename degree_type>
ModuleOp<matrix_index, algebroid> *cofree_comodule<algebroid,degree_type>::modOpers;
template<typename algebroid, typename degree_type>
std::function<degree_type(matrix_index)> cofree_comodule<algebroid,degree_type>::cofree_degree;
