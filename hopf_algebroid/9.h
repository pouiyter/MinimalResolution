// adjoint to the projection to the n-th coordinate
template<typename ring, typename algebroid>
template<typename degree_type>
cofree_comodule<algebroid,degree_type> Hopf_Algebroid<ring,algebroid>::adjoint(const CoModule<algebroid,degree_type> *X, int n, matrix<ring> *adjoint_map,int shift){
	// construct a cofree comodule of rank 1
	cofree_comodule<algebroid,degree_type> result;
	result.generators.rank = 1;
	result.generators.degree.push_back(X->degree(n));
	result.position_of_gens.push_back(0);
	result.total_rank = ranksBelowDeg(maxDeg-result.underlyingDeg(X->degree(n)));
	//compute the adjoint map
	std::function<vectors<matrix_index,ring>(int)> ad_map_rows = [this, &X, n,shift] (int i) {
		return algebroid2vector(algebroidModuleOper->component(n,X->coaction(i)),shift);
	};
	
	if(adjoint_map!=NULL)
		adjoint_map->construct(X->rank(),ad_map_rows);

	return result;
}

// adjoint to projection on many coordinates, the i-th component
template<typename ring, typename algebroid>
template<typename degree_type>
vectors<matrix_index,ring> Hopf_Algebroid<ring,algebroid>::adjoint(const CoModule<algebroid,degree_type> *X, std::vector<int> const& gens, std::vector<uint32_t> const& pos, int i, int shift){
	//read the coaction
	auto v = X->coaction(i);
	vectors<matrix_index,ring> result;
	for(unsigned w=0; w<gens.size(); ++w){
		//tranform the comonent into an element in the cofree comodule
		auto ns = algebroid2vector(algebroidModuleOper->component(gens[w],v),pos[w]+shift); 
		result.direct_sum(ns,0);
	}
	return result;
}

