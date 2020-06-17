//embed into a cofree comodule
template<typename ring, typename algebroid>
template<typename degree_type>
cofree_comodule<algebroid,degree_type> Hopf_Algebroid<ring,algebroid>::embed2cofree(const CoModule<algebroid,degree_type> *X, matrix<ring> *inj, curtis_table<ring> *table, std::vector<int> *gens, matrix<ring> *new_map, std::vector<int> *basis_order, std::iostream &tablefile){
	std::cout << X->rank() << std::flush;
	std::vector<int> gns;
	if(gens==NULL) gens = &gns;
	
	//initialize the map to a zero comodule
	inj->set2zero(X->rank());
	table->clear();
	cofree_comodule<algebroid,degree_type> res;
//	std::cout << "KK" << std::flush;
	//the set of basis
	std::vector<int> bas;
	for(int i=0; i<X->rank(); ++i)
		bas.push_back(i);
	
	//sort the basis under degree
	auto i_deg = [X] (int m){ 
		return cofree_comodule<algebroid,degree_type>::underlyingDeg(X->degree(m)); };
	std::cout << "sorting..." << std::flush;
	std::stable_sort(bas.begin(),bas.end(),[i_deg](int m,int n) { return i_deg(m)<i_deg(n); });
	
	if(basis_order != NULL)
		*basis_order = bas;
	
//	int current_deg = 0;
	for(int j=0;j<X->rank();++j){
		std::cout << "\r" << j << "/" << X->rank() << std::flush;
	//	std::cout << j << "PP" << std::flush;
		//look at the i-th row of the current map
		int i = bas[j];
			
		//load the i-th row
		std::cout << "i" << std::flush;
		vectors<matrix_index, ring> irow =  adjoint(X, *gens, res.position_of_gens, i, 0);
		//trace what is inside the current vector
		std::cout << "m" << std::flush;
		vectors<matrix_index, ring> sc = moduleOper->singleton(i);
		//simplify the i-th row using the current table, saving the homotopy in sc
		std::cout << res.rank() << "s" << X->rank() << std::flush;
		matrix_index pos = table->simplify_to_led(res.rank(),X->rank(),irow,sc);
	//	matrix_index pos = table->symplify_to_led(irow,sc);

		if(pos!=curtis_table<ring>::Boundary){
			//this means the i-th row maps to a nontrivial element in the cofree one
			std::cout << "e" << std::flush;
			table->insert(pos,i,irow,sc);
			//update inj
			std::cout << "u" << std::flush;
			inj->insert(i,irow);
		}
		else{
			//we need a new generator
			std::cout << "n" << std::flush;
			gens->push_back(i);
			
			//record the current rank
			int old_rank = res.rank();
			
			//compute the new summand of the cofree comodule
			std::cout << "n" << std::flush;
			cofree_comodule<algebroid,degree_type> new_sumd = adjoint(X,i,NULL,old_rank);

			//add the new summand
			std::cout << "d" << std::flush;
			res.direct_sum(new_sumd);

			//add new summand to the maps
			std::cout << "a" << std::flush;
			irow = adjoint(X, *gens, res.position_of_gens, i, 0);
			std::cout << "j" << std::flush;
			inj->insert(i,irow);
			
			//update the table
			std::cout << "t" << std::flush;
			table->insert(old_rank,i,irow,sc);
		}
	}
//	table->save(tablefile);
	return res;
}

//embed into a cofree one, using a model
template<typename ring, typename algebroid>
template<typename degree_type>
cofree_comodule<algebroid,degree_type> Hopf_Algebroid<ring,algebroid>::embed2cofree_modeled(const CoModule<algebroid,degree_type> *X, matrix<ring> *inj, std::vector<int> *gens, matrix<ring>* new_map){
	//initialize
	std::cout << "a" << std::flush;
	inj->set2zero(X->rank());
	std::cout << "o" << std::flush;
	cofree_comodule<algebroid,degree_type> res;
	new_map->clear();
	cofree_comodule<algebroid,degree_type> new_sumd; 

	for(int i = 0; i<(int)gens->size(); ++i){
		std::cout << i << std::flush;
		int old_rank = res.rank();
		//get the new summand
		new_sumd = adjoint(X,(*gens)[i],new_map,0);
		//insert the new summand
		res.direct_sum(new_sumd);
		inj->direct_sum(*new_map,old_rank);
		new_map->clear();
	}
	return res;
}

//embed into a cofree one, using a model
template<typename ring, typename algebroid>
template<typename degree_type>
cofree_comodule<algebroid,degree_type> Hopf_Algebroid<ring,algebroid>::embed2cofree_modeled(const CoModule<algebroid,degree_type> *X, matrix<ring> *inj, std::vector<int> *gens){
	//initialize
	cofree_comodule<algebroid,degree_type> res;
	cofree_comodule<algebroid,degree_type> new_sumd; 
	std::cout << "a" << gens->size() << std::flush;
	//construct the new comodule
	for(int i = 0; i<(int)gens->size(); ++i){
		std::cout << i << "   " << std::flush;
		//get the new summand
		new_sumd = adjoint(X, (*gens)[i], NULL,0);
		//insert the new summand
		res.direct_sum(new_sumd);
	}
	
	//construct the new map
	std::function<vectors<matrix_index,ring>(int)> map_rows = [this, X, &res, gens](int i){
		std::cout << "\r" << i << std::flush;
		return adjoint(X, *gens, res.position_of_gens, i, 0);
	};
	inj->construct(X->rank(), map_rows);
	
	return res;
}
