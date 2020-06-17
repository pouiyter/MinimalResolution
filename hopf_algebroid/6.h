//resolution with a model
template<typename ring, typename algebroid>
template<typename degree_type, typename table_type>
cofree_comodule<algebroid,degree_type> Hopf_Algebroid<ring,algebroid>::resolvor_modeled(comodule_generic<algebroid,degree_type> &X, matrix<ring> *inj, matrix<ring> *quot, matrix<ring> *indj, matrix<ring> *new_map, curtis_table<table_type> *table, std::vector<int> *gens, std::vector<int> *nex_gens, std::function<vectors<matrix_index,ring>(vectors<matrix_index,table_type> const&)> transformer, std::iostream &maps_file){
	//embed using the model
	std::cout << "embedding...\n" << std::flush;
	auto F = embed2cofree_modeled(&X,inj,gens);

	//prepare the gaussian elimination
	std::cout << "Preparing Gausing...\n" << std::flush;
	indj->clear();
	//compute the cycles
	auto inj_ind = table->cycle_matrix(*indj,X.rank(),transformer,inj);
	auto quot_inds = matrix<ring>::quot_index(inj_ind,F.rank());
	//get the irrelevant indices
	std::set<int> to_del;
	for(int j : quot_inds.first)
		to_del.insert(j);
	for(int k: *nex_gens)
		to_del.erase(quot_inds.first[k]);
	//get the collumns which form an invertible sub
	std::vector<std::pair<matrix_index,matrix_index>> gs;
	for(unsigned i=0; i<inj_ind.size(); ++i)
		gs.push_back(std::make_pair((matrix_index) i, (matrix_index) inj_ind[i]));
	sort_deg(gs,X.base_module.degree);
	
	//Gaussian ellimiation
	std::cout << "Gausings...\n" << std::flush;
	indj->del_and_gaussian(gs,to_del);
	
	//compute the quotient comodule
	std::cout << "quotieting...\n" << std::flush;
	indj->make_quotient(quot_inds.first, quot_inds.second, inj_ind, quot);
	quotient(&F, quot, quot_inds.first, X);
	return F;
}

//get the short exact sequences using a model
template<typename ring, typename algebroid>
template<typename degree_type, typename table_type>
void Hopf_Algebroid<ring,algebroid>::pre_resolution_modeled( comodule_generic<algebroid,degree_type>& resolved, string filename_maps, string filename_generators, int resolution_length, string filename_table, curtis_table<table_type> *table, std::vector<std::vector<int>> &gens,std::function<vectors<matrix_index,ring>(vectors<matrix_index,table_type> const&)> transformer, matrix<ring> *inj, matrix<ring> *qut, matrix<ring> *indj, matrix<ring> *new_map, string back_up_file_name, int res_start){
	maps_file_name = filename_maps;
	
	//open the files for the generators and maps
	auto mode = res_start>=0 ? std::fstream::app : std::fstream::out;
	std::fstream maps_file(filename_maps, mode | std::fstream::binary);
	std::fstream gens_file0(filename_generators, mode | std::fstream::binary);
	std::fstream table_file(filename_table, std::ios::in |std::fstream::binary);
	
	if(!table_file.is_open())
		std::cerr << "fail to open" << filename_table;
	
	int32_t M_rank;
	//load the partial resolution
	if(res_start>=0){
		std::fstream back_file(back_up_file_name + std::to_string(res_start), std::fstream::in);
		resolved.coaction_matrix->load(back_file);
		back_file.close();
	}
	
	//get the rank of the starting comodule
	M_rank = resolved.rank();  
	if(res_start<0){
		gens_file0.write((char*)&M_rank, 4);
	}
	gens_file0.close();
	maps_file.close();
	
	for(int i=0; i<=resolution_length; ++i) {
		inj->clear();
		qut->clear();
		
		//open the file for the maps
		std::fstream maps_file(filename_maps + std::to_string(i), std::ios::out | std::fstream::binary);
		
		//do resolution
		std::cout << "ready to resolve" << i << "/" << resolution_length << "\n" << std::flush;
		table->load(table_file);
 //		std::cout << table->output() << std::flush;
		auto F = resolvor_modeled(resolved, inj, qut, indj, new_map, table, &gens[i], &gens[i+1],transformer, maps_file);
		
		//save the data
		std::fstream gens_file(filename_generators + std::to_string(i), std::ios::out | std::fstream::binary);

		//save the quotient to the new comodule
		inj->save(maps_file);
		qut->save(maps_file);
		//save generators
		F.save(gens_file);
		
		//update the rank
		M_rank = resolved.rank();  
		gens_file.write((char*)&M_rank, 4);
		
		maps_file.close();
		gens_file.close();
		
		//write the back-up file
		if(back_up_file_name!=""){
			std::fstream back_file(back_up_file_name + std::to_string(i), std::fstream::out | std::fstream::binary);
			
			std::cout << back_file.is_open();
		//	back_file.write((char*)&M_rank, 4);
			resolved.coaction_matrix->save(back_file);
			back_file.close();
			
			
			back_file.open(back_up_file_name + std::to_string(i), std::fstream::in | std::fstream::binary);
			
			std::cout << back_file.is_open();
		//	back_file.read((char*)&M_rank, 4);
			resolved.coaction_matrix->load(back_file);
			M_rank = resolved.coaction_matrix->rank;
			back_file.close();
		}
	}
	
}

//combing the data of generators into one file
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::gens_file_combiner(string filename_generators, int resolution_length, comodule_generic<algebroid,degree_type>&){
	std::fstream gens_file(filename_generators, std::ios::in | std::ios::binary);

	int32_t M_rank;
	gens_file.read((char*)&M_rank,4);
	gens_file.close();
	
	gens_file.open(filename_generators, std::ios::out | std::ios::trunc |std::ios::binary);
	gens_file.write((char*)&M_rank,4);

	for(int i=0; i<=resolution_length; ++i){
		std::fstream f(filename_generators + std::to_string(i), std::ios::in | std::ios::binary);
		//load the generators
		cofree_comodule<algebroid, degree_type> F;
		F.load(f);
		f.read((char*)&M_rank,4);
		
		//save the generator in the new file
		F.save(gens_file);
		gens_file.write((char*)&M_rank,4);
	}
}

//load the data for generators
template<typename ring, typename algebroid>
void Hopf_Algebroid<ring,algebroid>::load_gens(std::vector<std::vector<int>> &gens, string gens_data){
	std::fstream genfile(gens_data, std::ios::in | std::ios::binary);
	std::cout << gens_data;
	std::cout << genfile.is_open() << std::flush;
	int32_t sz;
	genfile.read((char*)&sz, 4);
	gens.resize(sz);
	for(int i=0; i<sz; ++i){
		int32_t ss;
		genfile.read((char*)&ss, 4);
		gens[i].resize(ss);
		for(int j=0; j<ss; ++j){
			int32_t as;
			genfile.read((char*)&as, 4);
			gens[i][j] = as;
		}
	}
}