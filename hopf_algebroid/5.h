//embed a comodule into a cofree one and take the quotient
template<typename ring, typename algebroid>
template<typename degree_type>
cofree_comodule<algebroid,degree_type> Hopf_Algebroid<ring,algebroid>::resolvor(comodule_generic<algebroid,degree_type> &X, matrix<ring> *inj, matrix<ring> *indj, matrix<ring> *quot, curtis_table<ring> *table, matrix<ring> *new_map, std::vector<int> *gens, std::vector<int> *basis_order, std::iostream &tablefile){
	std::cout << "embedding...\n" << std::flush;
	//embed X into a cofree comodule F, with the injective map stored in inj, and the data for generators stored in gens
	auto F = embed2cofree(&X,inj,table,gens,new_map,basis_order,tablefile);

//	std::cout << inj->output();
	std::cout << "Preparing Gausing...\n" << std::flush;
	//prepare the injective map
	auto inj_ind = table->cycle_matrix(*indj,X.rank());
//	std::cout << table->output();
	table->save(tablefile);
//	table->clear();
	std::vector<std::pair<matrix_index,matrix_index>> gs;
	for(unsigned i=0; i<inj_ind.size(); ++i)
		gs.push_back(std::make_pair((matrix_index) i, (matrix_index) inj_ind[i]));

	std::cout << "Gausing...\n" << std::flush;

//	std::cout << indj->output();
	//do Gaussan ellimination on indj
	indj->gaussian(gs);
	
	std::cout << "quotieting...\n" << std::flush;
	//compute the quotient comodule, and replace X by this quotient
	auto quot_inds = matrix<ring>::quot_index(inj_ind,F.rank());

//	std::cout << indj->output();
//	for(auto k : quot_inds.first) std::cout << k << " " << std::flush;
//	for(auto k : quot_inds.second) std::cout << k.first << ";" << k.second << " " << std::flush;
	
	std::cout << "constructing the quotient matrix\n" << std::flush;
	indj->make_quotient(quot_inds.first, quot_inds.second, inj_ind, quot);

	indj->clear();
	std::cout << "constructing quotient coactions...\n" << std::flush;
	quotient_p(&F, quot, quot_inds.first, X);
	return F;
}

//construct the short exact sequences in a resolution
//resolved is the input comodule to be resolved. It will be the container for the intermediate steps. So it will be re-valued into the n-th "connected cover" of the orinigal one
//the filenames are the files to store the data of the resolution
//result will be the container to store the intermediate curtis tables to be used as the model for the lifting to other Hopf algebroids
//gens will be the data of the generators
//inj and qut will be the containers for the imtermidiet matrices
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::pre_resolution_tab( comodule_generic<algebroid,degree_type>& resolved, string filename_maps, string filename_generators, int resolution_length, std::vector<curtis_table<ring>*>& result, std::vector<std::vector<int>> &gens, matrix<ring> *inj, matrix<ring> *indj, matrix<ring> *qut, matrix<ring> *new_map, std::string tablename, std::vector<std::vector<int>> *basis_orders){
	//open the file for the maps
	std::fstream maps_file(filename_maps, std::fstream::out | std::fstream::binary);
	//open the file for the generators
	std::fstream gens_file(filename_generators, std::fstream::out | std::fstream::binary);
	//file for tables
	std::cout << tablename << "\n";
	std::fstream tablefile(tablename, std::fstream::out | std::fstream::binary);
	//check the file
	if(maps_file.is_open() && gens_file.is_open())
		std::cout << "\nfiles succesfully opened\n";
	else std::cerr << "\nfailed to open the files";
	if(!tablefile.is_open())
		std::cerr << "fail to open";
	
	//get the rank of the current comodule
	int32_t M_rank = resolved.rank();  gens_file.write((char*)&M_rank, 4);
	
	gens.resize(resolution_length+1);
	for(int i=0; i<=resolution_length; ++i){
		std::cout << "doing the " << i << "/" << resolution_length << "-th step\n" << std::flush;
		//the injection to a cofree one
		inj->clear();
		indj->clear();
		//the quotient map to the new comodule
		qut->clear();

		//take the resolution step
		std::vector<int> *basis_order;
		if(basis_orders == NULL) basis_order = NULL;
		else basis_order = &(*basis_orders)[i];
		auto F = resolvor(resolved, inj, indj, qut, result[i], new_map, &gens[i], basis_order, tablefile);
		
		std::cout << "saving data...\n" << std::flush;

		//save the injective map
		inj->save(maps_file);
		//save the quotient map
		qut->save(maps_file);
		//save the generators
		F.save(gens_file);
		
		//update the rank of current comodule
		M_rank = resolved.rank();  gens_file.write((char*)&M_rank, 4);
	}
}
