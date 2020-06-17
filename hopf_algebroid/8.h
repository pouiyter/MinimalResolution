//combing pre-resolutions into resolutions
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::resolution(string filename_maps, string filename_generators, string filename_resolution, int resolution_length, comodule_generic<algebroid,degree_type>& tmp, matrix<ring> *inj, matrix<ring>* qut, matrix<ring>* composed, std::ostream *cout){
	//open the files
	std::fstream maps_file(filename_maps, std::fstream::in | std::fstream::binary);
	std::fstream gens_file(filename_generators, std::fstream::in | std::fstream::binary);
	std::fstream res_file(filename_resolution, std::fstream::out | std::fstream::binary);
	if(maps_file.is_open() && gens_file.is_open() && res_file.is_open())
		std::cout << "\nfiles succesfully opened\n";
	else std::cerr << "\nfailed to open files\n";
	
	//read the rank of the first comodule
	int32_t M_rank;
	gens_file.read((char*)&M_rank, 4);
	gens_file.close();
	maps_file.close();
	
	//read the generators for the first cofree comodule
	cofree_comodule<algebroid,degree_type> F1;
	gens_file.open(filename_generators + "0", std::ios::in | std::ios::binary);
	maps_file.open(filename_maps + "0", std::ios::in | std::ios::binary);
	F1.load(gens_file);
	
	//load the first injective map
	inj->clear();
	qut->clear();
	inj->load(maps_file);
	
	for(int i=0; i<resolution_length; ++i){
		//the rank of the current comodule
		gens_file.read((char*)&M_rank, 4);
		gens_file.close();
		//load the next generators
		gens_file.open(filename_generators + std::to_string(i+1), std::ios::in | std::ios::binary);
		cofree_comodule<algebroid,degree_type> F2;
		F2.load(gens_file);
		//load the quotient map
		qut->load(maps_file);
		if(cout!=NULL)
			*cout << "qut:" << qut->output();
		maps_file.close();
		
		//load the injective map
		maps_file.open(filename_maps + std::to_string(i+1), std::ios::in | std::ios::binary);
		//the rule to modify the entries
		std::function<vectors<matrix_index,ring>(vectors<matrix_index,ring> const&)> rule = [&F2,this](vectors<matrix_index,ring> const &v){
			//replace the index using the free module data
			std::function<matrix_index(matrix_index)> ri = [&F2](matrix_index n){ 
				return F2.find_index(n); };
				matrix_index invalid = cofree_comodule<algebroid,degree_type>::invalid_pos;
				return moduleOper->filtered_reindex(ri, v, invalid);
		};
		
		inj->load_modify(maps_file, rule);
		if(cout!=NULL)
			*cout << "inj:" << inj->output() << std::flush;

		//compose the injective map with the quotient map
		composed->clear();
		inj->compose(qut, composed);
		if(cout!=NULL)
			*cout << "com:" << composed->output() << std::flush;
		composed->save(res_file);
		
		//go to the next step
		F1=F2;
	}
}

//the type of complexes
template<typename degree_type,typename ring>
class complex{
public:
	std::vector<modules<degree_type>> terms;
	std::vector<matrix<ring>*> maps;
	
	int size() const{ return terms.size(); }
};

