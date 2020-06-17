//steenrod_init.cpp
#include"steenrod_init.h"

//the initialization
SteenrodInit::SteenrodInit(int prime, int max_deg, int res_length, string delta_data, bool IorO) : steenrod_oper(max_deg, prime, &deltaTable){
	//the delta data file does not exist
	if(!IorO){
		//open the file for writing
		std::fstream delta_file(delta_data, std::ios::out | std::ios::binary | std::ios::trunc);
		
		//construct the delta data
		std::cout << "constructing the co-product for the dual Steenrod algebra...\n" << std::flush;
		steenrod_oper.make_delta(delta_file);
		delta_file.close();
		
		resolution_length = res_length;
		
		std::cout << "maximal degree:" << steenrod_oper.maxDeg << "\n";
		std::cout << "resolution length:" << res_length << "\n";
	}
	
	//open the data for delta
	std::fstream delta_file(delta_data, std::ios::in | std::ios::binary);
	steenrod_oper.initialize(delta_file);
	delta_file.close();
	
	//initialize matrix class
	matrix<Fp>::moduleOper = &steenrod_oper.FpMod_opers;
	matrix<P>::moduleOper = &steenrod_oper.PMod_opers;
	
	//initialize the comod to a trivial one with one generator at degree 0
	steenrod_oper.set_to_trivial(comod,0);
	
	//initialize the curtis tables
	curtis_table<Fp>::ModOper = &steenrod_oper.FpMod_opers;
	ResolutionTables.resize(resolution_length+2);
	for(int i=0; i<resolution_length+2; ++i)
		resolutionTables.push_back(&ResolutionTables[i]);
}

//initialize a generic comodule
ComodInit::ComodInit() : SteenrodCoMod_generic(&coaction_matrix){
}

//do resolutions
void SteenrodInit::resolve(string director, std::vector<std::vector<int>> *basis_orders){
	steenrod_oper.pre_resolution_tab(comod, director + "maps", director + "gens", resolution_length, resolutionTables, gens, &inj, &indj, &qut, &new_map, director + "BPtables");
}


//save the resolution tables
void SteenrodInit::saveResolutionTables(string table_data){
	std::fstream tables(table_data, std::ios::out | std::ios::binary);
	for(unsigned i=0; i<ResolutionTables.size(); ++i){
		ResolutionTables[i].save(tables);
	}
}

//save the generators, 4 bytes size, then for each of the generator sets
void SteenrodInit::save_gens(string gens_data){
	std::fstream genfile(gens_data, std::ios::out | std::ios::binary);
	int32_t sz = gens.size();
	genfile.write((char*)&sz, 4);
	for(int i=0; i<sz; ++i){
		int32_t ss = gens[i].size();
		genfile.write((char*)&ss, 4);
		for(int j=0; j<ss; ++j)
			genfile.write((char*)&gens[i][j], 4);
	}
}
