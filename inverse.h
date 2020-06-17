//inverse.h
#pragma once
#include"matrices.h"

//find the partial inverse of an injective matrix using a table
template<typename ring>
void inverse(int source_rank, curtis_table<ring> &table, matrix<ring> &result, std::vector<int> &basis_order){
	std::cout << "rank = " << source_rank << std::flush;
	result.clear();

	for(int j=0; j<source_rank; ++j){
	//	int i = basis_order[j];
		int i = source_rank - j - 1;
		//compute the i-th row of the inverse matrix;
		//those rows not in the table goes to zero
		if(!table.is_member(i)){
			result.insert(i, table.ModOper->zero());
		//	std::cout << "\n" << i << std::flush;
			continue;
		}
		//otherwise use the table, such that the cycle of the i-th entry maps to the tag.
		auto et = table.search(i);
	//	std::cout <<"\n" << i << ";" << table.ModOper->output(et.full_cycle) << std::flush;
		//find the other terms in the cycle
		auto others = table.ModOper->minus(et.full_cycle);
		others = table.ModOper->add(table.ModOper->singleton(i), others);
		//find the image of the other terms
		auto oim = result.maps_to(others);
		//subtract this to the full_tag
		auto row = table.ModOper->add(et.full_tag, oim);
		result.insert(i, row);
	}
}

//compute the inverse maps
template<typename ring>
void inverse(std::vector<curtis_table<ring>*> tables, std::vector<int> res_rank, string filename_splitting, matrix<ring> *sp_mat, std::vector<std::vector<int>> &basis_orders){
	for(unsigned i=0; i<res_rank.size(); ++i){
		std::fstream split_file(filename_splitting + std::to_string(i), std::ios::out | std::ios::binary);
		inverse(res_rank[i], *(tables[i]), *sp_mat, basis_orders[i]);
		sp_mat->save(split_file);
	}
}

#include"hopf_algebroid.h"
//combing pre-resolutions into resolutions, return the vector of the rank of the resolutions
template<typename ring, typename algebroid, typename degree_type>
std::vector<int> resolution(Hopf_Algebroid<ring, algebroid> &HA_oper, string filename_maps, string filename_generators, string filename_resolution, int resolution_length, matrix<ring> *inj, matrix<ring>* qut, matrix<ring>* composed, std::ostream *cout, cofree_comodule<algebroid,degree_type>&){
	std::vector<int> res_rank;
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
	
	//read the generators for the first cofree comodule
	cofree_comodule<algebroid,degree_type> F1;
	F1.load(gens_file);
	res_rank.push_back(F1.rank());
	
	//load the first injective map
	inj->clear();
	qut->clear();
	inj->load(maps_file);
	
	for(int i=0; i<resolution_length; ++i){
		//the rank of the current comodule
		gens_file.read((char*)&M_rank, 4);
		
		//load the next generators
		cofree_comodule<algebroid,degree_type> F2;
		F2.load(gens_file);
		res_rank.push_back(F2.rank());
		//load the quotient map
		qut->load(maps_file);
		if(cout!=NULL)
			*cout << "qut:" << qut->output();
		
		//load the injective map
		inj->load(maps_file);
		if(cout!=NULL)
			*cout << "inj:" << inj->output() << std::flush;
		
		//compose the injective map with the quotient map
		composed->clear();
	//	inj->compose(qut, composed);
		if(cout!=NULL)
			*cout << "com:" << composed->output() << std::flush;
		composed->save(res_file);
		
		//go to the next step
		F1=F2; 
		
		std::cout << "\r"  << i << "/" << resolution_length << std::flush;
	}
	
	return res_rank;
}

//check the filename_splitting
template<typename ring>
bool check_splitting(int resolution_length,matrix<ring> *inj, matrix<ring> *split, matrix<ring> *qut, string filename_maps, string filename_splitting){
	std::fstream maps_file(filename_maps, std::ios::in | std::ios::binary);

	bool res = true;
	for(int i=0; i<resolution_length; ++i){
		std::fstream split_file(filename_splitting + std::to_string(i), std::ios::in | std::ios::binary);
		split->load(split_file);
		inj->load(maps_file);
		qut->load(maps_file);
		split->compose(inj, qut);
		
		inj->set2unit(qut->rank);
		res = res && inj->equal(*qut);
	}
	return res;
}
