//lift.h
#pragma once
#include"hopf_algebroid.h"

//compute the map to the co-generators F->cogen(F)
template<typename ring, typename algebroid, typename  degree_type>
vectors<matrix_index,ring> cogens_map(vectors<matrix_index,ring> &v, cofree_comodule<algebroid, degree_type> &F, ModuleOp<matrix_index, ring> &Module_oper){
	std::function<matrix_index(matrix_index)> rule = [](matrix_index i){
		matrix_index h = F.find_index(i);
		if(h==cofree_comodule<algebroid, degree_type>::invalid_pos) 
			return ModuleOp<matrix_index, ring>::invalid;
		return h;
	};
	return Module_oper->filtered_reindex(rule, v);
}

//compute the adjoind map of a map to cogenerators, return the i-th row
template<typename ring, typename algebroid, typename  degree_type>
vectors<matrix_index,ring> adjoint(matrix<ring> *X, std::vector<uint32_t> const& pos, int i, Hopf_Algebroid<ring, algebroid> &HA_oper){
	// psi(e_i)
	auto ci = F.coaction(i);
	// (1\otimes X) (psi (e_i))
	std::function<vectors<matrix_index,algebroid>(const algebroid&, const vectors<matrix_index,ring>&)> right_mult = [&HA_oper] (const algebroid& A, const vectors<matrix_index,ring> &V){
		return HA_oper.right_scalor_mult(A,V); };
	auto av = X->maps_to(ci,right_mult, HA_oper.algebroidModuleOper); 
	//change to vector notation
	vectors<matrix_index, ring> result;
	for(int k=0; k<F.pos.size(); ++k){
		algebroid kp = HA_oper->algebroidModuleOper->component(i);
		auto nd = HA_oper->algebroid2vector(kp, pos[i]);
		result.direct_sum(nd, 0);
	}
	return result;
}

//lift a map to its one-step resolution
template<typename ring, typename algebroid, typename degree_type>
void lift_resolvor(cofree_comodule<algebroid, degree_type> &F2, matrix<ring> &inj2, matrix<ring> &qut2, matrix<ring> &spliter1, std::vector<matrix_index> &inv_ind1, matrix<ring> &starting_map, matrix<ring> &result, Hopf_Algebroid<ring, algebroid> &HA_oper, matrix<ring> *lg, matrix<ring> *cofree_map){
	result.clear();
	//compute the i-th term of the lifting F1->cogens(F2)
	std::function<vectors<matrix_index, ring>(int i)> cogens_lifting = [&spliter1, &starting_map, &inj2, &HA_oper, &F2](int i){
		//get the i-th row of the spliting map 
		auto spv = spliter1.find(i);
		//the image in X2
		auto pva = starting_map.maps_to(spv);
		//find the image in F2
		auto pb = inj2.maps_to(pva);
		//find the map to the cogenerators
		return cogens_map(pb, F2, HA_oper);
	}
	lg->construct(spliter1.rank, cogens_lifting);
	
	if(cofree_map!=NULL){
		std::function<vectors<matrix_index, ring>(int i)> F1F2 = []{
			return adjoint(lg, F2.position_of_gens, i, HA_oper); };
		cofree_map->construct(lg->rank, F1F2);
	}
	
	//compute the i-th term of the map F1->M2
	std::function<vectors<matrix_index, ring>(int i)> F1M2 = [&cogens_map, &qut2, &inv_ind1](int i){
		auto F1F2 = adjoint(lg, F2.position_of_gens, inv_ind1[i], HA_oper);
		return qut2.maps_to(F1F2); 
	};
	
	//construct the matrix
	result.construct(inv_ind1.size(), F1M2);
}

//lift a map between resolutions
template<typename ring, typename algebroid>
void resolution_lift(int resolution_length, matrix<ring> *current_map, string splitfile1, string gensfile2, string mapsfile2, std::vector<std::vector<matrix_index>> &inv_ind1, matrix<ring> *next_map, Hopf_Algebroid<ring,algebroid> &HA_oper, string liftedfile, matrix<ring> *inj2, matrix<ring> *qut2, matrix<ring> *spliter1, matrix<ring> *liftedmap, std::iostream *cofree_map_file = NULL, matrix<ring> *cofree_map){
	//generators for the resolution of the second comodule
	std::fstream gens_file2(gensfile2, std::ios::in | std::ios::binary);
	//maps in the resolution of the second comodule
	std::fstream maps_file2(mapsfile2, std::ios::in | std::ios::binary);
	for(int i=0; i<resolution_length; ++i){
		//file for the splitting of the resolution of first comodule
		std::fstream split_file1(splitfile1 + std::to_string(i), std::ios::in | std::ios::binary);
		//file for the lifted maps
		std::fstream lift_file(liftedfile + std::to_string(i), std::ios::out | std::ios::binary);
		
		//read the rank of second resolution
		int32_t M_rank;
		gens_file2.read((char*)&M_rank, 4);
		
		//read the generators for second resolution
		cofree_comodule<algebroid,degree_type> F2;
		F2.load(gens_file);
		
		//load the maps
		inj2->clear();
		qut2->clear();
		inj2->load(maps_file);
		qut2->load(maps_file);
		
		//load the spliting map for the first resolution
		spliter1->load(split_file1);
		
		//lift the map
		lift_resolvor(F2, *inj2, *qut2, *spliter1, inv_ind1[i], *current_map, *next_map, HA_oper, liftedmap, cofree_map);
		
		//save the maps
		liftedmap->save(lift_file);
		if(cofree_map_file!=NULL)
			cofree_map->save(*cofree_map_file);
		
		//renew the comodule maps
		current_map->construct(*next_map);
	}
}

//make the quotient index for an injective map X->F, with data from a table
template<typename ring>
std::vector<matrix_index> make_quot_index(curtis_table<ring> *table, int Frank, int Xrank){
	auto inj_ind = table->cycle_matrix(Xrank());
	auto quot_inds = matrix<ring>::quot_index(inj_ind, Frank);
}
