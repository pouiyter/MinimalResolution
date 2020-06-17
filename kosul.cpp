//kosul.cpp
#include"hopf_algebroid.h"
#include"steenrod.h"
#include"steenrod_init.h"
#include"inverse.h"

//construct the double of a comodule M
void double_comod(CoModule<P,int> &M, SteenrodCoMod_generic &result, Steenrod_Op *Stoper, int max_deg){
	//double the degree
	std::vector<matrix_index> gens;
	for(int i=0; i<M.rank(); ++i){
		if(M.degree(i) * 2 > max_deg)
			continue;
		result.base_module.degree.push_back(M.degree(i) * 2);
		gens.push_back(i);
	}
	result.base_module.rank = gens.size();
	//double the coaction
	result.coaction_matrix->clear();
	std::function<P(const P&)> square = [Stoper](const P& x){
		return Stoper->P_opers.multiply(x, x); };
	std::function<vectors<matrix_index, P>(int)> ca = [Stoper, square, &M, &gens](int i){
		return Stoper->PMod_opers.termwise_operation(square, M.coaction(gens[i])); };
	result.coaction_matrix->construct(result.rank(), ca);
}

//construct the double of a cofree comodule
void cofree_double(SteenrodCoMod_generic &result, Steenrod_Op *Stoper){
	// construct a cofree comodule of rank 1
	FreeSteenrodCoMod cf1;
	cf1.generators.rank = 1;
	cf1.generators.degree.push_back(0);
	cf1.position_of_gens.push_back(0);
	cf1.total_rank = Stoper->mon_index.number_of_all_mons();
	//double
	double_comod(cf1, result, Stoper, Stoper->mon_index.max_degree);
}

//main function to construct the kozsul resolution
int main(int argc, char** argv){
	string filename = argv[1];
	filename += "_";
	
	int maxdeg = std::atoi(argv[1]);
	int length = std::atoi(argv[2]);
	
	SteenrodInit st(2,maxdeg, length, filename + "steenrod_coaction.data");
	
	//set the double of a cofree comodule
	cofree_double(st.comod, &st.steenrod_oper);
	
	//compute the resolution
	std::vector<std::vector<int>> basis_orders(length+2);
	st.resolve(filename, &basis_orders);
	st.saveResolutionTables(filename + "ResTables");
	st.save_gens(filename + "gens_data");
	
	FreeSteenrodCoMod cc;
	
	//compute the resolution
	std::vector<int> res_rank = resolution(st.steenrod_oper, filename + "maps", filename + "gens", filename + "res", length-1, &st.inj, &st.qut, &st.indj, NULL, cc);
	
	std::cout << "resolution complete\n" <<std::flush;
 	
	//compute the splitting
	inverse(st.resolutionTables, res_rank, filename + "splitting", &st.indj, basis_orders);
	
	//check the splitting
	if(check_splitting(length, &st.inj, &st.indj, &st.qut, filename + "maps", filename + "splitting"))
		std::cout << "success";
	else std::cout << "fail";
	
	return 0;
}
