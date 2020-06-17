//mot_main.cpp
#include "hopf_algebroid.h"
#include "mot_steenrod.h"
#include "matrices_stream.h"

int main(int agrc, char **argv){
	//get the maximal degree and the length of the resolution
	int max_deg = std::atoi(argv[1]);
	int resolution_length = std::atoi(argv[2]);
	string directory(argv[1]);
	directory += "_";
	string director = directory;
  
	//initialize matrix class
	matrix<tauPoly>::moduleOper = &tau_module_oper;
	matrix<motSteenrod>::moduleOper = &motSteenrod_module_oper;
	//the matrix of the coactions
	matrix_mem<motSteenrod> coa;
	//the operator for motivic dual steenrid algebra
	MotSteenrodOp MOP(&coa, max_deg);
	//initialize the list of monomials
	MOP.init_mon_array(directory + "ex2poly_index");
	std::cout << MOP.output_monomials();

	MOP.generate_cofree_coaction(directory + "mot_deltas", directory + "poly_exponents");
	// std::cout << coa.output();
	
	//the comodule to be resolved
	matrix_file<motSteenrod> coaction_matrix(director + "mot_coac_d");
	comodule_generic<motSteenrod ,MotDegree> comod(&coaction_matrix);
	//set comod to be the homology of the sphere
	MotDegree zo = {0,0};
	MOP.set_to_trivial(comod,zo);
	
	//operations over F2
	Fp_Op F2oper(2);
	ModuleOp<matrix_index,F2> F2ModuleOper(&F2oper);
	
	//set the curtis table
	curtis_table<F2>::ModOper = &F2ModuleOper;
	curtis_table_mem<F2> ctable;
	
	//transform a F2 vector to a tau vector
	std::function<vectors<matrix_index,tauPoly>(const vectors<matrix_index,F2>&)> tfm = [&MOP] (const vectors<matrix_index,F2>& v){ 
		return MOP.lift(v); };
		
	//the matrices for the intermediat steps
	matrix_file<tauPoly> inj(director + "motinj"), qut(director + "motqut"), indj(director + "motindj"), new_map(director + "motnewm");
	
	//the generators for a resolution
	std::vector<std::vector<int>> gens;
	
	//load the generators
	MOP.load_gens(gens, director + "gens_data_ctau");
	
	//construct the pre-resolution
	MOP.pre_resolution_modeled(comod, director + "mot_maps", director + "mot_gens", resolution_length, director+"extables", &ctable, gens, tfm, &inj, &qut, &indj, &new_map, director + "back");
	
	//combining generator files
	MOP.gens_file_combiner(director + "mot_gens", resolution_length, comod);

	//construct the resolution
	std::fstream outfile("maps.txt", std::ios::out);
	MOP.resolution(director + "mot_maps", director + "mot_gens", director + "mot_res" , resolution_length, comod, &inj, &qut, &indj, &outfile);
}