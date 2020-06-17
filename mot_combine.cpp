//mot_combine.cpp
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

	//the operator for motivic dual steenrid algebra
	MotSteenrodOp MOP(NULL, max_deg);
	//initialize the list of monomials
	MOP.init_mon_array(directory + "ex2poly_index");
	std::cout << MOP.output_monomials();
	
	//the comodule to be resolved
	matrix_file<motSteenrod> coaction_matrix(director + "mot_coac_d");
	comodule_generic<motSteenrod ,MotDegree> comod(&coaction_matrix);
	//set comod to be the homology of the sphere
	MotDegree zo = {0,0};
	MOP.set_to_trivial(comod,zo);
		
	//the matrices for the intermediat steps
	matrix_file<tauPoly> inj(director + "motinj"), qut(director + "motqut"), indj(director + "motindj");
		
	//combining generator files
	MOP.gens_file_combiner(director + "mot_gens", resolution_length, comod);
		
	//construct the resolution
	std::fstream outfile("maps.txt", std::ios::out);
	MOP.resolution(director + "mot_maps", director + "mot_gens", director + "mot_res" , resolution_length, comod, &inj, &qut, &indj, &outfile);
}