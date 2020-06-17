//taubocmain.cpp
#include"tao_bockstein.h"

//the main function to construct the tau-Bockstein
int main(int agrc, char **argv){
	//get the maximal degree and the length of the resolution
	int max_deg = std::atoi(argv[1]);
	int resolution_length = std::atoi(argv[2]);
	string directory(argv[1]);
	directory += "_";
	string director = directory;
	
	//initialize matrix class
	matrix<tauPoly>::moduleOper = &tau_module_oper;
	
	//the operator for motivic dual steenrid algebra
	MotSteenrodOp MOP(NULL, max_deg);
	//initialize the list of monomials
	MOP.init_mon_array(directory + "ex2poly_index");
	std::cout << MOP.output_monomials();
	
	//output the set of generators
	output_generators(resolution_length, director + "mot_gens", director + "mot_gen.txt");
	
	//the complex of the primitives
	motComplex Complex;
	
	//load the maps
	Complex.load(resolution_length, director + "mot_gens", director + "mot_res");
	
	//output the complex
	Complex.output(director + "mot_maps.txt");
	
	//computet the tau-Bockstein
	auto tables = make_table(Complex);
	
	//output the tables
	std::fstream taufile(director + "tau_bockstein.txt", std::ios::out);
	taufile << output_tables(tables);
	
	return 0;
}