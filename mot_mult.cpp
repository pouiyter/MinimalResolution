//mot_mult.cpp
#include "hopf_algebroid.h"
#include "mot_steenrod.h"
#include"tao_bockstein.h"
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
	//the complex of the primitives
	motComplex Complex;
	//load the maps
	Complex.load(resolution_length, director + "mot_gens", director + "mot_res");
	
	//computet the tau-Bockstein
	auto tables = make_table(Complex);
	//output the tables
	std::fstream taufile(director + "tau_bockstein.txt", std::ios::out);
	taufile << output_tables(tables);
	
	//load the new generators with diagonize the boundary maps
	std::ofstream ngf(director + "new_generators");
	std::vector<cycle_data> cycles_table(tables.size());
	make_cycle_tables(tables,cycles_table);
	for(int i=0;i<tables.size();++i){
		ngf << output(cycles_table[i], i);
	}
	
	//compute the h0 multiplication
	multiplication_table(MOP.hi(0), 1, director + "mot_gens", director + "mot_res", director + "h0.txt", resolution_length, MOP, cycles_table); 
	multiplication_table(MOP.hi(1), 2, director + "mot_gens", director + "mot_res", director + "h1.txt", resolution_length, MOP, cycles_table); 
	multiplication_table(MOP.hi(2), 4, director + "mot_gens", director + "mot_res", director + "h2.txt", resolution_length, MOP, cycles_table); 
	multiplication_table(MOP.hi(3), 8, director + "mot_gens", director + "mot_res", director + "h3.txt", resolution_length, MOP, cycles_table); 
}