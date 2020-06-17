//BPmain.cpp
#include"BP.h"
#include"BP_init.h"
#include<cstring>

int main(int argc, char** argv){
	//set the maximal degree
	int max_degree = std::atoi(argv[1]);
	int resolution_length = std::atoi(argv[2]);
	//initialize the monomial index
	monomial_index mon_index(max_degree);
	
	string filename = argv[1];
	string filename0 = filename + "_";
	filename += "_BP";
	
	//construct the operators for BP
	BPInit BPoper(max_degree, resolution_length, filename0 + "etaL", filename0 + "R2L", filename0 + "delta", filename);
	
	//load the resolution table for BP/I
//	BPoper.loadResolutionTables(filename0 + "ResTables");
	//load the generators
	BPoper.load_gens(filename0 + "gens_data");
	
	std::cout << "starting resolution..." << std::flush;
	//do the resolution
	BPoper.resolve();
	
	//construct the resolution
	BPoper.resolution();
	
	//compute algebraic Novikov
	BPoper.make_algNov();
	
	//compute Bocstein table
	BPoper.make_Boc();
	
	//compute the multiplicative structure
	//BPoper.mult_table();
	BPoper.mult_theta();
	
	return 0;
}
