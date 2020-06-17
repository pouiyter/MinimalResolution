//BP_init.h
#pragma once
#include"BP.h"
#include"matrices_mem.h"
#include"matrices_stream.h"
#include"algNov.h"
#include"multiplication.h"
#include"Boc.h"

//initialize a generic comodule
class BPComodInit : public BPCoMod_generic{
	matrix_file<BPBP> coaction_matrix;
public:
	//the constructor
	BPComodInit(string);
};

//the initialization of the data for BP
class BPInit{
public:
	//the maximal degee
	int max_degree;
	
	//the length of the resolution
	int resolution_length;
	
	//the director for the data
	string director;
	
	//the operations on Z2
	Z2_Op Z2_oper;
	//the operations on BP
	BP_Op BP_oper;
	//operations on F2-modules
	ModuleOp<matrix_index,F2> F2Mod_opers;
	
	//the matrices for the structure data
	matrix_file<BP> etaL_matrix, R2L_matrix;
	matrix_file<BPBP> delta_matrix;
	
	//the complex of primitives
	BPComplex Complex;
	//algebraic Novikov tables
	std::vector<algNov_table> AAN_table;
	algNov_tables AANtables;
	//Bockstein tables
	std::vector<Boc_table> B_table;
	Boc_tables Btables;
	
	//the container for the matrix of the injection to a cofree comodule, and the quotient to the next comodule
	matrix_mem<BP> inj;
	matrix_file<BP> indj, qut, new_map;
	//the matrices for the resolutions
	std::vector<matrix_mem<Z2>> mapses;
	matrix_file<Z2> mm;
	
	//the curtis table for the resolutions
	std::vector<curtis_table_mem<F2>> ResolutionTables;
	curtis_table_mem<F2> ctable;
	
	//the generators for a resolution
	std::vector<std::vector<int>> gens;
	
	//a comdule which is initialization to the trivial comodule of rank 1
	BPComodInit comod;
	
	//the operations for multiplicative structures
	multiplication multp;
	
	//load the curtis table data
	void loadResolutionTables(string table_data);
	
	//load the data for generators
	void load_gens(string gens_data);
	
	//the constructor
	BPInit(int max_deg, int resolution_length, string etaL_data, string delta_data, string R2L_data, string dirname);
	
	//do resolutions
	void resolve();
	
	//construct resolution
	void resolution();
	
	//make algebraic Novikov table
	void make_algNov();
	//make Bockstein table
	void make_Boc();
	
	//make multiplication table
	void mult_table(BPBP const &, int, string);
	void mult_table();
	//make multiplication table for top theta on the Moore spectrum
	void mult_theta();
};
