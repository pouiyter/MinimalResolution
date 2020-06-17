//ex_steenrod_init.h
#pragma once
#include"ctau_steenrod.h"
#include"matrices_mem.h"
#include"matrices_stream.h"
//#include"curtis.h"
#include<memory>

//initialize a generic comodule
class ComodInit : public SteenrodCoMod_generic{
	matrix_file<P> coaction_matrix;
public:
	//the constructor
	ComodInit(string);
};

//the initialization of the data for the steenrod algebra
class SteenrodInit{
public:
	//the length of the resolution
	int resolution_length;
	
	//the operations on the dual steenrod algebra
	Steenrod_Op steenrod_oper;
	
	//the matrix for the co-multiplication
	matrix_mem<P> deltaTable;
	
	//the container for the matrix of the injection to a cofree comodule, and the quotient to the next comodule
	matrix_file<Fp> inj; 
	matrix_file<Fp> indj, qut, new_map;
	
	//the curtis table for the resolutions
	std::vector<curtis_table<Fp>*> resolutionTables;
	std::vector<curtis_table_mem<Fp>> ResolutionTables;
//	std::vector<std::unique_ptr<curtisTable_stream<Fp>>> ResTables;
//	std::vector<std::unique_ptr<std::fstream>> TableDatas;
	
	//the generators for a resolution
	std::vector<std::vector<int>> gens;
	
	//a comdule which is initialization to the trivial comodule of rank 1
	ComodInit comod;
	
	//the constructor
	SteenrodInit(int prime, int max_deg, int resolution_length, string delta_data, string director, bool IorO=false);
	
	//do resolutions
	void resolve(string director, std::vector<std::vector<int>> *basis_orders = NULL);
	
	//save the resolution tables
	void saveResolutionTables(string);
	
	//save the generators
	void save_gens(string);
};
