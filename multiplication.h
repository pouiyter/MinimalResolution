//multiplication.h
#pragma once

#include "algNov.h"

//the entries in a multiplacation table
template<typename cycle_name>
struct multiplication_table_entry {
	//the name of the original element
	cycle_name original_name;
	//the name of the result by multiplacation
	std::pair<std::vector<cycle_name>, std::vector<cycle_name>> multiplied_names;
};
//the multiplacation table
template<typename cycle_name>
using  multiplication_table = std::vector<multiplication_table_entry<cycle_name>>;


//the class to produce the multiplicative structure
class multiplication{
	//operators for BP
	BP_Op *BPoper;
	
public:
	//constructor
	multiplication(BP_Op*);
	
	//compute the multiplication table for x->v\eta_R(x)
	void make_eta_R_multiplier(BPBP const&,matrix<BP>* result, int deg);
	
	//construct the multiplication matrix
	void multily_matrix(int max_deg, matrix<BP> *multiplier, matrix<BP> *resolv_map, FreeBPCoMod*, primitive_data*, matrix<Z2> *result);
	
	template<typename cycle_name, typename ring>
	static multiplication_table<cycle_name> mult_extension(matrix<ring>*, SS_table<cycle_name,ring>&, SS_table<cycle_name,ring>&, SS_table<cycle_name,ring>&, int);
	
	//make multiplacation table from the algebraic Novikov table of a complex
	std::vector<multiplication_table<cycle_name>> mult_extension(matrix<BP>*, int, int, std::vector<FreeBPCoMod>&, string, BPComplex&, algNov_tables, matrix<BP>*, matrix<Z2>*, int pric, bool fixedpric=false);
	
	//computing the table for multiplication by two
	template<typename cycle_name, typename ring>
	multiplication_table<cycle_name> two_extension(SS_table<cycle_name,ring>& cur_table, int pric);
	
	std::vector<multiplication_table<cycle_name>> two_extension(int resolution_length, BPComplex& Cm, algNov_tables Tb, int pric);
	
	//compute the multiplacation by two table in the whole algebraic Novikov table
	std::vector<multiplication_table<cycle_name>> two_extension(algNov_tables, BPComplex&, int pric);
	
	//output the table
	template<typename cycle_name, typename ring>
	string output_multiplication_table(multiplication_table<cycle_name> const&, int, int, int pric, SS_table<cycle_name, ring>&);
	
	//output the tables
	string output_multiplication_table(std::vector<multiplication_table<cycle_name>> const&, int, int pric);
};
