//tao_bockstein.h
#include "mot_steenrod.h"
#include"matrices_mem.h"
#include<stack>
#include<map>
#pragma once

//make the resolution 
void resolution(string filename_resolution, string filename_maps, string filename_generators, int resolution_length);
void resolution_truns(string filename_res_truns, string filename_maps, string filename_generators, int resolution_length);

//make the multiplication table
void multiplication_table(motSteenrod const&, int deg, string filename_generators, string filename_res_truns, string filename_outs, int resolution_length, MotSteenrodOp&); 
    
//output the set of generators
void output_generators(int resolution_length, string filename_generators, string filename_outs);

//complex over F2[tau]
class motComplex : public complex<MotDegree, tauPoly>{
	std::vector<matrix_mem<tauPoly>> maps_data;
public:
	//load the matrix from a resolution
	void load(int resolution_length, string filename_generators, string filename_resolution);
	//output to a file
	void output(string filename_outs);
};

class cycle_name{
	int i;
};

//the entries in a tau-Bockstein table
class tau_table_entry{
public:
	//the tag and cycle of an entry
	int tag, cycle;
	//the full data
	vectors<matrix_index, tauPoly> full_tag, full_cycle;
	//the length of the differential
	int diff_length;
		
	//the index for an unvalid taf
	static constexpr int Invalid = -1;
	//output the entry	
	string output(int) const;
};

//the data for the non-trivial cycles
typedef std::map<int,vectors<matrix_index,tauPoly>> cycle_data;

//tau-Bockstein table
class tau_table{
public:
	//the array of entries
	std::vector<tau_table_entry> table;
	//the index for the tags and cycles
	std::map<int,int> tag_index, cycle_index;
	//the leading term of a vector over F2[tau]
	static int leading_term(vectors<matrix_index, tauPoly> const&);
	//construct the bockstein table 
	static void make_pretable(std::stack<int> &pot, matrix<tauPoly> *M, tau_table &result, tau_table &pre_table);
	//insert a new entry
	void insert(int tag, int cycle, vectors<matrix_index, tauPoly> const &full_tag, vectors<matrix_index, tauPoly> const &full_cycle, int diff_length);
	//modify the table, the corresponding old entry in the n-th place is erased.
	void modify(int n, int tag, int cycle, vectors<matrix_index, tauPoly> const &full_tag, vectors<matrix_index, tauPoly> const &full_cycle, int diff_length);
	//make the pot for the next table
	std::stack<int> pot_maker(int rank) const;
	//output the table
 	string output(int nm) const;
	//get the data for the cycles
	cycle_data get_cycles() const;
	//get the data for the tags
	cycle_data get_tags() const;
};

//make the tables from the complex
std::vector<tau_table> make_table(motComplex&);
//output the tables
string output_tables(const std::vector<tau_table>&);

//get the cycle data from a table
void make_cycle_tables(std::vector<tau_table> const &tab, std::vector<cycle_data> &result);
//output the data of the cycles
string output(cycle_data const&, int);

//the mutliplication table
void multiplication_table(motSteenrod const&, int deg, string filename_generators, string filename_res_truns, string filename_outs, int resolution_length, MotSteenrodOp&, std::vector<cycle_data>&); 
