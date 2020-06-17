//BPcomplex.h
#pragma once
#include"BP.h"

//the data for an Zp generator of primitives, of the form v^coeficient[gen_pos]
class prim_entry{
public:
	//the position of the generator
	matrix_index gen_pos;
	//the exponent of the coeficient
	exponent coeficient;
	
	//the comparison
	bool operator<(const prim_entry) const;
	
	//output of entries
	string output() const;
};

//the data for the primitives, i.e. the generators of the primitives over Zp
class primitive_data : public std::vector<prim_entry>{
public:
	//BP operations
	static BP_Op *BPoper;
	//Zp module operations
	static ModuleOp<matrix_index,Z2> *Z2Mod_oper;
	//index for searching
	std::map<prim_entry,unsigned> prim_index;
	//the shift of the generators, i.e. the place of the primitives in the whole complex
	std::vector<matrix_index> gen_shift;
	//degree of generators
	std::vector<int> gen_deg;
	
	//set the operators, remember to run this before doing anything
	static void set_oper(BP_Op*);
	
	//add a new entry
	void add(prim_entry,int,int);
	
	//make the primitive data using the list of monimials
	static primitive_data make_primitives(int degree, int pos, int, monomial_index*);
	
	//combining two datas
	void direct_summed(const primitive_data&);
	
	//make the data of primitives from a cofree comodule
	static primitive_data primitive_maker(cofree_comodule<BPBP,int> &, monomial_index*);

	vectors<matrix_index, Z2> expand(const vectors<matrix_index,BP>&) const;
	
	//output primitive data
	string output();
};

//the class for complex of BP_* modules
class BPComplex{
public:
	std::vector<primitive_data> Prims;
	std::vector<matrix<Z2>*> Maps;
	
	//construct the complex of primitives from a resolution
	void load(int resolution_length, string generator_filename, string maps_filename, matrix<BP>*, matrix<BP>*, std::function<matrix<Z2>*(int)>);
	
	//IO operations
	void save_matrix(string);
	void load_matrix(int resolution_length, string generator_filename, string matrix_filename);
	std::vector<int> load_prim(int resolution_length, std::fstream &gens_file);
        
	//the size of the compex
	unsigned size();
	
	//construct the map of primitives, using the data source ->^{X} target
	static void make_prim_map(const primitive_data& source, const primitive_data& target, matrix<BP>* X, matrix<BP> *, matrix<Z2> *);
	
	//get the data for the generators
	static std::vector<FreeBPCoMod> get_generator(int resolution_length, string generator_filename);
};
