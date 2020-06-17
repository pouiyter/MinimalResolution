//algNov.h
#pragma once
#include"BP.h"
#include"SS.h"
#include"BPcomplex.h"

//the type of the name of cycles
//((filtration, v0, v1, v2, \dots), gen_ind)
typedef std::pair<std::vector<int>,int> cycle_name; // gen_ind, v0, v1, ...

//the class of algebraic Novikov spectral sequences
class algNov_table : virtual public SS_table<cycle_name, Z2>{
public:
	//basic operations
	static Z2_Op *Z2Oper;
	static ModuleOp<matrix_index, Z2> *Modop;
	static BP_Op *BPoper;
	
	//set the operations
	static void set_op(BP_Op*);
	
	//the space for the tags
	primitive_data *Ptag;
	//the space for the cycles
	primitive_data *Pcyc;

	//the filtration by number of v's
	virtual int v_valuation(Z2, exponent);
	
	//filtration of a cycle
	int filtration(cycle_name);
	
	//transform a term into a cycle name with given primitive data
	cycle_name naming(const Z2&, matrix_index, primitive_data&);
	
	//transform a term into a cycle name, using Pcyc
	cycle_name naming(const Z2&, matrix_index);
	
	//computing the leading term of a cycle vector
	unsigned leading_term(SS_entry<cycle_name, Z2>::value_type&);
	
	//constructors
	algNov_table();
	
	//output of cycle names
	string output(cycle_name,int);
	
	//output entries
	std::pair<std::pair<int,int>,string> output(SS_entry<cycle_name, Z2>&, int);
	
	//output the table
	virtual std::set<std::pair<std::pair<int,int>,string>> output_table(int, int);
	
	//IO of cycle names
 	void save(cycle_name, std::iostream&);
	cycle_name load(std::iostream&);
	
	//IO of entries
	void save(SS_entry<cycle_name, Z2>&, std::iostream&);
	void load(SS_entry<cycle_name, Z2>&, std::iostream&);
	
	//construct a vector using the primitive data
	static SS_entry<cycle_name, Z2>::value_type make_vec(cycle_name, primitive_data&, ModuleOp<matrix_index,Z2>*, Z2_Op*);
	//construct a full tag from the leading term
	SS_entry<cycle_name, Z2>::value_type get_tag(cycle_name);
	
	//construct the cycle pot. pot has additional leading number for sorting reasons...
	std::set<cycle_name> cycle_pot(int);
	
	//check the entry is a tag or not
	bool tagged(const SS_entry<cycle_name,Z2>&);
	
	//check if a name is valid
	bool valid(const cycle_name&);
	
	//the invalid cycle name
	cycle_name invalid();
	
	//the degree of an entry
	int degree(cycle_name);
};
   
    
class algNov_tables{
public:
	//the pointer to tables
	std::vector<algNov_table*> tables;
	//set the data of tables
	void set_table(std::vector<algNov_table>*);
	//construct the table from a complex of primitives
	void table_of_complex(BPComplex&,int pric);
	//output the table
	string output_tables();
	//IO operations
	void save(std::iostream&);
	void load(std::iostream&);
	void save(string);
	void load(string);
};
    
  
