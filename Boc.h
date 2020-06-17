//Boc.h
#pragma once
#include"SS.h"
#include"algNov.h"
#include"multiplication.h"

//the class for the Bockstein table
class Boc_table : public algNov_table{
public:
	//constructor
	Boc_table();
	
	//re-define the valuation so that it is given by the number of v0's
	int v_valuation(Z2, exponent) override;
	
	//number of v's
	int num_v(cycle_name);
	
	//output entries
	std::pair<std::pair<int,int>,string> output(SS_entry<cycle_name, Z2>&, int);
	
	//output table 
	std::set<std::pair<std::pair<int,int>,string>> output_table(int, int) override;
	
	//transform Bockstein name into Algebraic Novikov name
	multiplication_table<cycle_name> Bname2Aname(algNov_table &Atable, int pric, bool simple=true);
};

class Boc_tables : public algNov_tables{
	//this is a warning to prevent misuse
	void set_table(std::vector<algNov_table>*);
	std::vector<Boc_table*> Btables;
public:
	//initialize the table
	void set_table(std::vector<Boc_table>*);
	
	//transform Bockstein names into Algebraic Novikov names
	std::vector<multiplication_table<cycle_name>> Bname2Anames(int resolution_length, algNov_tables&, int pric);
};
