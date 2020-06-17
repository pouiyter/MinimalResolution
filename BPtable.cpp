//BPtable.cpp
#include"BP.h"
#include"matrices_mem.h"
#include<cstring>

int main(int argc, char** argv){
	//set the maximal degree
	int max_degree = std::atoi(argv[1]);
	//initialize the monomial index
	monomial_index mon_index(max_degree);
	//set the maximal generator
	int max_var = mon_index.max_var;
	
	string filename = argv[1];
	filename += "_";
	
	//construct rational operators
	Q2_int Q2intoper;
	Q2_Op Q2_oper;
	BPQ_Op BPQ_oper(&Q2_oper,&Q2intoper,max_var);
	
	//output the structure data
	std::fstream data_gens(filename + "structures", std::ios::out);
	data_gens << BPQ_oper.show_li();
	data_gens << BPQ_oper.show_vn();
	data_gens << BPQ_oper.show_etaR();
	data_gens.close();
	
	//construct the binary structure tables
	std::fstream R2Lfile(filename + "R2L_gen", std::ios::out | std::ios::binary);
	std::fstream deltafile(filename + "delta_gen", std::ios::out | std::ios::binary);
	BPQ_oper.output_R2L(R2Lfile);
	BPQ_oper.output_delta(deltafile);
	R2Lfile.close();
	deltafile.close();
	
	//construct the integral operators
	Z2_Op Z2_oper;
	matrix_mem<BP> etaL_matrix;
	BP_Op BP_oper(max_degree, &Z2_oper, &etaL_matrix, NULL, NULL);
	matrix<BP>::moduleOper = &BP_oper.BPMod_opers;
	
	//construct the tables for the complete structure maps
	std::fstream ouput(filename + "Ls", std::ios::out);
	BP_oper.make_tables(max_var, filename + "R2L_gen", filename + "delta_gen", filename + "etaL", filename + "R2L", filename + "delta", ouput);
}
