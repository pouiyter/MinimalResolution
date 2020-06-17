//stmain.cpp
#include"ex_steenrod_init.h"

int main(int argc, char** argv){
	string filename = argv[1];
	filename += "_";
	
	int maxdeg = std::atoi(argv[1]);
	int length = std::atoi(argv[2]);
	
	SteenrodInit st(2,maxdeg, length, filename + "ctau_steenrod_coaction.data", filename);
	
//	std::cout<<st.steenrod_oper.delta_table->output();
	
	st.resolve(filename);
	st.saveResolutionTables(filename + "ResTables_ctau");
	st.save_gens(filename + "gens_data_ctau");
	
	std::cout << st.steenrod_oper.output_resolution(filename+"gens",length);
	
	return 0;
}
