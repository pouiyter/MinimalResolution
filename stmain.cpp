//stmain.cpp
#include"steenrod_init.h"

int main(int argc, char** argv){
	string filename = argv[1];
	filename += "_";
	
	int maxdeg = std::atoi(argv[1]);
	int length = std::atoi(argv[2]);
	
	SteenrodInit st(2,maxdeg, length, filename + "steenrod_coaction.data");
	
//	std::cout<<st.steenrod_oper.delta_table->output();
	
	st.resolve(filename);
	st.saveResolutionTables(filename + "ResTables");
	st.save_gens(filename + "gens_data");
	
	std::cout << st.steenrod_oper.output_resolution(filename+"gens",length);
	
	return 0;
}
