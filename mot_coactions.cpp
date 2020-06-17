//mot_coactions.cpp
#include"mot_steenrod.h"
#include"mon_index.h"

//r t_n2^e2 \otimes t_n1^e1
mStmSt biMon(int n1, int e1, int n2, int e2, tauPoly r){
	//r t_n2^e2
	motSteenrod inner = motSteenrod_oper.monomial(singleVar(n2,e2),r);
	//inner x_n1^e1
	return mStmSt_oper.monomial(singleVar(n1,e1),inner);
}

//2^k
int power(int k){ return 1<<k; }

//construct the data for the comultiplications on the motivic dual Steenrod algebra
void make_delta_table(string deltafilename, monomial_index &mon_index){
	std::fstream deltafile(deltafilename, std::ios::out | std::ios::binary);
	//the co-multiplication of the generators
	std::function<mStmSt(int)> deltaXn = [](int n){
		// [x_n|1] + [1|x_n]
		auto res = mStmSt_oper.add(biMon(n,1,0,0,0), biMon(0,0,n,1,0));
		// tau^{-2^{i-1}} [x_i|x_{n-i}^{2^i}]. Note that \xi_i =  \tau^{-1} \tau_i^2
		for(int i=1;i<n;++i)
			res = mStmSt_oper.add(res, biMon(n-i,power(i),i,1,-power(i-1)));
		return res;
	};
	//construct the table
	mon_index.substitution_table(deltaXn, deltafile, &mStmSt_oper);
}

//output the coactions table
string output_coactions(string coaction_filename, string expoents_filename){
	string result;
	
	mStmSt_oper.output_term = [](exponent e, motSteenrod y){ 
		return "(" + motSteenrod_oper.output(y) + ")" + ::output(e); };
	
	std::fstream coaction_file(coaction_filename, std::ios::in | std::ios::binary);
	std::fstream expoents_file(expoents_filename, std::ios::in | std::ios::binary);
	
	int32_t mon_num;
	expoents_file.read((char*)&mon_num, 4);
	for(int i=0; i<mon_num; ++i){
		exponentArry eps = load_expArry(expoents_file);
		result += ::output(pack(eps));
		result += ":";
		mStmSt ca = mStmSt_oper.load(coaction_file);
		result += mStmSt_oper.output(ca);
		result += "\n";
	} 
	return result;
}

//constructing the data for the coaction and the list of monomials
int main(int argv, char** argc){
	int max_deg = std::atoi(argc[1]);
	monomial_index mon_index(max_deg);
	mon_index.init_mon_array();
	string filename = std::to_string(max_deg) + "_";
	//construct the 
	make_delta_table(filename + "mot_deltas", mon_index);
	//construct the list of exponents
	std::fstream monf(filename + "poly_exponents", std::ios::out | std::ios::binary);
	int32_t mon_num = mon_index.mon_array.size();
	monf.write((char*)&mon_num, 4);
	for(auto e : mon_index.mon_array){
		auto eps = unpack(e);
		save_expArry(eps,monf);
	}
	monf.close();
	
//	std::cout << output_coactions(filename + "mot_deltas", filename + "poly_exponents");
	return 0;
}