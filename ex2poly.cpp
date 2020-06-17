//ex2poly.cpp
#include"ex_exponents.h"
#include"ex_index.h"
#include<fstream>


//output the exponent
string output_exponentArry(exponentArry eps){
	string res;
	for(auto n :eps)
		res += std::to_string(n);
	return res;
}

//convert $z_i x_j$ into $x_i x_j^2$
exponentArry doubling_ex(exponent e){
	auto eps = unpack(e);
	exponentArry res;
//	std::cout << output_exponentArry(eps.first);
	for(int i=0; i<maxVar; i++)
		res[i] = eps.first[i]*2 + eps.second[i];
//	std::cout << output_exponentArry(res);
	return res;
}

//construct the list of exponents using the order of the ex-poly ex_index
void ex2poly(string filename, monomial_index &mon_index){
	std::fstream fs(filename, std::ios::out | std::ios::binary);
	int32_t nums = mon_index.mon_array.size();
	fs.write((char*)&nums, 4);
	for(auto e:mon_index.mon_array){
		exponentArry eps = doubling_ex(e);
	//	std::cout << output(e) << ":" << output_exponentArry(eps) << "\n";
		fs.write((char*)&eps, sizeof(exponentArry));
	}
}

//construct the list
int main(int argc, char** argv){
	int max_deg = std::atoi(argv[1]);
	monomial_index mon_index(max_deg);
	mon_index.init_mon_array();
	string filename = argv[1] ;
	ex2poly(filename  + "_"+ "ex2poly_index", mon_index);
}