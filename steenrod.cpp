//steenrod.cpp
#include "steenrod.h"

//constructor
Steenrod_Op::Steenrod_Op(int maxdeg, int prime, matrix<P> *deltaTable) : Fp_Op(prime), P_opers(this), PP_opers(&P_opers), FpMod_opers(this), PMod_opers(&P_opers), mon_index(maxdeg){
	//initialize the monomial set
	mon_index.init_mon_array();
	
	this->maxDeg = maxdeg;

	Hopf_Algebroid<Fp,P>::ringOper = this;
	algebroidRingOper = &P_opers;
    
	moduleOper = &FpMod_opers;
	algebroidModuleOper = &PMod_opers;
	
	delta_table = deltaTable;

	//p-th power
	power_p = [prime](int i){
		int res = 1;
		for(int k=0;k<i;++k)
			res *= prime;
		return res;
	};
}

//constructor of ring operations on P
P_Op::P_Op(Steenrod_Op *St_op) : ModuleOp<exponent,Fp>(St_op), PolynomialOp_Para<Fp>(St_op){}

//constructor of ring operations on P\otimes P
PP_Op::PP_Op(P_Op *P_op) : ModuleOp<exponent,P>(P_op), PolynomialOp_Para<P>(P_op){
	P_opers=P_op; }

//load the file of the co-multiplication
void Steenrod_Op::load_delta(std::iostream &reader){
	delta_table->rank = mon_index.number_of_all_mons();
	for(int i=0; i<mon_index.number_of_all_mons(); ++i){
		auto dels = PP_opers.load(reader);
		
		//transform the polynomial on the dual Steenrod algebra into a vector
		std::function<matrix_index(exponent)> rd = [this](exponent e){
			return mon_index.mon_index[e]; };
		delta_table->insert(i, PP_opers.re_index(rd, dels));
	}
}

//construct the data of co-multiplication
void Steenrod_Op::make_delta(std::iostream &writer){
	//the co-multiplication of the generators
	std::function<PP(int)> deltaXn = [this](int j){
		auto result = PP_opers.add(PP_opers.construct(0,0,j,1),PP_opers.construct(j,1,0,0));
		for(int k=1; k<j; ++k)
			result = PP_opers.add(std::move(result),PP_opers.construct(k,1,j-k,power_p(k)));

		return result;
	};
	//construct the table
	mon_index.substitution_table(deltaXn, writer, &PP_opers);
}

//convert an element in the steenrod algebra into a vector
vectors<matrix_index, Fp> Steenrod_Op::algebroid2vector(const P& x,int shift){
	std::function<matrix_index(exponent)> rd = [this,shift](exponent e){
		return mon_index.mon_index[e]+shift; };
	return P_opers.re_index(rd,x);
}

//convert a vector to an element in the dual steenrod algebra
P Steenrod_Op::vector2algebroid(const vectors<matrix_index, Fp>& v){
	std::function<exponent(matrix_index)> rd = [this](matrix_index n){
		return mon_index.mon_array[n]; };
	return FpMod_opers.re_index(rd,v); 
}

//the co-multiplication
vectors<matrix_index, P> Steenrod_Op::delta(matrix_index n){
	return delta_table->find(n); }

//the left unit
P Steenrod_Op::etaL(const Fp& x){
	return P_opers.monomial(0,x); }

//the right unit
P Steenrod_Op::etaR(const Fp& x){
	return P_opers.monomial(0,x); }

//rank up to a degree
unsigned Steenrod_Op::ranksBelowDeg(unsigned n){ return mon_index.ranksBelow[n]; }

//initialization
void Steenrod_Op::initialize(std::iostream &delta_data){
	load_delta(delta_data);
	std::cout << "Steenrod co-product data loaded\n" << std::flush;
	
	std::function<int(matrix_index)>  cofree_degs = [this](int n){
		exponent e = mon_index.mon_array[n];
		return total_deg(e,xnDeg); };
	
	this->init_cofree_data(cofree_degs);
}

//construct x_k^i
P P_Op::construct(int k, int i){
	return this->monomial(singleVar(k,i));
}

//construct x_{k1}^{i1} \otimes x_{k2}^{i2}
PP PP_Op::construct(int k1, int i1, int k2, int i2) {
	return this->monomial(singleVar(k2,i2), P_opers->construct(k1,i1));
}

//the degree of comodules
template<>
int FreeSteenrodCoMod::underlyingDeg(int i){
	return i; }
	
//add degrees
template<>
int FreeSteenrodCoMod::add_degree(int a, int b){
	return a+b; }

//output resolution
string Steenrod_Op::output_resolution(std::iostream &gens_file, int length){
	std::vector<FreeSteenrodCoMod> gens(length);
	for(int i=0; i<length-1; ++i){
		int32_t M_rank;  
		gens_file.read((char*)&M_rank, 4);
		gens[i].load(gens_file);
	}
	
	string result = "\n";

	for(int i=length-1; i>=0; --i)
		result += gens[i].generators.output() + "\n";
	
	result += "\n";
	
	std::vector<std::vector<int>> dims(maxDeg, std::vector<int>(length));
	
	for(int i=0; i<length; ++i)
		for(auto d: gens[i].generators.degree)
			dims[d-i][i]++;
	
	for(int i=length-1; i>=0; --i){
		for(int j=0; j<(int)maxDeg; ++j)
			result += std::to_string(dims[j][i]) + "\t";
		result += "\n";
	}
	
	return result;
}

//output resolution
string Steenrod_Op::output_resolution(string gens_file, int length){
	std::fstream gf(gens_file, std::ios::in | std::ios::binary);
	
	return output_resolution(gf,length);
}