//Qp.cpp
#include"Qp.h"
#include<iostream>
#include<sstream>
#include<fstream>

//the i-th power of p
num_type Qp_Op::power_p(unsigned i){
	num_type res = 1;
	for(unsigned k=0;k<i;++k)
		res = safe_mult(res,prime());
	return res;
}

//i-th power of p
int Qp_Op::power_p_int(unsigned i){
	int res = 1;
	for(unsigned k=0;k<i;++k)
		res = res*prime();
	return res;
}

//addition of the numerator
inline num_type Qp_Op::safe_add(num_type x, num_type y){
    return x+y;
}

//multiplication for the numerator
inline num_type Qp_Op::safe_mult(num_type x, num_type y) {
    return x*y;
}

//subtract the power of p's
inline void Qp_Op::simplify(Qp &x) {
	//convention for 0
	if(x.numerator==0) {
		x.valuation = 0;
		return;
	}

	//return when the numerator is not divisible by p
	if(x.numerator % ((num_type)prime()) != (num_type)0) return;

	//divide by p
	x.numerator /= (num_type) prime();
	++x.valuation;

	simplify(x);
}

//output
string Qp_Op::output(Qp x) { 
    std::ostringstream fr;
    fr << x.numerator;
    return fr.str() + "(" + std::to_string(x.valuation) + ")"; 
}

//additions
inline Qp Qp_Op::add(Qp const &x, Qp const &y){
	//we need to make the denominators equal before addition
	if(x.valuation<y.valuation)
		return {safe_add(x.numerator, safe_mult(y.numerator,power_p(y.valuation-x.valuation))), x.valuation};
	
	if(x.valuation>y.valuation)
		return {safe_add(y.numerator, safe_mult(x.numerator,power_p(x.valuation-y.valuation))), y.valuation};
	
	//if the denominators are already equal
	Qp z = {safe_add(x.numerator,y.numerator), x.valuation};
	//check if the numerator is divisible by p
	simplify(z);
	return z;
}

//addition
inline Qp Qp_Op::add(Qp &&x, Qp &&y){ 
	return add(x,y); }

//zero elememnt
inline Qp Qp_Op::zero() { 
	return {(num_type)0,0}; }
	
//check if equal to zero
inline bool Qp_Op::isZero(Qp const &x){ 
	return x.numerator==(num_type)0; }
	
//the negative
inline Qp Qp_Op::minus(Qp const &x){ 
	return {-x.numerator, x.valuation}; }

//multiplication
inline Qp Qp_Op::multiply(Qp const &x, Qp const &y){ 
	return {safe_mult(x.numerator,y.numerator), (int16_t)(x.valuation+y.valuation)}; }

//unit map
inline Qp Qp_Op::unit(int x){ 
	Qp result = {(num_type)x,0};
	simplify(result);
	return result;
}

//check if an element is invertible
inline  bool  Qp_Op::invertible(Qp const &x){ 
	return !isZero(x); }

//the inverse operation is not needed in this project
inline Qp Qp_Op::inverse(Qp const &x) 
{ 
	std::cerr << "not implemented!"; abort(); }

//IO operations
inline void Qp_Op::save(Qp const &x, std::iostream& writer){
	//transform the numerator
	auto vr = save2vector(x.numerator);
	//write the valuation
	writer.write((char*)&x.valuation,2);
	//write the numerator
	int16_t l = vr.size();
	writer.write((char*)&l,2);
	writer.write((char*)vr.data(),l);
}

//IO operations
inline Qp Qp_Op::load(std::iostream& reader){
	Qp result;
	//read the valuation
	reader.read((char*)&result.valuation,2);
	//read the numerator
	int16_t l;
	reader.read((char*)&l,2);
	std::vector<unsigned char> vr(l);
	reader.read((char*)vr.data(),l);
	//transform the numerator
	result.numerator = load_vector(vr);
	return result;
}

//the characteristic of the residue field for Q2
inline int Q2_Op::prime(){ 
	return 2;}

//the constructor
Qp Qp_Op::construct(int num, int16_t val){ 
	return {(num_type)num,val}; }

//transform an element into an interger
uint64_t Qp_Op::int_part(Qp x){
	if(x.valuation<0){
		std::cout << output(x);
		std::cerr << "not integral!";
		abort();
	}
	num_type ip = x.numerator * power_p(x.valuation);
	//the maximal number for an 64-bit unsigned integer
	static const num_type t64 = power_p(64);
	//truncate the lower part
	ip = ip % t64;
	if(ip<0) ip+=t64;
	
	return (uint64_t) ip.get_ui();
}

//transform the numerator to a char*
inline std::vector<unsigned char> Qp_Op::save2vector(num_type x){
	std::vector<unsigned char> result;
	if(x==0){
		result.push_back(0);
		return result;
	}
	
	//save the sign
	if(x>0)
		result.push_back(1);
	else
		result.push_back(2);
	
	//save the 256-adic expression
	while(x!=0){
		num_type ns = x%256;
		result.push_back((unsigned char) ns.get_ui());
		x/=256;
	}
	
	return result;
}

//transform a char* into the numerator
inline num_type Qp_Op::load_vector(std::vector<unsigned char> const &rv){
	num_type res=0;
	
	num_type vl = 1;
	for(int i=1;i<(int)rv.size();++i){
		res += vl*(num_type)rv[i];
		vl *= 256;
	}
    
	//check the signe
	if(rv[0]==2)
		res=-res;
	return res;
}

void Qp_int::save(Qp const &x, std::iostream& writer){
	uint64_t w = this->int_part(x);
	writer.write((char*)&w, 8);
}

string Q2_int::output(Qp x){
	uint64_t w = this->int_part(x);
	return std::to_string((int)w);
}
