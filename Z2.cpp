//Z2.cpp
#include"Z2.h"
#include<sstream>

//the constructor
Z2_Op::Z2_Op() : F2_opers(2) {}

//the residual characteristic
inline int Z2_Op::prime(){ return 2; }

//addition
inline Z2 Z2_Op::add(const Z2& x, const Z2& y){ return x+y; }

//addition
inline Z2 Z2_Op::add(Z2&& x, Z2&& y){ return x+y; }

//multiplication
inline Z2 Z2_Op::multiply(const Z2& x, const Z2& y){ return x*y; }

//unit
inline Z2 Z2_Op::unit(int x){ return x; }

//check if zero
inline bool Z2_Op::isZero(const Z2& x){ return x==(Z2)0; }

//the zero element
inline Z2 Z2_Op::zero(){ return 0; }

//negation
inline Z2 Z2_Op::minus(const Z2& x){ return -x;}

//output using hex form
string Z2_Op::output(Z2 x){
	std::stringstream res;
	res << std::hex << x;
	return res.str();
}

//write to a 8-byte. This would depend on the cpu archetecture
inline void Z2_Op::save(const Z2& x, std::iostream& writer){ writer.write((char*)&x, 8); } 

//load from a 8-byte
inline Z2 Z2_Op::load(std::iostream& reader){ 
	Z2 result;
	reader.read((char*)&result,8);
	return result;
}

//if invertible
bool Z2_Op::invertible(const Z2& x){ return (x&(Z2)1) != 0; }

//the inverse, using a recursion formula
Z2 Z2_Op::inverse(const Z2& x){
	Z2 u = 1;
	for(int i=0;i<6;++i)
		u=u*(2-x*u);
	return u;
}

//the 2-adic valuation
unsigned Z2_Op::valuation(Z2 x){
	//convention for 0
	static constexpr int max_val = 35536;
	if(isZero(x)) return max_val;
	
	//valuation of invertible is 0
	if(invertible(x)) return 0;
	//for non-invertible elements
	return valuation(x/prime()) + 1;
}

//p^n
Z2 Z2_Op::power_p(int n){
	if(n==0) return unit(1);
	Z2 p = unit(prime());
	if(n>0)
		return power_p(n-1)*p;
	//undefined if n<0
	std::cerr << "not integral!";
	return 0;
}

//divide by p. The missing term casued by truncation is guessed by the sign
Z2 Z2_Op::divide(const Z2& x, int n){
	int64_t y=x;
	return y/(int64_t)power_p(n);
}

//lift from residue field
Z2 Z2_Op::lift(F2 x){
    return x;
}
