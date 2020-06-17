//Fp.cpp
#include "Fp.h"

//constructor
Fp_Op::Fp_Op(int p) : prime(p){
	std::cout << "Fp operations initialized with p=" << prime << "\n" << std::flush;
}
    
//addition
Fp Fp_Op::add(const Fp& x, const Fp& y){ 
	return (x+y)%prime; }

//addition
Fp Fp_Op::add(Fp&& x, Fp&& y){
	return (x+y)%prime; }
    
//multiplication
	Fp Fp_Op::multiply(const Fp& x, const Fp& y){ 
		return (x*y)%prime; }
    
//uint from integers
Fp Fp_Op::unit(int n){ 
	if(n>=0) return n%prime;
		else return prime - unit(-n);
}

//check if it is zero
bool Fp_Op::isZero(const Fp& x){
	return x%prime==0;}
	
//the zero element
	Fp Fp_Op::zero(){ 
		return 0; }
		
//the negation
Fp Fp_Op::minus(const Fp& x){ 
	if(isZero(x)) 
		return 0; 
	else 
		return prime - (x%prime); 
}
    
//output into string
string Fp_Op::output(Fp x){ 
	return std::to_string(x); }
	
//save to stream
void Fp_Op::save(const Fp& x, std::iostream& writer){ 
	writer.write((char*)&x, sizeof(Fp)); }

//load from stream 
Fp Fp_Op::load(std::iostream& reader){ 
	Fp res; 
	reader.read((char*)&res,sizeof(Fp)); 
	return res; 
}
    
//check if invertibel
bool Fp_Op::invertible(const Fp& x){ 
	return !isZero(x); }

//inverse of an invertible element, only for p=2
Fp Fp_Op::inverse(const Fp& x){
	if(prime==2) return x;
	std::cerr << "not implemented!";
	return 0;
}
