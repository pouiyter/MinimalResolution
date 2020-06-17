//Fp.h
#pragma once
#include "algebra.h"

//an element in Fp is stored in a unsigned 16-bit integer
typedef ushort Fp;

//class for operations on Fp
class Fp_Op : public virtual RingOp<Fp>{
public:
	//the characteristic
	int const prime;
	
	//constructor
	Fp_Op(int);
    
	//additions
	Fp add(const Fp&, const Fp&);
	Fp add(Fp&&, Fp&&);
    
	//multiplication
	Fp multiply(const Fp&, const Fp&);
    
	//the unit map
	Fp unit(int);
	
	//check if it is 0
	bool isZero(const Fp&);
	//the zero element
	Fp zero();
	
	//the negation
	Fp minus(const Fp&);
    
	//IO operations
	string output(Fp);
	void save(const Fp&, std::iostream&);
	Fp load(std::iostream&);
    
	//check if invertibel
	bool invertible(const Fp&);
	//the inverse
	Fp inverse(const Fp&);
};
