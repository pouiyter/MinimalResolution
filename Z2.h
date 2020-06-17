//Z2.h
#pragma once
#include"algebra.h"
#include"Fp.h"

//we use 64-bit unsigned integer to denote 2-adic numbers
typedef uint64_t Z2;

//the residue field
typedef Fp F2;

class Z2_Op : public virtual RingOp<Z2>{
public:
	//the constructor
	Z2_Op();
	
	//the charateristic of the residue field
	int prime();
	
	//the operations on F2
	Fp_Op F2_opers;
	
	//the valuation
	unsigned valuation(Z2);
	
	//n-th power of p
	Z2 power_p(int n);
	
	// divide by p^n
	Z2 divide(const Z2&, int n); 
	
	//lift from the residue field
	Z2 lift(F2);
    
	//addition
	Z2 add(const Z2&, const Z2&);
	Z2 add(Z2&&, Z2&&);
    
	//multiplication
	Z2 multiply(const Z2&, const Z2&);
    
	//the unit
	Z2 unit(int);
    
	//check if it is zero, or rather divisible by 2^64
	bool isZero(const Z2&);
	
	//the zero element
	Z2 zero();
	
	//negation
	Z2 minus(const Z2&);
    
	//IO operations
	string output(Z2);
	void save(const Z2&, std::iostream&);
	Z2 load(std::iostream&);
    
	//check if invertible
	bool invertible(const Z2&);
	//inverse of an invertible element
	Z2 inverse(const Z2&);
};
