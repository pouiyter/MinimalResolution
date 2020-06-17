//Qp.h
#pragma once
#include<vector>
#include<gmpxx.h>
#include"algebra.h"

//the type of the numerators
typedef mpz_class num_type;

//an p-rational number has a numerator and a valuation serving as the power of the denominator
typedef struct {
	num_type numerator;
	int16_t valuation;
} Qp;

//operations on Qp
class Qp_Op : public virtual RingOp<Qp>{
	//operations on the numerators
	num_type safe_add(num_type,num_type);
	num_type safe_mult(num_type,num_type);
public:
	//the residue characteristic
	virtual int prime()=0;
	
	//i-th power of p
	num_type power_p(unsigned i);
	int power_p_int(unsigned i);

	//kill off powers of p in numerator
	void simplify(Qp &x);
        
	//IO operations
	string output(Qp x);
	void save(Qp const &x, std::iostream& writer);
	Qp load(std::iostream& reader);
	
	//ring operations
	Qp add(Qp const &x, Qp const &y);
	Qp add(Qp &&x, Qp &&y);
	Qp zero();
	bool isZero(Qp const &x);
	Qp minus(Qp const &x);
	Qp multiply(Qp const &x, Qp const &y);
	Qp unit(int x);
	bool invertible (Qp const &x);
	Qp inverse (Qp const &x);
        
	//the constructor
	Qp construct(int num, int16_t val);
    
	//tranform of the numerator between char*
	typedef unsigned char uchar;
	std::vector<uchar> save2vector(num_type x);
	num_type load_vector(std::vector<uchar> const&);
        
	//tranform into an unsigned integer
	uint64_t int_part(Qp);
};
  
//spectialize to p=2
class Q2_Op : virtual public Qp_Op{
    int prime();
};

//the Qp operations with integral IO
class Qp_int : virtual public Qp_Op{
public:
	void save(Qp const &x, std::iostream& writer);
};

//Q2 operations with integral IO
class Q2_int : virtual public Qp_int, virtual public Q2_Op{
public:
	string output(Qp x);
};
