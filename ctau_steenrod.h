//ctau_steenrod.h
#pragma once
#include "Fp.h"
#include "polynomial.h"
#include "matrices.h"
#include "ex_index.h"
#include "hopf_algebroid.h"

//The dual Ctau Steenrod algebra
typedef poly<ex_poly,Fp> P;
//The tensor square of the Steenrod algebra
typedef poly<ex_poly,P> PP;

//The Hopf algebra of the dual steenrod algebra
class Steenrod_Op;

//the ring operations on the dual steenrod algebra
class P_Op : public ExPolyOp_Para<Fp>{
public:
	//constructor
	P_Op(Steenrod_Op*);
	//construct a monomial
	P construct(int k, int i);
	//construt a monomial on exterior
	P construct_ex(int);
};

//ring operations on the tensor square of the dual steenrod algebra
class PP_Op : virtual public ExPolyOp_Para<P>{
public:
	//constructor
	PP_Op(P_Op*);
	//operations on the dual steenrod algebra
	P_Op *P_opers;
	//construct a monomial
	PP construct(int k1, int i1, int k2, int i2);
	//construct exr \otimes poly
	PP construct_expoly(int,int,int);
};
    
//the class of the Hopf algebra of dual steenrod algebra
class Steenrod_Op : virtual public Hopf_Algebroid<Fp,P>, virtual public Fp_Op{
public:
	//ring operations
	P_Op P_opers;
	//ring operations on P\otimes P
	PP_Op PP_opers;
        
	//operations on Fp-modules
	ModuleOp<matrix_index,Fp> FpMod_opers;

	//operations on P-modules
	ModuleOp<matrix_index,P> PMod_opers;
        
	//compute the p-th power
	std::function<int(int)> power_p;

	//the index of the monomials ordered in total degree
	monomial_index mon_index;

	//the constructor
	Steenrod_Op(int maxdeg, int, matrix<P>*);
        
	//load the data of the co-multiplications
	void load_delta(std::iostream&);
	
	//construct the data of the co-multiplications
	void make_delta(std::iostream&);
	
	//the data for the co-multiplications
	matrix<P> *delta_table;
                    
	//transform an element of P to a vector over Fp
	vectors<matrix_index, Fp> algebroid2vector(const P&,int);
	//the inverse of the transformation
	P vector2algebroid(const vectors<matrix_index, Fp>&);
    
	//the left unit
	P etaL(const Fp&);
	//the right unit
	P etaR(const Fp&);
    
	//the co-multiplications
	vectors<matrix_index, P> delta(matrix_index);

	//the number of monomials below a certain degree
	unsigned ranksBelowDeg(unsigned); 
        
	//the initialization
	void initialize(std::iostream &delta_data);
	
	//output the generators of a resolution
	string output_resolution(std::iostream &gens_file, int length);
	
	//output the generators of a resolution
	string output_resolution(string gens_file, int length);
	
private:
	static int id_int(int);
};

//the class of cofree comoules over P
typedef cofree_comodule<P,int> FreeSteenrodCoMod;

//a comodule over P
typedef comodule_generic<P,int> SteenrodCoMod_generic;
