 #pragma once

#include "hopf_algebroid.h"
#include "exponents.h"
#include "Fp.h"

typedef Fp F2;

//a homogenious polynomial in tau
typedef int16_t tauPoly;
//operations on F2[tau]  
class tauOper : public virtual RingOp<tauPoly>{
	//the zero polynomial
	constexpr static tauPoly internal_zero = (1<<15);
public:
	inline bool isZero(tauPoly const &x);
public:
	inline tauPoly add(tauPoly const &x, tauPoly const &y);
	inline tauPoly add(tauPoly &&x, tauPoly &&y);
	inline tauPoly zero();
	inline tauPoly minus(tauPoly const &x);
	inline tauPoly multiply(tauPoly const &x, tauPoly const &y);
	inline tauPoly unit(int n);
	inline bool invertible(tauPoly const &x);
	inline tauPoly inverse(tauPoly const &x);
    
	//IO operations
	void save(tauPoly const&, std::iostream&);
	tauPoly load(std::iostream&);
	inline string output(tauPoly x);
	
    //the number of taus
	int tau_valuation(tauPoly);
	//shift the exponent of tau
	void shift(tauPoly&,int);
	//the power of tau
	tauPoly power_tau(int);
};
//operations on Fp[tau]  
extern tauOper tau_oper;

//the dual Steenrod algebra
typedef polynomial<tauPoly> motSteenrod;
//A \otimes A
typedef polynomial<motSteenrod> mStmSt;

//the operations on the motivic Steenrod algebra
class MotSteenrodRingOp : virtual public PolynomialOp<tauPoly>{
public:
	MotSteenrodRingOp(tauOper *tauop);
};

//ring operations
extern MotSteenrodRingOp motSteenrod_oper;
extern PolynomialOp<motSteenrod> mStmSt_oper;

//the topological degree and the weight
class MotDegree{
public:
	int deg, weight;
	//adding the bidegree
	static MotDegree add(MotDegree x, MotDegree y);
    //the bidegree of the varaibles
	static int varDeg(int n);
	static int varWeight(int n);
    //degree of the monomials
	static int expoDeg(exponent e);
	static int expoWeight(exponent e);
    //output the monomial
	string output();
	//compare
	bool operator<(const MotDegree);
};

//cofree comodules
typedef cofree_comodule<motSteenrod,MotDegree> FreeMotCoMod;
  
//operations on modules
extern ModuleOp<matrix_index,motSteenrod> motSteenrod_module_oper;  
extern ModuleOp<matrix_index,tauPoly> tau_module_oper;

//the Hopf algebroid of the motivic Steenrod algebra
class MotSteenrodOp : virtual public Hopf_Algebroid<tauPoly,motSteenrod>{
	//the list of monmials
    std::vector<exponent> mon_array;
	//the index of monoials
	std::map<exponent, matrix_index> mon_index;
public:
	//the matrix for the coactions
	matrix<motSteenrod> *cofree_coaction;
    
	//constructor
	MotSteenrodOp(matrix<motSteenrod> *cofa, int max_degree);
    
	//the tau-adic valuation in motivic Steenrod algebra
	int tauVal(exponent e);
    
	//the degree of the cofree comodule
	MotDegree cofree_degs(matrix_index n);
    
	//turn a element in A to a vector
	vectors<matrix_index,tauPoly> algebroid2vector(motSteenrod const &x,int);
    
	//transform a vector back to A
	motSteenrod vector2algebroid(vectors<matrix_index,tauPoly> const &x);
	//transform a polynomial in A into a vector in A
	vectors<matrix_index,motSteenrod> algebroid2vectorSt(mStmSt const &x);
    
	//generate data for the coactions
	void generate_cofree_coaction(string coaction_filename, string expo_filename);

	//left and right units
	motSteenrod etaR(const tauPoly &x);
	motSteenrod etaL(const tauPoly &x);
    
	//number of monomials
	unsigned ranksBelowDeg(unsigned n);
	std::vector<unsigned> ranksbelow;
	
	//co-multiplication
	vectors<matrix_index, motSteenrod> delta(matrix_index n);
  
	//construct the elements hi
    motSteenrod hi(int);
	
	//the comultiplication
	std::function<mStmSt(exponent)> delta_rule;
	
	//initialize the list of monoials
	void init_mon_array(string filename);
	
	//output the array of monomials
	string output_monomials();
	
	//lift elements in F2 to  F2[tau]
	tauPoly lift(F2 x);
	
	//lift vectors over F2 to vectors over F2[tau]
	vectors<matrix_index,tauPoly> lift(const vectors<matrix_index,F2>&v);
};
