//BP.h
#pragma once
#include"Z2.h"
#include"polynomial.h"
#include"matrices.h"
#include"mon_index.h"
#include"hopf_algebroid.h"
#include"BPQ.h"
#include<map>

//BP_* = Zp[v_1,v_2,\dots]
typedef polynomial<Z2> BP;
//BP_*BP = BP_*[t_1,t_2,\dots]
typedef polynomial<BP> BPBP;
typedef polynomial<BPBP> BPBPBP;
    
//operations for BP, with BPBP as a polynomial algebra over BP via the Right unit.
class BP_Op;

//operations on BP_*BP    
class BPBP_Op : public PolynomialOp_Para<BP>{
public:
	BPBP_Op(BP_Op*);
};
    
//operations on BPBPBP
class BPBPBP_Op : public PolynomialOp_Para<BPBP>{
public:
	BPBPBP_Op(BPBP_Op*);
};
    
//structur of Hopf algebroid on (BP_*,BP_*BP)
class BP_Op : virtual public Hopf_Algebroid<BP,BPBP>, public PolynomialOp_Para<Z2>{
public:
	//ring operations
	Z2_Op *Z2_oper;
	BPBP_Op BPBP_opers;
	BPBPBP_Op BPBPBP_opers;
        
	//module operations
	ModuleOp<matrix_index,BP> BPMod_opers;
	ModuleOp<matrix_index,BPBP> BPBPMod_opers;
	ModuleOp<matrix_index,Z2> Z2Mod_opers;

	//index of monomials
	monomial_index mon_index;
        
private:
	//table for eta_L
	matrix<BP> *etaL_table;
        
    //table for delta
	matrix<BPBP> *delta_table;
	
	//table for change from right to left and switch the place of vi and ti
	matrix<BP> *R2L_table;
        
public:
	BP_Op(int maxdeg, Z2_Op*, matrix<BP> *etaL_mat, matrix<BPBP> *delta_mat, matrix<BP> *R2L_mat);
        
	//load the structure maps from files
	void load_etaL(string);
	void load_R2L(string);
	void load_delta(string);

	//left unit
	BPBP etaL(const BP&);
	
	//change the notation from right unit presentation to left unit presentation
	BPBP etaL(const BPBP&);
        
	//right unit
	BPBP etaR(const BP&);
        
	//change the notation from right to left
	BPBP R2L(const BPBP&);
        
	//divide by p^n
	BP divide_power_p(const BP&,int);
	BPBP divide_power_p(const BPBP&,int);
	//divide powers of v1
	BPBP divide_v1(const BPBP&,int);
        
	//change from polynomial notation to vector notation
	vectors<matrix_index, BP> algebroid2vector(const BPBP&,int);
	
	//change from polynomial notation to vector notation
	vectors<matrix_index, BP> algebroid2vector(BPBP&&,int);
    
	//change from vector to polynoimial notation
	BPBP vector2algebroid(const vectors<matrix_index, BP>&);
    
	//the co-multiplication
	vectors<matrix_index, BPBP> delta(matrix_index);

	//number of generators below a degree
	unsigned ranksBelowDeg(unsigned); 
        
	//initilaizaiotn
	void initialize(string etaL_filename, string R2L_filename, string delta_filename);
        
	//lift elements in F2 to BP_*
	BP lift(F2);
	
	//lift vectors over F2 to vectors over BP_*
	vectors<matrix_index,BP> lift(const vectors<matrix_index,F2>&);
	
	//make the structure tables
	void make_tables(int maxVar, string R2Lfilename, string deltafilename, string etaL_filename, string R2L_filename, string delta_filename, std::ostream &outputfile=std::cout);
	
	//return the element h1
	BPBP h1();
	//return the element h2
	BPBP h2();
	//return the element h3
	BPBP h3();
	//return the element v1
	BP v1();
	//return the top thetas on the Moore spectrum
	std::vector<BPBP> thetas();
};
    
//cofree comodule for BP_*BP
typedef cofree_comodule<BPBP,int> FreeBPCoMod;
//generic comodules
typedef comodule_generic<BPBP,int> BPCoMod_generic;
    
