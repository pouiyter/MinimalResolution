//BPQ.h
#pragma once
#include"polynomial.h"
#include"Qp.h"
#include"exponents.h"

//The ring BP_*\otimes\mathbb{Q}
typedef polynomial<Qp> BPQ;
//BP_*BP\otimes\mathbb{Q}
typedef polynomial<BPQ> BPBPQ;
//BP_*BP\otimes BP_*BP\otimes\mathbb{Q}
typedef polynomial<BPBPQ> BPBPBPQ;

//the Hopf algebroid BP_*BP\otimes\mathbb{Q}
class BPQ_Op;
    
//ring operations on BPBPQ
class BPBPQ_Op : public PolynomialOp_Para<BPQ>{
	BPQ_Op* BPQ_opers;
public:
	BPBPQ_Op(BPQ_Op*);
	
	//construct c v_n^s1 t_k^s2
	BPBPQ construct(int n, int s1, int k, int s2, Qp coef);
};

//ring operations on BPBPBPQ
class BPBPBPQ_Op : public PolynomialOp_Para<BPBPQ>{
	BPBPQ_Op* BPBPQ_opers;
public:
	BPBPBPQ_Op(BPBPQ_Op*);

	//construct vn^s[t_k1^i1|t_k2^i2]
	BPBPBPQ construct(int n, int s, int k1, int i1, int k2, int i2, Qp coef);
};

//the Hopf algebroid BPQ
class BPQ_Op : virtual public PolynomialOp_Para<Qp>{
public:
	//ring operations
	Qp_Op *Qp_oper;
	BPBPQ_Op BPBPQ_oper;
	BPBPBPQ_Op BPBPBPQ_oper;
	
	Qp_int *Qpint_oper;
	PolynomialOp_Para<Qp> BPint_oper;
	PolynomialOp_Para<BPQ> BPBPint_oper;
	PolynomialOp_Para<BPBPQ> BPBPBPint_oper;
	
	//maximal number of varaibles
	int maxVar;
    
	//the expression for li in terms of vn
	std::vector<BPQ> li_table;

	//the expression for vn in terms of li
	std::vector<BPQ> vn_table;
        
	//the table for etaR
	std::vector<BPBPQ> etaR_table;
        
	//the table for the co-multiplication
	std::vector<BPBPBPQ> delta_table;

	//the constructor
	BPQ_Op(Qp_Op*, Qp_int*, int);
        
	//tansform between li expression and vi expression
	BPQ l2v(const BPQ&);
	BPQ v2l(const BPQ&);
        
	BPBPQ l2v(const BPBPQ&);
	BPBPQ v2l(const BPBPQ&);
        
	BPBPBPQ l2v(const BPBPBPQ&);
	BPBPBPQ v2l(const BPBPBPQ&);
	
	//left unit
	BPBPQ etaL(const BPQ&);
        
	//right unit in terms of li
	BPBPQ etaR_l(const BPQ&);
	//right unit in terms of vn
	BPBPQ etaR_v(const BPQ&);
        
	//the co-multiplication
	BPBPBPQ delta_l(exponent e);
	BPBPBPQ delta_v(exponent e);
	
	//output the formulas
	string show_li();
	string show_vn();
	string show_etaR();
	
	//output the formulas in binary file
	void output_R2L(std::iostream &writer);
	void output_delta(std::iostream &writer);
};
