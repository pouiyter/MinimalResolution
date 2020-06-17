//ex_exponents.h
//exponents for the ring F[x_1,x_2,\dots]\otimes\Lambda(\xi_1,\xi_2,\dots)
//the maximal valid degree is 253
#pragma once
#include"polynomial.h"
#include<array>

//maximal number of variables
constexpr int maxVar = 7;
//exponent for multivarible monomials
typedef std::array<int,maxVar> exponentArry;
//exponent for the exterior variables
typedef std::array<bool,maxVar> exArry;

//class for the exponents for ex-poly rings
class ex_poly{
public:
	exponent e;
public:
	//constructor
	ex_poly(exponent);
	ex_poly();
	//add exponents
	ex_poly operator+(const ex_poly) const;
	//compare
	bool operator<(const ex_poly) const;
	bool operator>(const ex_poly) const;
	bool operator==(const ex_poly) const;
	bool operator>=(const ex_poly) const;
	//add, but not used
	ex_poly operator+(const int) const;
};

string to_string(ex_poly);

//add two exponents
exponent add(exponent, exponent);

//check if two exponents have common exterior variable
bool common_ex(exponent, exponent);

//multiply all the elements in a list
exponent multiply_all(exponent lis[], int start, int end);

//find the hash value of the n-th polynomial variable
std::vector<exponent> get_pos(exponent lis[], int length);
//find the hash value of the n-th exterior variable
std::vector<exponent> get_ex(exponent lis[], int length);

//the degree of the n-th polynomial generator
extern int xnDegs[maxVar+1];
//the degree of the n-th exterior generator
extern int znDegs[maxVar+1];

//return the degree of the n-th exterior generator
int xnDeg(int n);
//return the degree of the n-th exterior generator
int znDeg(int n);

//the maximal exponent of the n-th generator
extern exponent xnMaxExpo[maxVar+1];
//maximal total degree
extern exponent totalMax;

//the hash value of the n-th polynomial variable
extern std::vector<exponent> xnPos;
//the hash value of the n-th exterior variable
extern std::vector<exponent> znPos;

//the exponent for xn in e
int xnVal(exponent e, int n);
//the exponent for zn in e
bool znVal(exponent e, int n);

//the multivarible exponent for e
std::pair<exponentArry,exArry> unpack(exponent e);
  
//the total degree of e
int total_deg(exponent e, std::function<int(int)> varDeg = xnDeg, std::function<int(int)> exDeg = znDeg);

//the exponent for a single varible
exponent singleVar(int n, int i);
//the exponent for the n-th variable
exponent vars(int n);
//the exponent for the n-th exterior variable
exponent ex(int n);

//construct an exponent
exponent pack(const int*, const bool*);
  
//output e
string output(exponent e, string name = "x", string exname = "z");

//the previous exponent
std::pair<int,exponent> previous(exponent e);

//substitute values in a monomial
template<typename ring>
ring substitute(std::function<ring(int)> values, std::function<ring(int)> exvalues, exponent e, RingOp<ring> *ringop){
	ring res = ringop->unit(1);
	for(int i=1; i<=maxVar; ++i){
		//the value of the i-th factor
		ring res_i;
		if(xnVal(e,i)!=0)
			res_i = ringop->power(values(i), xnVal(e,i));
		else
			res_i = ringop->unit(1);
		
		res = ringop->multiply(res, res_i);
	}
	for(int i=1; i<=maxVar; ++i)
		if(znVal(e,i))
			res = ringop->multiply(res, exvalues(i));
	return res;
}

//substitute value into a polynomial
template<typename value_ring, typename base_ring>
value_ring substitute(std::function<value_ring(int)> values, std::function<value_ring(int)> exvalues, std::function<value_ring(base_ring)> coef_rule, polynomial<base_ring> const &x, RingOp<value_ring>* ringop){
	std::function<value_ring(int)> summands = [&values, &exvalues, ringop, &coef_rule, &x](int k){
		//compute the value of an monomial
		auto a = substitute(values, exvalues, x.dataArray[k].ind, ringop);
		
		auto b = coef_rule(x.dataArray[k].coeficient);
		
		return ringop->multiply(a,b);
	};
	
	//sum up the summands
	return sum(summands, 0, x.dataArray.size(), ringop);
}


//substitute value into a polynomial
template<typename value_ring, typename base_ring>
value_ring substitute(std::vector<value_ring> const &values, std::vector<value_ring> const &exvalues, std::function<value_ring(base_ring)> coef_rule, polynomial<base_ring> const &x, RingOp<value_ring>* ringop){
	std::function<value_ring(int)> val = [&values](int i){
		return values[i]; };
	std::function<value_ring(int)> exval = [&exvalues](int i){
		return exvalues[i]; };
	return substitute(val, exval, coef_rule, x, ringop);
}

//substitute value into a polynomial
template<typename value_ring>
value_ring substitute(std::vector<value_ring> const &values, std::vector<value_ring> const &exvalues, polynomial<value_ring> const &x, RingOp<value_ring>* ringop){
	std::function<value_ring(value_ring)> id = [](value_ring x){
		return x; };

		return substitute(values, exvalues, id, x, ringop);
}

//class for polynomials tensored with exterior algebra
template<typename ring>
class ExPolyOp_Para : public PolyOp_Para<ex_poly,ring>{
public:
	//multiply with monomials
	poly<ex_poly,ring> mon_multiply(const poly<ex_poly,ring> &x, ex_poly e, const ring &r){
		poly<ex_poly,ring> result;
		
		for(unsigned i=0; i<x.size(); ++i){
			// z^2=0 for exterior variables
			if(common_ex(x.dataArray[i].ind.e, e.e))
				continue;
			//nontrivial multiply
			typename vectors<ex_poly,ring>::term new_term;
			new_term.ind = x.dataArray[i].ind + e;
			new_term.coeficient = this->ringOper->multiply(x.dataArray[i].coeficient, r);
			result.push(std::move(new_term));
		}
		return result;
	}
	
	//constructor
	ExPolyOp_Para(RingOp<ring> *ringop) : ModuleOp<ex_poly,ring>(ringop),  PolyOp_Para<ex_poly,ring>(ringop){}
};
