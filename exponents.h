//exponents.h
#pragma once
#include"polynomial.h"
#include<array>

//maximal number of variables
constexpr int maxVar = 7;
//exponent for multivarible monomials
typedef std::array<int,maxVar> exponentArry;
  
//multiply all the elements in a list
exponent multiply_all(exponent lis[], int start, int end);

//find the hash value of the n-th variable
std::vector<exponent> get_pos(exponent lis[], int length);
  
//the degree of the n-th polynomial generator
extern int xnDegs[maxVar+1];
//return the degree of the n-th generator
int xnDeg(int n);

//the maximal exponent of the n-th generator
extern exponent xnMaxExpo[maxVar+1];
//maximal total degree
extern exponent totalMax;

//the hash value of the n-th variable
extern std::vector<exponent> xnPos;

//the exponent for xn in e
int xnVal(exponent e, int n);

//the multivarible exponent for e
exponentArry unpack(exponent e);
  
//the total degree of e
int total_deg(exponent e, std::function<int(int)> varDeg = xnDeg);

//the exponent for a single varible
exponent singleVar(int n, int i);
//the exponent for the n-th variable
exponent vars(int n);

//construct an exponent
exponent pack(const int*);
exponent pack(const exponentArry&);
  
//output e
string output(exponent e, string name = "x");

//the previous exponent
std::pair<int,exponent> previous(exponent e);

//save the list of exponents
void save_expArry(exponentArry,std::iostream&);
//load the list of exponents
exponentArry load_expArry(std::iostream&);

//substitute values in a monomial
template<typename ring>
ring substitute(std::function<ring(int)> values, exponent e, RingOp<ring> *ringop){
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
	return res;
}

//substitute value into a polynomial
template<typename value_ring, typename base_ring>
value_ring substitute(std::function<value_ring(int)> values, std::function<value_ring(base_ring)> coef_rule, polynomial<base_ring> const &x, RingOp<value_ring>* ringop){
	std::function<value_ring(int)> summands = [&values, ringop, &coef_rule, &x](int k){
		//compute the value of an monomial
		auto a = substitute(values, x.dataArray[k].ind, ringop);
		
		auto b = coef_rule(x.dataArray[k].coeficient);
		
		return ringop->multiply(a,b);
	};
	
	//sum up the summands
	return sum(summands, 0, x.dataArray.size(), ringop);
}


//substitute value into a polynomial
template<typename value_ring, typename base_ring>
value_ring substitute(std::vector<value_ring> const &values, std::function<value_ring(base_ring)> coef_rule, polynomial<base_ring> const &x, RingOp<value_ring>* ringop){
	std::function<value_ring(int)> val = [&values](int i){
		return values[i]; };
	return substitute(val, coef_rule, x, ringop);
}

//substitute value into a polynomial
template<typename value_ring>
value_ring substitute(std::vector<value_ring> const &values, polynomial<value_ring> const &x, RingOp<value_ring>* ringop){
	std::function<value_ring(value_ring)> id = [](value_ring x){
		return x; };

		return substitute(values, id, x, ringop);
}
