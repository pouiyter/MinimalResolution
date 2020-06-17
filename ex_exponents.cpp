//ex_exponents.cpp
#include "ex_exponents.h"

//add two exponents
exponent add(exponent a, exponent b){
	//the polynomial part = (a+b) mod 2^25
	unsigned poly_part = (a+b) & ((1 << 25) - 1);
	//the exterior part = (a xor b) & (- 2^25)
	unsigned ex_part = (a^b) & (-(1 << 25));
	
	return poly_part + ex_part;
}

//check if two exponents have common exterior variable
bool common_ex(exponent a, exponent b){
	a = a & (-(1 << 25));
	b = b & (-(1 << 25));
	return a & b;
}

//multiply all the elements in a list
exponent multiply_all(exponent lis[], int start, int end){ 
	exponent res = 1;
	for(int i=start;i<end;++i) res*=lis[i];
	return res;
}

//find the hash value of the n-th variable
std::vector<exponent> get_pos(exponent lis[], int length){
	std::vector<exponent> res(length+1,1);
	for(int i=1;i<length;i++)
		res[i] = res[i-1]*lis[i];
	return res;
}

//find the hash value of the n-th exterior variable
std::vector<exponent> get_ex(int length){
	std::vector<exponent> res(length+1,1);
	for(int i=1; i<length; ++i)
		res[i] = 1 << (24+i);
	return res;
}


//the degree of the n-th variable
int xnDegs[maxVar+1] = {0,2,6,14,30,62,126,254};
int xnDeg(int n){ 
	return xnDegs[n]; }
	
//the degree of the n-th exterior generator
int znDegs[maxVar+1] = {0,1,3,7,15,31,63,127};
int znDeg(int n){
	return znDegs[n]; };

//the maximal value of the n-th exponent
exponent xnMaxExpo[maxVar+1] = {0, 128, 43, 19, 9, 5, 3, 1};

//the maximal total degree
exponent totalMax = multiply_all(xnMaxExpo,1,maxVar+1);

//the hash-value of the n-th varaible
std::vector<exponent> xnPos = get_pos(xnMaxExpo,maxVar+1);
//the hash value of the n-th exterior variable
std::vector<exponent> znPos = get_ex(maxVar+1);

//the value of the exponent on the n-th variable
int xnVal(exponent e, int n){ 
	e = e & ((1 << 25) - 1);
	return (e/xnPos[n-1]) % xnMaxExpo[n]; }
	
//the exponent for zn in e
bool znVal(exponent e, int n){
	return e & znPos[n]; }

//unpack the informations contained in e
std::pair<exponentArry,exArry> unpack(exponent e){
	exponentArry respol;
	exArry resex;
	for(int i=0; i<maxVar; i++)
		respol[i] = xnVal(e,i+1);
	for(int i=0; i<maxVar; ++i)
		resex[i] = znVal(e,i+1);
	return std::make_pair(respol, resex);
}

//the degree of e
int total_deg(exponent e, std::function<int(int)> varDeg, std::function<int(int)> exDeg){
	int res = 0;
	for(int i=1; i<=maxVar; i++)
		res+=varDeg(i)*xnVal(e,i);
	for(int i=1; i<=maxVar; ++i)
		if(znVal(e,i))
			res+=exDeg(i);
	return res;
}  
  
//the exponent for a single variable
exponent singleVar(int n, int i){ 
	if(n==0) return 0;
	return xnPos[n-1]*i; 
}

//the exponent for the n-th varaible
exponent vars(int n){ 
	if(n==0) return 0;
	return xnPos[n-1]; 
}

//the exponent for the n-th exterior variable
exponent ex(int n){
	return znPos[n];
}

//output an exponent
string output(exponent e, string name, string exname){
	auto arr = unpack(e);
	auto arrpol = arr.first;
	auto arrex = arr.second;
	string res = "";
	for(int i=0;i<maxVar;++i)
		if(arrex[i])
			res+= exname + "_" + std::to_string(i+1);
	for(int i=0;i<maxVar;++i)
		if(arrpol[i]!=0)
			res+= name + "_" + std::to_string(i+1) + "^" + std::to_string(arrpol[i]);
	return res;
}

//the previous exponent, by diminishing the first non-trivial variable
std::pair<int,exponent> previous(exponent e){
	for(int i=0; i<maxVar; i++)
		if(xnVal(e,i+1)>0)
			return std::pair<int,exponent>(i+1,e-vars(i+1));
	for(int i=0; i<maxVar; ++i)
		if(znVal(e,i+1))
			return std::pair<int,exponent>(-(i+1),e-ex(i+1));
	return std::pair<int,exponent>(0,0);
}
  
//find the hash-valus of some multi-variable exponent
exponent pack(const int *eps, const bool *exs){
	exponent res = 0;
	for(int i=1; i<=maxVar; ++i)
		res+=singleVar(i,eps[i-1]);
	for(int i=1; i<=maxVar; ++i)
		if(exs[i-1])
			res+=ex(i);
	return res;
}

//constructor
ex_poly::ex_poly(exponent e0){
	e = e0; }

//constructor
ex_poly::ex_poly(){
	e = 0; }

//add exponents
ex_poly ex_poly::operator+(const ex_poly y) const{
	return add(e,y.e); }
  
//compare
bool ex_poly::operator<(const ex_poly y) const{
	return e<y.e; }
	
//compare
bool ex_poly::operator>(const ex_poly y) const{
	return y<*this; }
  
//compare
bool ex_poly::operator==(const ex_poly y) const{
	return e==y.e; }

//compare
bool ex_poly::operator>=(const ex_poly y) const{
	return e>=y.e; }
	
//add, but not used
ex_poly ex_poly::operator+(const int y) const{
	std::cerr << "adding on Hash values;\n";
	return {e+y}; 
}
  
string to_string(ex_poly e){
	return std::to_string(e.e);
}
  
  
  
