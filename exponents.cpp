//exponents.cpp
#include "exponents.h"

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

//the degree of the n-th variable
int xnDegs[maxVar+1] = {0,1,3,7,15,31,63,127};
int xnDeg(int n){ 
	return xnDegs[n]; }

//the maximal value of the n-th exponent
exponent xnMaxExpo[maxVar+1] = {0, 255, 85, 37, 17, 9, 5, 3};

//the maximal total degree
exponent totalMax = multiply_all(xnMaxExpo,1,maxVar+1);

//the hash-value of the n-th varaible
std::vector<exponent> xnPos = get_pos(xnMaxExpo,maxVar+1);

//the value of the exponent on the n-th variable
int xnVal(exponent e, int n){ 
	return (e/xnPos[n-1]) % xnMaxExpo[n]; }

//unpack the informations contained in e
exponentArry unpack(exponent e){
	exponentArry res;
	for(int i=0; i<maxVar; i++)
		res[i]=xnVal(e,i+1);
	return res;
}

//the degree of e
int total_deg(exponent e, std::function<int(int)> varDeg){
	int res = 0;
	for(int i=1; i<=maxVar; i++)
		res+=varDeg(i)*xnVal(e,i);
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

//output an exponent
string output(exponent e, string name){
	auto arr = unpack(e);
	string res = "";
	for(int i=0;i<maxVar;++i)
		if(arr[i]!=0)
			res+= name + "_" + std::to_string(i+1) + "^" + std::to_string(arr[i]);
	return res;
}

//the previous exponent, by diminishing the first non-trivial variable
std::pair<int,exponent> previous(exponent e){
	for(int i=0; i<maxVar; i++)
		if(xnVal(e,i+1)>0)
			return std::pair<int,exponent>(i+1,e-vars(i+1));
	return std::pair<int,exponent>(0,0);
}
  
//find the hash-valus of some multi-variable exponent
exponent pack(const int *eps){
	exponent res = 0;
	for(int i=1; i<=maxVar; ++i)
		res+=singleVar(i,eps[i-1]);
	return res;
}

//find the hash-valus of some multi-variable exponent
exponent pack(const exponentArry &eps){
	exponent res = 0;
	for(int i=1; i<=maxVar; ++i)
		res+=singleVar(i,eps[i-1]);
	return res;
}


//save the list of exponents
void save_expArry(exponentArry eps,std::iostream &writer){
	writer.write((char*)&eps, sizeof(exponentArry));
}
//load the list of exponents
exponentArry load_expArry(std::iostream &reader){
	exponentArry res;
	reader.read((char*)&res,sizeof(exponentArry));
	return res;
}