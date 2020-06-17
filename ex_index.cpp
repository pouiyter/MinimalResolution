//mon_index.cpp
#include "ex_index.h"

//constructor
monomial_index::monomial_index(int maxdeg){ 
	//the monomials in degree 0
	mon_sets.emplace(0,std::set<exponent>({0})); 
	//set the maximal degree
	max_degree = maxdeg;
	
	//initialize the largest generator
	max_var = 0;
	for(int i=1; i<=maxVar; ++i)
		if(xnDeg(i)<max_degree)
			max_var = i;
	max_ex = 0;
	for(int i=1; i<=maxVar; ++i)
		if(znDeg(i)<max_degree)
			max_ex = i;
}

//generate the monimials in degre n
std::set<exponent> monomial_index::generate_monomials(int n){
	std::set<exponent> res;
	
	//run through all the variables
	for(int i=1;i<=maxVar;i++){
		//add one to the exponent of the i-th variables
		auto prev_set = monomial_sets(n-xnDeg(i));
		//suppose the previous sets of monomials are already computed
		for(auto it=prev_set.begin(); it!=prev_set.end(); it++)
			res.emplace(*it+vars(i));
	}
	//for the exterior variables
	for(int i=1; i<=maxVar; ++i){
		//add one exterior variables
		auto prev_set = monomial_sets(n-znDeg(i));
		//construct the current n
		for(auto &it : prev_set)
			if(!znVal(it,i))
				res.emplace(it+ex(i));
	}
	return res;
}  

//the number of all monomials
int monomial_index::number_of_all_mons(){ 
	return mon_array.size(); }

//initialization
void monomial_index::init_mon_array(){
	mon_array.clear();
	mon_index.clear();
      
	unsigned ranks = 0;
	//construct the set of monomials of each degree
	for(int d=0; d<=max_degree; d++){
		//construct all the monomials in degree d
		auto m_set = monomial_sets(d);
		//save all the monomials
		for(auto e: m_set)
			mon_array.push_back(e);
		//compute the number of monomials under degree d
		ranks += m_set.size();
		ranksBelow.push_back(ranks);
	}
	
	//construct the index of the monomials
	for(int i=0; i<(int)mon_array.size(); ++i)
		mon_index.emplace(mon_array[i],i);
}

//return the set of monomials in degree n
std::set<exponent> monomial_index::monomial_sets(int n){
	if(n<0) return std::set<exponent>();
	//try to find the data in mon_sets
	if(mon_sets.find(n)!=mon_sets.end()) return mon_sets[n];
	else{
		auto res = generate_monomials(n);
		mon_sets.emplace(n,res);
		return res;
	}
}
