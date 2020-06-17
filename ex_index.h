//mon_index.h
#pragma once
#include"ex_exponents.h"
#include"matrices.h"
#include<set>
#include<map>
#include<functional>

//data for the set of monomials in the dual steenrod algebra
class monomial_index{
	//set of monomials in a fixed degree
	std::map<int,std::set<exponent>> mon_sets;
public:
	//the maximal degree
	int max_degree;
	//array of the monomials ordered in degree
	std::vector<exponent> mon_array;
	//the index of the monoials
	std::map<exponent,matrix_index> mon_index;
	//number of monoials under a fixed degree
	std::vector<unsigned> ranksBelow;
    
	//constructor
	monomial_index(int maxdeg);
    
	//the set of monomials in a fixed degree
	std::set<exponent> monomial_sets(int);
	
	//generat all the monomials
	std::set<exponent> generate_monomials(int n);
    
	//the number of all monoials
	int number_of_all_mons();
	
	//the largest generator
	int max_var, max_ex;
    
	//initialization
	void init_mon_array();
	
	//compute the table for some substitution rule
	template<typename ring, typename output_type>
	void substitution_table(std::function<ring(int)> values, std::function<ring(int)> exvalues, std::iostream &writer, RingOp<ring> *ringop, std::function<output_type(const ring&)> rule, std::function<void(output_type const&, std::iostream&)> outputter){
		//the values of the powers of the varaibles
		std::vector<std::vector<ring>> monos(maxVar+1);
		for(int i=1; i<=maxVar; ++i){
			//compute the i-th variables
			monos[i].resize(xnMaxExpo[i]);
			monos[i][0] = ringop->unit(1);
			if(xnDeg(i)>max_degree) continue;
			ring val = values(i);
			for(int j=1; j*xnDeg(i)<=max_degree; ++j)
				monos[i][j] = ringop->multiply(monos[i][j-1], val);
		}
		
		std::cout << "computations of single factors complete\n" << std::flush;
		
		//compute all the monomials
		int counter = 0;
		for(auto e: mon_array){
			ring ev = ringop->unit(1);
			for(int k=1; k<=maxVar; ++k){
				//load the power of the value of the k-th varaible
				auto fct = xnVal(e,k);
				if(fct!=0)
					ev = ringop->multiply(ev,monos[k][fct]);
			}
			for(int k=1; k<=maxVar; ++k)
				//load the values of the k-th ex
				if(znVal(e,k))
					ev = ringop->multiply(ev, exvalues(k));
			
			std::cout << "\r" << counter++ << "/" << mon_array.size() << std::flush;
			//save the result
			outputter(rule(ev),writer);
		}
		std::cout << "table complete\n" << std::flush;
	}
	
	//compute the table for some substitution rule
	template<typename ring>
	void substitution_table(std::function<ring(int)> values, std::function<ring(int)> exvalues, std::iostream &writer, RingOp<ring> *ringop){
		std::function<ring(const ring&)> id = [](const ring &x){ 
			return x; };
		std::function<void(const ring&, std::iostream&)> outputter = [ringop](const ring &x, std::iostream &writer){
			ringop->save(x,writer); };
		
		substitution_table(values, exvalues, writer, ringop, id, outputter);
	}
	
	
	//change to vector notation
	template<typename base_ring>
	vectors<matrix_index,base_ring> poly2vec(polynomial<base_ring> const &x, ModuleOp<exponent,base_ring> *modoper){
		std::function<matrix_index(exponent)> rule = [this](exponent e){
			return mon_index[e]; };
		return modoper->re_index(rule,x);
	}
};
