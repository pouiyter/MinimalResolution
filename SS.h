//SS.h
#pragma once
#include"matrices.h"

//class for the entries in a spectral sequence
template<typename cycle_name, typename ring>
class SS_entry{
public:
	//type of the full cycles
	typedef vectors<matrix_index, ring> value_type;
	
	//leading term of an entry
	cycle_name tag, cycle;
	//full cycle and full tags
	value_type full_tag, full_cycle;
};
    
//the class for a spectral sequence
template<typename cycle_name, typename ring>
class SS_table : public std::vector<SS_entry<cycle_name, ring>>{
public:
	//operations on the vectors
	ModuleOp<matrix_index, ring> *Modop;
        
	//the index of tags and cycles
	std::map<cycle_name, unsigned> tag_index, cycle_index;
	
	//the filtration function
	virtual int filtration(cycle_name)=0;
	
	//the length of the differential
	int diff_length(SS_entry<cycle_name, ring> et){
		return filtration(et.cycle) - filtration(et.tag);
	}
	
	//transform a term into a cycle name
	virtual cycle_name naming(const ring &coeficient, matrix_index ind)=0;
	
	//computing the leading term of a cycle vector
	virtual unsigned leading_term(typename SS_entry<cycle_name, ring>::value_type&)=0;

	//inset an entry
	void insert(SS_entry<cycle_name,ring> et){
		auto i = this->size();
		//add the new entry in the end
		this->push_back(et);
		//update the indeces
		tag_index.emplace(et.tag,i);
		cycle_index.emplace(et.cycle,i);
	}
        
	//types of cycles
	enum cycle_type {Boundary, Untagged, NonFound};
        
	//the default invalid entry
	int const Invalid;
	virtual cycle_name invalid()=0;
	
	//the constructors
	SS_table(ModuleOp<matrix_index, ring> *modoper) : Invalid(-65537){
		Modop = modoper;	}
		
	//check the entry is a tag or not
	virtual bool tagged(const SS_entry<cycle_name,ring>&)=0;
	//check if an entry is valid
	virtual bool valid(const cycle_name&)=0;
	
	//simplify a cycle by subtraction boundaries using the data in the table
	cycle_type simplify(typename SS_entry<cycle_name, ring>::value_type& cyc, typename SS_entry<cycle_name, ring>::value_type& htpy, int pric){
		//if cycle is zero, then it is a boundary
		if(Modop->isZero(cyc)) return Boundary;
		
		//compute the leading term
		unsigned led = leading_term(cyc);
		cycle_name led_name = naming(cyc.dataArray[led].coeficient, cyc.dataArray[led].ind);
		
		//if the filtration exceeds the required pricision, then regard it as a boundary
		if(filtration(led_name)>pric) 
			return Boundary;
		
		//if not trivial, try to find it in the table
		auto it = cycle_index.find(led_name);
		if(it==cycle_index.end()) return NonFound;
		
		//if the searched entry is not tagged, then it is a nontrivial cycle
		auto pas = it->second;
		if(this->at(pas).tag.second==Invalid) return Untagged;
		
		//then we subtract the tagged entry from the cycle
		auto nc = Modop->minus(this->at(pas).full_cycle);
		auto nt = Modop->minus(this->at(pas).full_tag);
		
		cyc = Modop->add(std::move(cyc),std::move(nc));
		htpy = Modop->add(std::move(htpy),std::move(nt));
		
		//then simplify the subtracted cycle
		return simplify(cyc,htpy,pric);
	}
        
	//clear the table
	void clearing(){
		tag_index.clear();
		cycle_index.clear();
		this->clear();
	}
       
	//output entry
	virtual std::pair<std::pair<int,int>,string> output(SS_entry<cycle_name, ring>&, int)=0;
	
	//output of cycle names
	virtual string output(cycle_name,int)=0;
	
	//IO of entries
	virtual void save(SS_entry<cycle_name, ring>&, std::iostream&)=0;
	virtual void load(SS_entry<cycle_name, ring>&, std::iostream&)=0;
	
	//IO operations
	void save(std::iostream& fl){
		unsigned ls = this->size();
		fl.write((char*)&ls,sizeof(unsigned));
		for(auto &tm: *this){
			save(tm,fl);
		}
	}
	
	//IO operations
	void load(std::iostream& fl){
		clearing();
		unsigned ls;
		fl.read((char*)&ls, sizeof(unsigned));
		SS_entry<cycle_name, ring> tm;
		for(unsigned i=0; i<ls; ++i){
			load(tm,fl);
			insert(tm);
		}
	}
	
	//output the table, k indicates the homological degree
	virtual std::set<std::pair<std::pair<int,int>,string>> output(int k, int pric){
		std::set<std::pair<std::pair<int,int>,string>> lst;
		for(auto tm : *this){
			if((tm.tag.second==Invalid || diff_length(tm)!=0) && filtration(tm.cycle)<pric)
				lst.emplace(output(tm,k));
		}
		return lst;
	}
	
	//construct a full tag from the leading term
	virtual typename SS_entry<cycle_name, ring>::value_type get_tag(cycle_name)=0;

	//construct the cycle pot. pot has additional leading number for sorting reasons...
	virtual std::set<cycle_name> cycle_pot(int pric)=0;
	
	//filter the pot to delete those boundaries
	void filter_pot(std::set<cycle_name>  &pot, SS_table<cycle_name,ring> &prev_table){
		for(auto &et : prev_table){
		if(tagged(et))
			pot.erase(et.cycle);
		}
	}
	
	//construct the table form a set of candidate of tags
	void make_table(int pric, matrix<ring>* M, SS_table& T){
		//construct the cycle_pot
		auto pot = cycle_pot(pric);
		//filter the cycle_pot, using the previoud table
		filter_pot(pot, T);
		
		//construct the Curtis table. Note we start with tags of highest filtration
		for(auto it = pot.rbegin(); it != pot.rend(); it++){
			//get the canditate of the cycle
			auto vs = get_tag(*it);
			//get the boundary of the canditate
			auto bs = M->maps_to(vs);

			//try to add lower terms to make it a cycle
			auto BT = simplify(bs,vs,pric);
			
			//make the new entry
			SS_entry<cycle_name, ring> neE;
			//if the candidate can be made a cycle, then insert into the previous table
			if(BT==Boundary){
				neE = {invalid(), *it, Modop->zero(), vs}; 
				T.insert(neE);
			} //if we find a nontrivial tag, then insert into the current table
			else{
				//find the leading term if the boundary
				unsigned led = leading_term(bs);
				//make the new entry
				auto bn = naming(bs.dataArray[led].coeficient, bs.dataArray[led].ind);
				neE = {*it, bn, vs, bs};
				insert(neE);
			}
		}
	}
	
	//Return the set of cycles modulo filtration pric boundaries and cycles, original vector should be discarded after
	std::tuple<std::vector<cycle_name>, std::vector<cycle_name>, std::vector<cycle_name>> name_of_cycle(typename SS_entry<cycle_name, ring>::value_type& v, SS_table *next_table, int pric){
		//the set of tags, boundaries and cycles in v
		std::vector<cycle_name> tgs, bdrs, cycs;
		while(!Modop->isZero(v)){
			//find the leading term
			unsigned led = leading_term(v);
			cycle_name leading_name = naming(v.dataArray[led].coeficient, v.dataArray[led].ind);
			//if arrive at the desired pricision then return
			if(filtration(leading_name)>=pric) break;
			
			//find the entry for the leading term,
			auto it = cycle_index.find(leading_name);
			
			if(it == cycle_index.end()){
				//when the leading term is not a cycle, try it with the next table
				if(next_table == NULL)
					std::cerr << "invalid table!";
				else{
					//try to find the leading term as a tag
					auto ot = next_table->tag_index.find(leading_name);
					if(ot == next_table->tag_index.end())
						std::cerr << "invalid table!";
					else{
						if(filtration(next_table->at(ot->second).cycle)<pric)
							std::cerr << "invalid table!";
						//subtract the full tag
						v = Modop->add(v, Modop->minus(next_table->at(ot->second).full_tag));
						//add the information of the tag
						tgs.push_back(leading_name);
					}
				}
			}
			else {
				//when find it in the cycle list, modify the vector
				v = Modop->add(v, Modop->minus(this->at(it->second).full_cycle));
				//if it is untagged, mark as cycle, otherwise mark as boundary
				if(!tagged(this->at(it->second))) 
					cycs.push_back(leading_name);
				else 
					bdrs.push_back(leading_name);
			}
		}
		return std::make_tuple(tgs,bdrs,cycs);
	}
    
    //combine the tag and cycles
    std::pair<std::vector<cycle_name>, std::vector<cycle_name>> combine_cycles(std::tuple<std::vector<cycle_name>, std::vector<cycle_name>, std::vector<cycle_name>> tc){
		std::vector<cycle_name> cycs = std::get<0>(tc);
		for(auto tm : std::get<2>(tc))
			cycs.push_back(tm);
		return std::make_pair(std::get<1>(tc), cycs);
	}
    
	//Return the set of boundaries and cycles, original vector should be discarded after
	std::pair<std::vector<cycle_name>, std::vector<cycle_name>> name_of_cycle(typename SS_entry<cycle_name, ring>::value_type& v, int pric=64){
		auto cycls = name_of_cycle(v,NULL,pric);
		return std::make_pair(std::get<1>(cycls), std::get<2>(cycls));
	}
};
