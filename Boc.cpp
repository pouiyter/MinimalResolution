//Boc.cpp
#include"Boc.h"

int Boc_table::v_valuation(Z2 r, exponent){
	return Z2Oper->valuation(r);
}

void Boc_tables::set_table(std::vector<algNov_table>*){
	std::cerr << "action not allowed!\n";
}

void Boc_tables::set_table(std::vector<Boc_table>* tbs){
	tables.resize(tbs->size());
	Btables.resize(tbs->size());
	for(unsigned i=0; i<tbs->size(); ++i){
		tables[i] = &tbs->at(i);
		Btables[i] = &tbs->at(i);
	}
}

int Boc_table::num_v(cycle_name cyc){
	int res = 0;
	for(unsigned i=1; i<cyc.first.size(); ++i)
		res += cyc.first[i];
	return res;
}

std::pair<std::pair<int,int>,string> Boc_table::output(SS_entry<cycle_name, Z2>& et, int k){
	string res;
	//if it is a untagged cycle
	if(!tagged(et)) 
		res = this->algNov_table::output(et.cycle,k);
	//if it is tagged, output cyc <- tag
	else 
		res = this->algNov_table::output(et.cycle,k) + "\t<-\t" + this->algNov_table::output(et.tag,k-1) + "\t|d" + std::to_string(diff_length(et));
	
	std::pair<int,int> degs = std::make_pair(degree(et.cycle)*2-k, k+filtration(et.cycle));
	return std::make_pair(degs, res + "\t|deg=" + std::to_string(degs.first));
}

std::set<std::pair<std::pair<int,int>,string>> Boc_table::output_table(int k, int pric){
	std::set<std::pair<std::pair<int,int>,string>> lst;
	for(auto tm : *this){
		if(filtration(tm.cycle)<pric){
			//skip those entries with v0
			if(!tagged(tm)){
				if(filtration(tm.cycle)>0)
					continue;
			}
			else if(filtration(tm.tag)>0)
				continue;
			//output those esential entries
			lst.emplace(output(tm,k));
		}
	}
	return lst;
}


//transform Bockstein name into Algebraic Novikov name
multiplication_table<cycle_name> Boc_table::Bname2Aname(algNov_table &Atable, int pric, bool simple){
	//initialize the resulting table
	multiplication_table<cycle_name> res;
	//go through all the entries
	for(auto &tm : *this){
		//skip those tagged entries
		if(tagged(tm)) continue;
		//skip those entries divizible by v0 if we just want simple tables
		if(simple && tm.cycle.first[0]>0) continue;
		//read the full cycle
		auto v = tm.full_cycle;
		//construct the new entry
		multiplication_table_entry<cycle_name> nm = {tm.cycle, Atable.name_of_cycle(v, pric)};
		res.push_back(nm);
	}
	return res;
}

//transform Bockstein names into Algebraic Novikov names
std::vector<multiplication_table<cycle_name>> Boc_tables::Bname2Anames(int resolution_length, algNov_tables& AT, int pric){
	std::vector<multiplication_table<cycle_name>> result;
	for(int i=0; i<resolution_length; ++i){
		result.push_back(Btables[i]->Bname2Aname(*(AT.tables[i]), pric-i)); 
	}
	return result;   
}

//constructor
Boc_table::Boc_table() : SS_table(&algNov_table::BPoper->Z2Mod_opers){}
