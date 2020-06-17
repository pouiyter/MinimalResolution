//algNov.cpp
#include"algNov.h"

//the operators
Z2_Op *algNov_table::Z2Oper;
ModuleOp<matrix_index, Z2> *algNov_table::Modop;
BP_Op *algNov_table::BPoper;

//filtration of a cycle
int algNov_table::filtration(cycle_name cyc){
	return cyc.first[0];
}

//the filtration by number of v's
int algNov_table::v_valuation(Z2 r, exponent e){
	return Z2Oper->valuation(r) + total_deg(e, [](int){return 1;});
}
	
//transform a term into a cycle name
cycle_name algNov_table::naming(const Z2 &coeficient, matrix_index ind, primitive_data &P){
	//get the exponent
	exponent e = P[ind].coeficient;
	//construct the cycle
	std::vector<int> eo;
	//first compute the algebraic Novikov filtration
	eo.push_back(v_valuation(coeficient,e));
	//then get the number of v0's
	eo.push_back(Z2Oper->valuation(coeficient));
	//then the exponent
	for(int i=1; i<=maxVar; ++i)
		eo.push_back(xnVal(e,i));
	
	//the cycle
	return {eo, P[ind].gen_pos};
}

//transform a term into a cycle name, using Pcyc
cycle_name algNov_table::naming(const Z2 &coeficient, matrix_index ind){
	return naming(coeficient, ind, *Pcyc);
}
	
//computing the leading term of a cycle vector
unsigned algNov_table::leading_term(SS_entry<cycle_name, Z2>::value_type& v){
	//case of zero
	if(v.size()==0) return -1;
	
	//the current leading term
	unsigned led = 0;
	//the name of the i-th term
	auto cn = [&v,this](int i){
		return naming(v.dataArray[i].coeficient, v.dataArray[i].ind); };
	//find the leading term with the smallest cycle name
	for(unsigned i=1; i<v.size(); ++i){
		if(cn(i) < cn(led))
			led = i;
	}
	return led;
}
	
//set the operations
void algNov_table::set_op(BP_Op *BP_oper){
	BPoper = BP_oper;
	Modop = &BPoper->Z2Mod_opers;
	Z2Oper = BPoper->Z2_oper;
}
	
//constructors
algNov_table::algNov_table() : SS_table(&BPoper->Z2Mod_opers){}

//check the entry is a tag or not
bool algNov_table::tagged(const SS_entry<cycle_name,Z2> &et){
	return et.tag.second != Invalid;
}

cycle_name algNov_table::invalid(){
	std::vector<int> vd;
	return std::make_pair(vd,Invalid);
}

//construct the cycle pot. pot has additional leading number for sorting reasons...
std::set<cycle_name> algNov_table::cycle_pot(int pric){
	std::set<cycle_name> res;
	for(unsigned i=0; i<Ptag->size(); ++i){
		//construct ((l,0,v1^e1,...),gen_ind)
		auto s = naming(Z2Oper->unit(1), i, *Ptag);
		//compute the filtration
		int l = filtration(s);
		for(int j=0; j+l<=pric; ++j){
			//construct ((l+j, j, v1^e1, ...), gen_ind)
			s.first[0] = j+l;
			s.first[1] = j;
			res.emplace(s);
		}
	}	
	return res;
}

//output of cycle names
string algNov_table::output(cycle_name sn, int k){
	string res;
	for(unsigned i=1; i<sn.first.size(); ++i)
		if(sn.first[i]!=0)
			res += "v" + std::to_string(i-1) + "^" + std::to_string(sn.first[i]);
	res += "[" + std::to_string(k) + "-" + std::to_string(sn.second) + "]";
	return res;
}

//IO of entries
std::pair<std::pair<int,int>,string> algNov_table::output(SS_entry<cycle_name, Z2>& et, int k){
	string res;
	//if it is a untagged cycle
	if(!tagged(et)) 
		res = output(et.cycle,k);
	//if it is tagged, output cyc <- tag
	else 
		res = output(et.cycle,k) + "\t<-\t" + output(et.tag,k-1) +"\t|d" + std::to_string(diff_length(et)+1);
	//degree = (2t-s, s+i)
	std::pair<int,int> degs = std::make_pair(degree(et.cycle)*2-k, k+filtration(et.cycle));
	
	return std::make_pair(degs, res + "\t|deg=(" + std::to_string(degs.first) + "," + std::to_string(degs.second)+ ")");
}

//save the cycle  name
void algNov_table::save(cycle_name nm, std::iostream& fl){
	unsigned ls = nm.first.size();
	fl.write((char*)&ls, sizeof(unsigned));
	fl.write((char*)nm.first.data(), ls*sizeof(int));
	fl.write((char*)&nm.second, sizeof(int));
}

//save the entry
void algNov_table::save(SS_entry<cycle_name, Z2>& et, std::iostream& fl){
	save(et.tag,fl);
	save(et.cycle,fl);
	Modop->save(et.full_tag,fl);
	Modop->save(et.full_cycle,fl);
}

//load the cycle name
cycle_name algNov_table::load(std::iostream& fl){
	unsigned ls;
	fl.read((char*)&ls, sizeof(unsigned));
	cycle_name res;
	res.first.resize(ls);
	fl.read((char*)res.first.data(), ls*sizeof(int));
	fl.read((char*)&res.second, sizeof(int));
	return res;
}

//load the entry
void algNov_table::load(SS_entry<cycle_name, Z2>& et, std::iostream& fl){
	et.tag = load(fl);
	et.cycle = load(fl);
	et.full_tag = Modop->load(fl);
	et.full_cycle = Modop->load(fl);
}

//construct a vector using the primitive data
SS_entry<cycle_name, Z2>::value_type algNov_table::make_vec(cycle_name cyc, primitive_data& P, ModuleOp<matrix_index,Z2> *Modop, Z2_Op *Z2Oper){ 
	//construct the exponent
	auto e = pack(cyc.first.data()+2);
	//find the entry
	auto n = P.prim_index[{(matrix_index)cyc.second, e}];
	//construct the term v0^e0v1^e1...[gen_pos]
	return Modop->singleton(n,Z2Oper->power_p(cyc.first[1]));
}

//construct a full tag from the leading term
SS_entry<cycle_name, Z2>::value_type algNov_table::get_tag(cycle_name tag){
	return make_vec(tag, *Ptag, Modop, Z2Oper); }

//the degree of cycle
int algNov_table::degree(cycle_name cyc){
	//get the exponent of cyc = ((fil, e0,e1,...),gen_pos)
	exponent e = pack(cyc.first.data()+2);
	//find the position of v^e[gen_pos]
	int n = Pcyc->prim_index[{(matrix_index)cyc.second,e}];
	//get the degre
	return Pcyc->gen_deg[n];
}
    
//construct the table using the complex    
void  algNov_tables::table_of_complex(BPComplex& Comp, int pric){
	//initialize primitive data
	tables[0]->Pcyc = &Comp.Prims[0];
	for(unsigned i=0; i<Comp.size(); ++i){
		tables[i]->Ptag = &Comp.Prims[i-1];
		tables[i]->Pcyc = &Comp.Prims[i];
	}
	//construct tables
	for(unsigned i=1; i<Comp.size(); ++i){
		std::cout << "\n" << "table:" << i << "\n" << std::flush;
		tables[i]->make_table(pric-i,Comp.Maps[i],*tables[i-1]);
	}
}

//out put the algebraic Novikov table
string algNov_tables::output_tables(){
	string res;
	std::set<std::pair<std::pair<int,int>,string>> lst;
	for(unsigned i=0; i<tables.size(); ++i){
		auto jst = tables[i]->output_table(i,tables.size()-i-1);
		for(auto sm: jst)
			lst.insert(sm);
	}
	for(auto sm : lst)
		res+=sm.second+"\n";
	return res;
}

//save the tables
void algNov_tables::save(std::iostream& fl){
	for(auto tm:tables)
		tm->SS_table<cycle_name,Z2>::save(fl);
}

//load the tables
void algNov_tables::load(std::iostream& fl){
	for(auto tm:tables)
		tm->SS_table<cycle_name,Z2>::load(fl);
}

//save to file
void algNov_tables::save(string filename){
	std::fstream fk(filename, std::ios::out | std::ios::binary);
	save(fk);
}

//load from filte
void algNov_tables::load(string filename){
	std::fstream fk(filename, std::ios::in | std::ios::binary);
	load(fk);
}
    
//set the data of tables
void algNov_tables::set_table(std::vector<algNov_table>* tbs){
	tables.resize(tbs->size());
	for(unsigned i=0; i<tbs->size(); ++i)
		tables[i] = &tbs->at(i);
}

//check if a name is valid
bool algNov_table::valid(const cycle_name& cyc){
	return cyc.second != Invalid;
}

std::set<std::pair<std::pair<int,int>,string>> algNov_table::output_table(int k, int pric){
	return SS_table<cycle_name,Z2>::output(k,pric);
}
