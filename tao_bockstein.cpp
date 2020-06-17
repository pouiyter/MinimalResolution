//tao_bockstein.cpp
#include"tao_bockstein.h"

//output the set of generators
void output_generators(int resolution_length, string filename_generators, string filename_outs){
//open the files
	std::fstream gens_file(filename_generators, std::fstream::in | std::fstream::binary);
	std::fstream out_file(filename_outs, std::fstream::out);
	
	for(int i=0; i<=resolution_length; ++i) {
		//read the number of generators
		int32_t M_rank; gens_file.read((char*)&M_rank, 4);
		//load the data of the cofree comodule
		cofree_comodule<motSteenrod, MotDegree> F0;
		F0.load(gens_file);
		
		//outpu t the degree of the generators
		for(unsigned k=0; k<F0.generators.rank; ++k)
			out_file << "[" + std::to_string(i) + "-" + std::to_string(k) +"]\t:\t" << "(" << F0.generators.degree[k].deg - i << "," << F0.generators.degree[k].weight << ")" << std::endl;
	}
}

//make the complex of the primiteves
void motComplex::load(int resolution_length, string filename_generators, string filename_resolution){
	//initialize the data
	terms.clear();
	maps.clear();
	maps_data.clear();
	terms.resize(resolution_length+1);
	maps.resize(resolution_length);
	maps_data.resize(resolution_length);
	for(unsigned i=0; i<resolution_length; ++i)
		maps[i] = &maps_data[i];
	
	//open the files
	std::fstream gens_file(filename_generators, std::fstream::in | std::fstream::binary);
	std::fstream res_file(filename_resolution, std::fstream::in | std::fstream::binary);

	//load the generators
	std::vector<cofree_comodule<motSteenrod, MotDegree>> F(resolution_length+1);
	for(int iee=0; iee<=resolution_length; ++iee){
		//load the data for the generators
		int32_t M_rank; 
		gens_file.read((char*)&M_rank, 4);
		cofree_comodule<motSteenrod, MotDegree> F0;
		F0.load(gens_file);

		terms[iee] = F0.generators;
		F[iee] = F0;
	}
	
	//construct the maps
	for(int iee=0; iee<resolution_length; ++iee){
		matrix_mem<tauPoly> M;
		M.load(res_file);
		maps[iee]->set_rank(F[iee].position_of_gens.size());
		for(unsigned k=0; k<F[iee].position_of_gens.size(); ++k){
			maps[iee]->insert(k,M.find(F[iee].position_of_gens[k]));
		}
	}
}

//output the complex
void motComplex::output(string filename_outs){
	std::fstream out_file(filename_outs, std::fstream::out);
	
	for(unsigned i=0; i<maps.size(); ++i)
		out_file << maps[i]->output(std::to_string(i), std::to_string(i+1));
}

//the leading term of a vector
int tau_table::leading_term(vectors<matrix_index, tauPoly> const& v){
	int res =0;
	for(unsigned k=1;k<v.size();k++)
		if(tau_oper.tau_valuation(v.dataArray[res].coeficient)>tau_oper.tau_valuation(v.dataArray[k].coeficient))
			res=k;
	return res;
}

//construct the bockstein table 
void tau_table::make_pretable(std::stack<int> &pot, matrix<tauPoly> *M, tau_table &result, tau_table &pre_table){
	//do nothing if the pot is already cleared
	if(pot.empty()) return;
	//deal with the current leading term
	int cur = pot.top();    
	pot.pop();
	
	//construct an elemwnt with the current leading term
	auto ft = tau_module_oper.singleton(cur);
	//find the current image 
	auto fc = M->maps_to(ft);

	//trying to simplify fc
	while(!tau_module_oper.isZero(fc)){
		//get the leading term of fc
		auto ld = leading_term(fc);
		auto li = fc.dataArray[ld].ind;
		//try to find the leading term of fc in the table
		auto it = result.cycle_index.find(li);
		
		//compute the filtration of fc
		int val = tau_oper.tau_valuation(fc.dataArray[ld].coeficient);
		
		//if the leading term of fc is not in the table, then update
		if(it == result.cycle_index.end()){
			//add the entry: li <- cur, diff_length=val
			result.insert(cur,li,ft,fc,val);
			//go the the next term
			make_pretable(pot,M,result,pre_table);
			return;
		}
		
		//if the leading_term term is found
		int ledind = it->second;
		//if the full tag of the founded entry in the current table has lower filtration than the cur, then the current table need to be modified
		if(val<result.table[ledind].diff_length || (val==result.table[ledind].diff_length && cur > result.table[ledind].tag )){
			//remove the tag of the found entry, and send it back to the pot
			pot.push(result.table[ledind].tag);
			//modify the table
			result.modify(ledind,cur,li,ft,fc,val);
			//deal with the next leading term
			make_pretable(pot,M,result,pre_table);
			return;
		}
		
		//otherwise, the current vector will be modified to have higher filtration
		//pv = tau^(val - diff_length)
		tauPoly pv = val - result.table[it->second].diff_length;
		//construct the tag to be subtracted
		auto mdf = tau_module_oper.scalor_mult(pv,result.table[ledind].full_tag);
		mdf = tau_module_oper.minus(mdf);
		//construct the cycle to be subtracted
		auto mfc = M->maps_to(mdf);
		//modify the current vectors
		ft = tau_module_oper.add(std::move(ft),std::move(mdf));
		fc = tau_module_oper.add(std::move(fc),std::move(mfc));
	}
	
	//if in the end fcbecomes zero, then fc is a cycle
	pre_table.insert(tau_table_entry::Invalid, cur,fc,ft, -1);
	
	//deal with next case
	make_pretable(pot,M,result,pre_table);
	return;
}

//insert a new entry
void tau_table::insert(int tag, int cycle, vectors<matrix_index, tauPoly> const &full_tag, vectors<matrix_index, tauPoly> const &full_cycle, int diff_length){
	int n=table.size();
	table.push_back({tag,cycle,full_tag,full_cycle,diff_length});
	tag_index.emplace(tag,n);
	cycle_index.emplace(cycle,n);
}

//modify the table
void tau_table::modify(int n, int tag, int cycle, vectors<matrix_index, tauPoly> const &full_tag, vectors<matrix_index, tauPoly> const &full_cycle, int diff_length){
	//delete the old entry
	cycle_index.erase(table[n].cycle);
	tag_index.erase(table[n].tag);
	
	//update the table
	table[n].tag = tag;
	table[n].cycle = cycle;
	table[n].diff_length = diff_length;
	table[n].full_cycle = full_cycle;
	table[n].full_tag = full_tag;
	//update the index
	tag_index.emplace(tag,n);
	cycle_index.emplace(cycle,n);
}

//make the pot for the next table
std::stack<int> tau_table::pot_maker(int rank) const{
	//construt the initial pot.
	std::set<int> psd;
	for(int i=0;i<rank;++i)
		psd.insert(i);
	//erase those which are already cycles
	for(auto k: table)
		psd.erase(k.cycle);
	//construct the pot
	std::stack<int> result;
	for(auto k : psd) 
		result.push(k);
	return result;
}

//make the tables from the complex
std::vector<tau_table> make_table(motComplex& CMP){
	//initial table
	std::vector<tau_table> res_table;
	res_table.resize(CMP.size());
	
	//construct the tables
	for(int i=1;i<CMP.size();++i){
		std::cout << i << std::flush;
		//make the pot
		auto pot = res_table[i-1].pot_maker(CMP.terms[i-1].rank);
		std::cout << pot.size() << "\n" << std::flush;
		//construt the table for the i-th map
		tau_table::make_pretable(pot, CMP.maps[i-1], res_table[i], res_table[i-1]);
	}
	return res_table;
}

//output the entries
string tau_table_entry::output(int nm) const{
	if(tag==Invalid)
		return "{" + std::to_string(nm) + "-" + std::to_string(cycle) + "}\t:cycle";
	return "t^"+ std::to_string(diff_length) + "{" + std::to_string(nm) + "-" + std::to_string(cycle) + "}\t<-\t{" + std::to_string(nm-1)+"-"+std::to_string(tag)+"}";
}

//output the table
string tau_table::output(int nm) const{
	string res;
	for(auto et: table)
		res += et.output(nm) + "\n";
	return res;
}

//output the tables
string output_tables(const std::vector<tau_table>& Tb){
	string res;
	for(unsigned i=0;i<Tb.size();++i)
		res+=Tb[i].output(i);
	return res;
}


//get the cycles in a table
cycle_data tau_table::get_cycles() const{
	cycle_data res;
	for(auto &itm : table){
		//include the nontrivial cycles
		if(itm.tag == tau_table_entry::Invalid)
			res.emplace(itm.cycle,itm.full_cycle);
		else{
			//modify the boundary by dividing powers of tau
			tauPoly fc = tau_oper.power_tau(-itm.diff_length);
			auto nt = tau_module_oper.scalor_mult(fc, itm.full_cycle);
			res.emplace(itm.cycle, nt);
		}
	}
	return res;
}


//get the tags 
cycle_data tau_table::get_tags() const{
	cycle_data res;
	for(auto &itm : table){
		if(itm.tag != tau_table_entry::Invalid)
			res.emplace(itm.tag,itm.full_tag);
	}
	return res;
}

//output the data of the cycles
string output(cycle_data const &cd, int nm){
	string res;
	for(auto &it : cd)
		res +=  "{" + std::to_string(nm) + "-" + std::to_string(it.first) + "}" + "\t=\t" + tau_module_oper.output(it.second,std::to_string(nm)) + "\n";
	return res;
}

//get the cycle data from a table
void make_cycle_tables(std::vector<tau_table> const &tab, std::vector<cycle_data> &result){
	result.resize(tab.size());
	for(unsigned i=0; i<tab.size(); ++i)
		result[i] = tab[i].get_cycles();
	cycle_data tg;
	for(unsigned i=1; i<tab.size(); ++i){
		tg = tab[i].get_tags();
		result[i-1].insert(tg.begin(),tg.end());
	}
}

//transform a cycle into a sum of generators in a cycle
vectors<matrix_index, tauPoly> find_cycle(cycle_data& table, vectors<matrix_index, tauPoly> v){
	vectors<matrix_index, tauPoly> res;
	while(!tau_module_oper.isZero(v)){
		//find the leading term of v
		auto ld = tau_table::leading_term(v);
		auto cf = v.dataArray[ld].coeficient;
		auto ind = v.dataArray[ld].ind;
		//subtract the generator from v
		v = tau_module_oper.add(v,tau_module_oper.scalor_mult(cf, table.at(ind)));
		//record the generator
		res = tau_module_oper.add(res, tau_module_oper.singleton(ind, cf));
	}
	return res;
}

//output the cycle in terms of 
string output_cycles(int nm, const vectors<matrix_index,tauPoly> &v){
	string res = "o";
	for(auto tm : v.dataArray){
		res += "+t^" + std::to_string(tm.coeficient) + "{" + std::to_string(nm) + "-" + std::to_string(tm.ind) + "}";
	}
	return res;
}

//construct the multiplication table
void multiplication_table(motSteenrod const& x, int x_deg, string filename_generators, string filename_res_truns, string filename_outs, int resolution_length, MotSteenrodOp& ms_oper, std::vector<cycle_data>& cyc_table){
	//construcut the basic multiplication table
	matrix_mem<tauPoly>  mult;
	ms_oper.make_multiplication_table(x, x_deg, &mult);
	
	//open the data files
	std::fstream resmaps_file(filename_res_truns, std::fstream::in | std::fstream::binary);
	std::fstream gens_file(filename_generators, std::fstream::in | std::fstream::binary);
	std::fstream out_file(filename_outs, std::fstream::out);
	
	for(int iee=0; iee<resolution_length-1; ++iee){
		//read the data of the generators
		int32_t M_rank; 
		gens_file.read((char*)&M_rank, 4);
		cofree_comodule<motSteenrod, MotDegree> F0;
		F0.load(gens_file);
		//construct the multiplication matrix
		matrix_mem<tauPoly>  mp,mmt;
		mp.load(resmaps_file);
		for(unsigned i=0; i<F0.generators.rank; ++i) 
			if((int)F0.generators.degree[i].deg + x_deg <= (int)ms_oper.maxDeg){
				auto des = mp.maps_to(F0.multiply_using_table(F0.position_of_gens[i],mult));
				mmt.insert(i,des);
			}
		//compute the multiplication in terms of chosen generators	
		for(unsigned i=0; i<F0.generators.rank; ++i) 
				if((int)F0.generators.degree[i].deg + x_deg <= (int)ms_oper.maxDeg){
					auto des = mmt.maps_to(cyc_table[iee].at(i));
					auto mres = find_cycle(cyc_table[iee+1],des);
					out_file << "{" + std::to_string(iee) + "-" + std::to_string(i) + "}\t->\t" + output_cycles(iee+1, mres) + "\n";
				}
	}
}
