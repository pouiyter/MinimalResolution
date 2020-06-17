//tao_bockstein.cpp
#include "tao_bockstein.h"
#include<assert.h>
#include<iostream>
#include<string>

void TauBockstein::multiplication_table(motSt const& x, int x_deg, string filename_generators, string filename_res_truns, string filename_outs, int resolution_length, MSOper& ms_oper){
        Matrices::matrix_array<MotSteentrod::tauPoly>  mult(1,&tau_module_oper, &tau_module_oper);
        ms_oper.make_multiplication_table(x, x_deg, &mult);
        
        std::fstream resmaps_file(filename_res_truns, std::fstream::in | std::fstream::binary);
        std::fstream gens_file(filename_generators, std::fstream::in | std::fstream::binary);
        std::fstream out_file(filename_outs, std::fstream::out);
        assert(resmaps_file.is_open() && gens_file.is_open() && out_file.is_open());
		
        int32_t M_rank; gens_file.read((char*)&M_rank, 4);
        HopfAlgebroid::cofree_comodule<motSt, MotDegree> F0;
        F0.load(gens_file);
        
        for(int iee=0; iee<resolution_length; ++iee) {
            
        Matrices::matrix_array<MotSteentrod::tauPoly>  mp(1,&tau_module_oper, &tau_module_oper);
 //       std::cout << F0.rank() << std::flush;
        mp.load(F0.rank(),resmaps_file);
           
        for(unsigned i=0; i<F0.generators.rank; ++i) 
                if((int)F0.generators.degree[i].deg + x_deg <= (int)ms_oper.maxDeg) {
                    auto des = mp.maps_to(F0.multiply_using_table(F0.position_of_gens[i],mult));
                    
                    out_file << "[" + std::to_string(iee) + "-" + std::to_string(i) + "]\t->\t" + tau_module_oper.output(des, std::to_string(iee+1)) + "\n";
           }
           
            gens_file.read((char*)&M_rank, 4);
            F0.load(gens_file);
        }
} 

void TauBockstein::output_generators(int resolution_length, string filename_generators, string filename_outs){
        std::fstream gens_file(filename_generators, std::fstream::in | std::fstream::binary);
        std::fstream out_file(filename_outs, std::fstream::out);
        assert(gens_file.is_open() && out_file.is_open());
        
        for(int i=0; i<=resolution_length; ++i) {
            int32_t M_rank; gens_file.read((char*)&M_rank, 4);
            HopfAlgebroid::cofree_comodule<motSt, MotDegree> F0;
            F0.load(gens_file);
            
            for(unsigned k=0; k<F0.generators.rank; ++k)
                out_file << "[" + std::to_string(i) + "-" + std::to_string(k) +"]\t:\t" << "(" << F0.generators.degree[k].deg - i << "," << F0.generators.degree[k].weight << ")" << std::endl;
        }
}

void motComplex::load(int resolution_length, string filename_generators, string filename_resolution) {
		terms.clear();
		maps.clear();
	
        std::fstream gens_file(filename_generators, std::fstream::in | std::fstream::binary);
        std::fstream res_file(filename_resolution, std::fstream::in | std::fstream::binary);
        assert(gens_file.is_open() && res_file.is_open());
	
		int32_t M_rank; gens_file.read((char*)&M_rank, 4);
		HopfAlgebroid::cofree_comodule<motSt, MotDegree> F0;
		F0.load(gens_file);
		
		terms.push_back(F0.generators);
		
		for(int iee=0; iee<resolution_length; ++iee){
			Matrices::matrix_array<tauPoly>  mp(1,&tau_module_oper, &tau_module_oper);
			mp.load(F0.generators.rank, res_file);
                        
   //                     std::cout << F0.generators.rank << std::flush;
			
			gens_file.read((char*)&M_rank, 4);
			F0.load(gens_file);
			
			maps.push_back(std::move(mp));
			terms.push_back(F0.generators);
		}
}

 string motComplex::output() {
            string result;
            for(auto &mp:maps)
                result += mp.output() + "\n";
            return result;
 }
 
 void motComplex::output(string filename_outs){
      std::fstream out_file(filename_outs, std::fstream::out);
      
      for(unsigned i=0; i<maps.size(); ++i)
          out_file << maps[i].output(std::to_string(i), std::to_string(i+1));
 }
 
int TauBockstein::dangerous_basket::leading_term(Modules::vectors<Matrices::matrix_index, MotSteentrod::tauPoly> const& v){
    int res =0;
    for(unsigned k=1;k<v.size();k++)
        if(tau_oper.tau_valuation(v.dataArray[res].coeficient)>tau_oper.tau_valuation(v.dataArray[k].coeficient))
            res=k;
    return res;
}

void TauBockstein::dangerous_basket::make_pretable(std::stack<int> &pot, Matrices::matrix<MotSteentrod::tauPoly> *M, tau_table &result, tau_table &pre_table){
    if(pot.empty()) return;
    int cur = pot.top();    pot.pop();
    auto ft = tau_module_oper.singleton(cur);
    auto fc = M->maps_to(ft);
 //   std::cout << "pp" << std::flush;
  //  std::cout << tau_module_oper.output(ft) + "\n";
  //  std::cout << tau_module_oper.output(fc) + "\n";
    while(!tau_module_oper.isZero(fc)){
 //   std::cout << "op";    
 //   std::cout << tau_module_oper.output(ft) + "\n";
 //   std::cout << tau_module_oper.output(fc) + "\n" << std::flush;
        auto ld = leading_term(fc);
        auto li = fc.dataArray[ld].ind;
        auto it = result.cycle_index.find(li);
    
        int val = tau_oper.tau_valuation(fc.dataArray[ld].coeficient);
        
        if(it==result.cycle_index.end()){
            result.insert(cur,li,ft,fc,val);
            make_pretable(pot,M,result,pre_table);
    //        std::cout << result.table.size();
            return;
        }
        
//std::cout << result.table[it->second].output(-1) << "\n";
    
        if(val<result.table[it->second].diff_length || (val==result.table[it->second].diff_length && cur > result.table[it->second].tag )){
            pot.push(result.table[it->second].tag);
            result.modify(it->second,cur,li,ft,fc,val);
            make_pretable(pot,M,result,pre_table);
            return;
        }
        
        auto mdf = tau_module_oper.singleton(result.table[it->second].tag, val - result.table[it->second].diff_length);
 //       mdf = tau_module_oper.minus(mdf);
        auto mfc = M->maps_to(mdf);
        
        ft = tau_module_oper.add(std::move(ft),std::move(mdf));
        fc = tau_module_oper.add(std::move(fc),std::move(mfc));
    }
    pre_table.insert(tau_table_entry::Invalid, cur,fc,ft, -1);
 //   std::cout << pre_table.output(-1);
    make_pretable(pot,M,result,pre_table);
    return;
}

void tau_table::insert(int tag, int cycle, int diff_length){
    static const Modules::vectors<Matrices::matrix_index, MotSteentrod::tauPoly> v0 = MotSteentrod::tau_module_oper.zero();
    int n = table.size();
    table.push_back({tag,cycle,v0,v0,diff_length});
    tag_index.emplace(tag,n);
    cycle_index.emplace(cycle,n);
}

void tau_table::insert(int tag, int cycle, Modules::vectors<Matrices::matrix_index, MotSteentrod::tauPoly> const &full_tag, Modules::vectors<Matrices::matrix_index, MotSteentrod::tauPoly> const &full_cycle, int diff_length){
    int n=table.size();
    table.push_back({tag,cycle,full_tag,full_cycle,diff_length});
    tag_index.emplace(tag,n);
    cycle_index.emplace(cycle,n);
}

void tau_table::modify(int n, int tag, int cycle, Modules::vectors<Matrices::matrix_index, MotSteentrod::tauPoly> const &full_tag, Modules::vectors<Matrices::matrix_index, MotSteentrod::tauPoly> const &full_cycle, int diff_length){
    cycle_index.erase(table[n].cycle);
    tag_index.erase(table[n].tag);
    
    table[n].tag = tag;
    table[n].cycle = cycle;
    table[n].diff_length = diff_length;
    table[n].full_cycle = full_cycle;
    table[n].full_tag = full_tag;
    
    tag_index.emplace(tag,n);
    cycle_index.emplace(cycle,n);
}

void tau_table::modify(int n, int tag, int cycle, int diff_length){
    cycle_index.erase(table[n].cycle);
    tag_index.erase(table[n].tag);
    
    table[n].tag = tag;
    table[n].cycle = cycle;
    table[n].diff_length = diff_length;
    
    tag_index.emplace(tag,n);
    cycle_index.emplace(cycle,n);
}

std::stack<int> TauBockstein::pot_maker(tau_table const& pre_table, int rank){
    std::set<int> psd;
    for(int i=0;i<rank;++i)
        psd.insert(i);
    for(auto k: pre_table.table)
        psd.erase(k.cycle);
    std::stack<int> result;
    for(auto k : psd) result.push(k);
    return result;
}

std::vector<tau_table> TauBockstein::make_table(motComplex& CMP){
    std::vector<tau_table> res_table;
    tau_table emp;
   res_table.resize(CMP.size(), emp);
    
    for(int i=1;i<CMP.size();++i)
    {
      //  std::cout << i << std::flush;
        auto pot=pot_maker(res_table[i-1],CMP.terms[i-1].rank);
     //   std::cout << pot.size() << "\n" << std::flush;
        dangerous_basket::make_pretable(pot, CMP.maps.data()+i-1, res_table[i], res_table[i-1]);
      //  std::cout << res_table[i].table.size();
  //      std::cout << res_table[i].output(i) << std::flush;
    //    std::cout << "table\n" << output_tables(res_table) << "??";
    }
    return res_table;
}

string tau_table_entry::output(int nm) const{
    if(tag==Invalid)
        return "{" + std::to_string(nm) + "-" + std::to_string(cycle) + "}\t:cycle";
    return "t^"+ std::to_string(diff_length) + "{" + std::to_string(nm) + "-" + std::to_string(cycle) + "}\t<-\t{" + std::to_string(nm-1)+"-"+std::to_string(tag)+"}";
}

string tau_table::output(int nm) const{
    string res;
    for(auto et: table)
        res += et.output(nm) + "\n";
    return res;
}

string TauBockstein::output_tables(const std::vector<tau_table>& Tb){
    string res;
    for(unsigned i=0;i<Tb.size();++i)
        res+=Tb[i].output(i);
    return res;
}

void cycle_data::get_cycles(tau_table table){
    for(auto &itm : table.table){
        if(itm.tag == tau_table_entry::Invalid)
            insert(itm.cycle,itm.full_cycle);
        else{
            std::function<tauP(const tauP&)> rl = [ &itm] (const tauP &r) {
                auto res = r;
                tau_oper.shift(res, -itm.diff_length);
                return res;
            };
            insert(itm.cycle,tau_module_oper.termwise_operation(rl,itm.full_cycle));
        }
    }
}

void cycle_data::get_tags(tau_table table){
    for(auto &itm:table.table){
        if(itm.tag != tau_table_entry::Invalid)
            insert(itm.tag,itm.full_tag);
    }
}

string cycle_data::output(int nm){
    string res;
    std::function<void(Modules::vectors<Matrices::matrix_index, MotSteentrod::tauPoly>&,int)> action = [nm,&res] (Modules::vectors<Matrices::matrix_index, MotSteentrod::tauPoly>& dt, int cl){
        res +=  "{" + std::to_string(nm) + "-" + std::to_string(cl) + "}" + "\t=\t" + tau_module_oper.output(dt,std::to_string(nm)) + "\n";
    };
    update_all(action);
    return res;
}

 void TauBockstein::make_cycle_tables(std::vector<tau_table> const& tab, std::vector<cycle_data>&result){
     result.resize(tab.size());
     for(unsigned i=0; i<tab.size(); ++i)
         result[i].get_cycles(tab[i]);
     for(unsigned i=1; i<tab.size(); ++i)
         result[i-1].get_tags(tab[i]);
 }

 new_vector TauBockstein::dangerous_basket::find_cycle(cycle_data& table, Modules::vectors<Matrices::matrix_index, MotSteentrod::tauPoly> v){
     new_vector res;
     while(!tau_module_oper.isZero(v)){
         auto ld = leading_term(v);
         auto cf = v.dataArray[ld].coeficient;
         auto ind = v.dataArray[ld].ind;
         
         v = tau_module_oper.add(v,tau_module_oper.scalor_mult(cf, table.search(ind)));
         res = tau_module_oper.add(res, tau_module_oper.singleton(ind, cf));
   //      std::cout << tau_module_oper.output(v) << "\n" << std::flush;
     }
     return res;
 }
 
 string TauBockstein::dangerous_basket::output_cycles(int nm, const new_vector&v){
     string res = "o";
     for(auto tm : v.dataArray){
         res += "+t^" + std::to_string(tm.coeficient) + "{" + std::to_string(nm) + "-" + std::to_string(tm.ind) + "}";
     }
     return res;
 }
 
 void TauBockstein::multiplication_table(motSt const& x, int x_deg, string filename_generators, string filename_res_truns, string filename_outs, int resolution_length, MSOper& ms_oper, std::vector<cycle_data>& cyc_table){
        Matrices::matrix_array<MotSteentrod::tauPoly>  mult(1,&tau_module_oper, &tau_module_oper);
        ms_oper.make_multiplication_table(x, x_deg, &mult);
        
        std::fstream resmaps_file(filename_res_truns, std::fstream::in | std::fstream::binary);
        std::fstream gens_file(filename_generators, std::fstream::in | std::fstream::binary);
        std::fstream out_file(filename_outs, std::fstream::out);
        assert(resmaps_file.is_open() && gens_file.is_open() && out_file.is_open());
		
        int32_t M_rank; gens_file.read((char*)&M_rank, 4);
        HopfAlgebroid::cofree_comodule<motSt, MotDegree> F0;
        F0.load(gens_file);
        
        for(int iee=0; iee<resolution_length-1; ++iee) {
            
        Matrices::matrix_array<MotSteentrod::tauPoly>  mp(1,&tau_module_oper, &tau_module_oper);
        Matrices::matrix_array<MotSteentrod::tauPoly> mmt(1,&tau_module_oper,&tau_module_oper);
 //       std::cout << F0.rank() << std::flush;
        mp.load(F0.rank(),resmaps_file);
           
        for(unsigned i=0; i<F0.generators.rank; ++i) 
                if((int)F0.generators.degree[i].deg + x_deg <= (int)ms_oper.maxDeg) {
                    auto des = mp.maps_to(F0.multiply_using_table(F0.position_of_gens[i],mult));
                    mmt.insert(i,des);
   //                 out_file << "[" + std::to_string(iee) + "-" + std::to_string(i) + "]\t->\t" + tau_module_oper.output(des, std::to_string(iee+1)) + "\n";
           }
           
           for(unsigned i=0; i<F0.generators.rank; ++i) 
                if((int)F0.generators.degree[i].deg + x_deg <= (int)ms_oper.maxDeg) {
                    auto des = mmt.maps_to(cyc_table[iee].search(i));
                    auto mres = dangerous_basket::find_cycle(cyc_table[iee+1],des);
                   out_file << "{" + std::to_string(iee) + "-" + std::to_string(i) + "}\t->\t" + dangerous_basket::output_cycles(iee+1, mres) + "\n";
           }
           
            gens_file.read((char*)&M_rank, 4);
            F0.load(gens_file);
        }
} 