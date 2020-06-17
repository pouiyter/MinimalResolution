     virtual vect maps_to(vect_ref x, std::function<vect(index)> rows){
         std::function<vect(int)> summands = [this, &x, rows] (int i){
             return scalor_mult(x.dataArray[i].coeficient, rows(x.dataArray[i].ind));};
        std::function<vect(vect&&,vect&&)> adding = [this] (vect&&x, vect&&y) { return add(std::move(x),std::move(y));};
        return Paralell::sum(summands,0,x.size(),adding);
     }

         template<typename modul1,typename modul2>
    static modul2 maps_to0(vect const &y, std::function<modul2(const R&, const modul1&)> scalor_mult, Algebra::AbGroupOp<modul2> *modop2,
		  unsigned start, unsigned stop, std::function<modul1(index)> find) {
     //   std::cout << start << " " << stop << "\n" << std::flush;
        auto x = y.summon();
      if(stop<=start) return modop2->zero();
      if(stop-start==1) return scalor_mult(x.dataArray[start].coeficient, find(x.dataArray[start].ind));
      unsigned mid = (start+stop) / 2; 
      return modop2->add(maps_to0(x,scalor_mult,modop2,start,mid,find), maps_to0(x,scalor_mult,modop2,mid,stop,find));
    }
    
        template<typename modul1,typename modul2>
    static modul2 maps_to_para(vect const &y, std::function<modul2(const R&, const modul1&)> scalor_mult, Algebra::AbGroupOp<modul2> *modop2,
        std::function<modul1(index)> find
    ) {
        auto x = y.summon();
      std::function<modul2(int)> enties = [scalor_mult, &x,find] (int i) {
        return scalor_mult(x.dataArray[i].coeficient, find(x.dataArray[i].ind)); };
      std::function<modul2(modul2&&,modul2&&)> adding = [modop2] (modul2 &&x, modul2 &&y) { return modop2->add(std::move(x), std::move(y)); };
      if(x.size()==0) return modop2->zero();
      return Paralell::sum(enties, 0, x.size(), adding);
    }
    
     template<typename index1>
    index1 symplify_to_led(vect &x, std::function<vect(index1)> enties_cycle, const std::map<index,index1> &cycle_index){
     //     auto &entries = datas.data;
          std::function<vect(vect&&,vect&&)> adder = [this] (vect &&x, vect &&y) { return this->add(std::move(x),std::move(y)); };
          std::function<vect(int)> smd = [this, &x, &enties_cycle, &cycle_index] (int i) { 
     //         ring c = this->ModOper->ringOper->minus(this->ModOper->component(entries[i].cycle,x));
              auto si = cycle_index.find(x.dataArray[i].ind);
              if(si==cycle_index.end())
                  return this->zero();
              auto c = this->ringOper->minus(x.dataArray[i].coeficient);
              return this->scalor_mult(c,enties_cycle(si->second)); 
          };
          if(x.size()!=0){
            auto sbt = Paralell::sum(smd,0,x.size(),adder);
        //    std::cout << this->ModOper->output(sbt) << "sss\n";
        //    std::cout << this->ModOper->output(x) << "se\n";
            x = this->add(std::move(x),std::move(sbt));
        //    std::cout << this->ModOper->output(x) << "sdf\n";
          }
          auto nst = this->first_invertible(x);
          if(nst>=x.size())
              return BOUNDARY_IND;
          else
              return x.dataArray[nst].ind;
      }

      
      template<typename modul1,typename modul2>
      static modul2 maps_to(vect const &x, std::function<modul2(const R&, const modul1&)> scalor_mult,
							Algebra::AbGroupOp<modul2> *modop2,
							std::function<modul1(index)> find) {
		  return maps_to0(x,scalor_mult,modop2,0,x.size(), find);
							}
							
							
							
							//scaling to make the term with index n have unit coeficient, if that is invertible
							virtual vectors<index,R> unify(index n, vectors<index,R> const &x){
								R po = ringOper->inverse(component(n,x));
								return scalor_mult(po,x);
							}
							
							//subtract a vect y from x to kill the term with index n. y should be unified beforehand
							virtual void take_away(index n, vectors<index,R> &x, vectors<index,R> const &y){
								R cp = ringOper->minus(component(n,x));
								vectors<index,R> smd = scalor_mult(cp,y);
								x = add(std::move(x),std::move(smd));
							}
///


template<typename R> 
class matrix_array_data : public MapTable::array_data<matrix_index, typename matrix<R>::vect> {
public:
	matrix_array_data(){}
	matrix_array_data(unsigned n) : MapTable::array_data<matrix_index, typename matrix<R>::vect>(n) {}
	unsigned source_rank() { return this->data.size(); }
};

template<typename R>
class matrix_array : virtual public matrix<R> {
	matrix_array_data<R> data;
public:
	matrix_array(unsigned rank, Modules::ModuleOp<matrix_index,R> *moop, Algebra::IOSet<typename matrix<R>::vect> *ioop) : 
	matrix<R>(moop, ioop, &this->data), data(rank)    {}
	
	matrix_array(unsigned rank, Modules::ModuleOp<matrix_index,R> *moop) : 
	matrix<R>(moop, moop, &this->data), data(rank)    {}
	
	matrix_array(stR, unsigned rank, Modules::ModuleOp<matrix_index,R> *moop) : 
	matrix<R>(moop, moop, &this->data), data(rank)    {}
	
	matrix_array(Modules::ModuleOp<matrix_index,R> *moop) : 
	matrix<R>(moop, moop, &this->data)    {}
	
	matrix_array(matrix_array<R> &&x) : matrix_array(0,x.moduleOper) {
		data = std::move(x.data);
	}
	
	int rank(){ return data.data.size(); }
	
	// row transfermations	
	void pre_gausian(matrix_index row, matrix_index col){
		//        std::cout << row << "," << col << std::flush;
		auto uv = this->moduleOper->unify(col,data.data[row]);
		//     std::cout << "unified" << std::flush;
		std::function<void(int)> action = [this, &uv, col, row] (int i) {
			if(i==row) this->data.data[i] = uv;
			else this->moduleOper->take_away(col, this->data.data[i], uv); };
			Paralell::para_for(0, (int) data.data.size(),action);
	}
	
	void gaussian(std::vector<std::pair<matrix_index,matrix_index>> const & row_cols){
		int counter = 0;
		for(auto rc: row_cols){
			#if defined TESTING
			if(Testing::interactiveMode)
				std::cout << "\r" << counter++ << "/" << row_cols.size() << std::flush;
			#endif
			pre_gausian(rc.first,rc.second);
		}
	}
	
	void double_gaussian(std::vector<std::pair<matrix_index,matrix_index>> const & row_cols, matrix<R> *B, unsigned shifting){
		std::function<void(int)> summing = [this,B,shifting] (int i) {
			this->data.data[i].direct_sum(B->find(i),shifting); };
			Paralell::para_for(0,(int)this->data.data.size(),summing);
			
			int counter = 0;
			for(auto rc: row_cols){
				std::cout << "\r" << counter++ << "/" << row_cols.size() << std::flush;
				pre_gausian(rc.first,rc.second);
			}
			
			std::function<Modules::vectors<matrix_index,R>(int)> un_suuming = [this,shifting] (int i) {
				return this->data.data[i].un_direct_sum(shifting); };
				B->construct(this->data.data.size(),un_suuming);
	}
	
	void logged_gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols, matrix<R> *logs){
		logs->set2unit(rank());
		unsigned shifting = (1<<30)-rank();
		double_gaussian(row_cols, logs, shifting);
	}
	
};

template<typename R>
class curtis_table{
public:
	Modules::ModuleOp<matrix_index,R> *ModOper;
public:
	typedef typename matrix<R>::vect vect;
	typedef struct { 
		matrix_index cycle, tag;
		vect full_cycle, full_tag;
	} entry;
	MapTable::data_base<matrix_index, entry> *data;
	
	curtis_table(Modules::ModuleOp<matrix_index,R> *moduleOper, MapTable::data_base<matrix_index,entry> *datB){
		ModOper = moduleOper;	data = datB;
	}
	
	void clear() { data->clear(); }
	
	static constexpr matrix_index Boundary = BOUNDARY_IND;
private:
	matrix_index symplify_to_led0(vect &x,matrix_index start){
		if(ModOper->isZero(x)) return Boundary;
		auto nst = ModOper->next_non_trvial_entry(start,x);
		if(nst>ModOper->last_entry(x))
			return Boundary;
		if(nst>start)
			return symplify_to_led(x,nst);
		if(data->is_member(start)){
			ModOper->take_away(start, x, data->search(start).full_cycle);
			return symplify_to_led(x,start+1);
		}
		if(ModOper->ROper->invertible(ModOper->component(start,x))) return start;
		else return symplify_to_led(x,start+1);
	}
	
	matrix_index symplify_to_led(vect &x,matrix_index start, vect &homotopy){
		//    std::cout << "start:" << start << ModOper->output(x)  <<"\n";
		if(ModOper->isZero(x)) return Boundary;
		auto nst = ModOper->next_non_trvial_entry(start,x);
		if(nst>ModOper->last_entry(x))
			return Boundary;
		if(nst>start)
			return symplify_to_led(x,nst,homotopy);
		if(data->is_member(start)){
			R cp = ModOper->ringOper->minus(ModOper->component(start,x));
			x = ModOper->add(std::move(x), std::move(ModOper->scalor_mult(cp,data->search(start).full_cycle)));
			homotopy = ModOper->add(std::move(homotopy), std::move(ModOper->scalor_mult(cp,data->search(start).full_tag)));
			return symplify_to_led(x,start+1,homotopy);
		}
		if(ModOper->ringOper->invertible(ModOper->component(start,x))) return start;
		else return symplify_to_led(x,start+1,homotopy);
	}
public:
	matrix_index symplify_to_led(vect &x) { return symplify_to_led0(x,0); }
	matrix_index symplify_to_led(vect &x, vect &homotopy) { return symplify_to_led(x,0,homotopy); }
	
	void modify_all(std::function<void(entry&)> action) {
		data->update_all(action);
	}
	
	void insert(matrix_index cyc, matrix_index tg, vect const &f_cyc, vect const &f_tg) {
		entry new_entry = {cyc, tg, f_cyc, f_tg};
		data->insert(cyc,new_entry);
	}
	
	std::vector<matrix_index> cycle_matrix(matrix<R> &result, unsigned source_rank) {
		std::vector<matrix_index> leading_term;
		leading_term.resize(source_rank);
		result.clear();
		std::function<void(entry&)> action = [this, &leading_term, &result] (entry &itm) {
			leading_term[itm.tag] = itm.cycle;
			result.insert(itm.tag, itm.full_cycle); };
			data->update_all(action);
			return leading_term;
	}
	
	
	std::vector<matrix_index> cycle_matrix(matrix<R> &cyc_mat, matrix<R> &tag_mat, unsigned source_rank) {
		std::vector<matrix_index> leading_term;
		leading_term.resize(source_rank);
		cyc_mat.clear();
		tag_mat.clear();
		std::function<void(entry&)> action = [this, &leading_term, &cyc_mat, &tag_mat] (entry &itm) {
			leading_term[itm.tag] = itm.cycle;
			cyc_mat.insert(itm.tag, itm.full_cycle);
			tag_mat.insert(itm.tag, itm.full_tag);
		};
		data->update_all(action);
		return leading_term;
	}
	
	std::vector<matrix_index> cycle_matrix(unsigned source_rank) {
		std::vector<matrix_index> leading_term;
		leading_term.resize(source_rank);
		std::function<void(entry&)> action = [this, &leading_term] (entry &itm) {
			leading_term[itm.tag] = itm.cycle;
		};
		data->update_all(action);
		return leading_term;
	}
	
	string output() {
		string res = "";
		std::function<void(entry&)> ous = [&res] (entry& et) { res += std::to_string(et.cycle) + "\t<-\t" + std::to_string(et.tag) + "\n"; };
		data->update_all(ous);
		return res;
	}
	
	template<typename R2>
	std::vector<matrix_index> cycle_matrix(matrix<R2> &result, unsigned source_rank, std::function<Modules::vectors<matrix_index,R2>(const vect&)> transformer, matrix<R2> *Map) {
		std::vector<matrix_index> leading_term;
		leading_term.resize(source_rank);
		//    result.clear();
		//    std::cout << "f";
		//  std::cout << output();
		std::function<void(entry&)> action = [this, &leading_term, &result, &transformer,Map] (entry &itm) {
			leading_term[itm.tag] = itm.cycle;
			auto new_tag = transformer(itm.full_tag);
			auto new_cyc = Map->maps_to(new_tag);
			result.insert(itm.tag, new_cyc);
			//     std::cout << itm.tag << "d";
		};
		data->update_all(action);
		return leading_term;
	}
};

template<typename R>
class curtis_table_maps : virtual public curtis_table<R> {
	MapTable::maps<matrix_index, typename curtis_table<R>::entry> datas;
public:
	curtis_table_maps(Modules::ModuleOp<matrix_index,R> *moduleOper) : curtis_table<R>(moduleOper, &this->datas) {}
};

template<typename R>
class quasi_table : public curtis_table<R> {
	MapTable::array_data<matrix_index,typename curtis_table<R>::entry> datas;
	std::map<matrix_index,matrix_index> cycle_index;
	typedef typename matrix<R>::vect vect;
	//   std::vector<typename curtis_table<R>::entry> entries;
public:
	quasi_table(Modules::ModuleOp<matrix_index,R> *moduleOper) : curtis_table<R>(moduleOper, &this->datas) {}
	
	matrix_index symplify_to_led(vect &x){
		std::function<vect(matrix_index)> enties_cyc = [this] (matrix_index i) {
			return this->datas.data[i].full_cycle;
		};
		return this->ModOper->symplify_to_led(x,enties_cyc,cycle_index);
	}
	
	matrix_index symplify_to_led0(vect &x){
		auto &entries = datas.data;
		std::function<vect(vect&&,vect&&)> adder = [this] (vect &&x, vect &&y) { return this->ModOper->add(std::move(x),std::move(y)); };
		std::function<vect(int)> smd = [this, &x, &entries] (int i) { 
			//         R c = this->ModOper->ringOper->minus(this->ModOper->component(entries[i].cycle,x));
			auto si = cycle_index.find(x.dataArray[i].ind);
			if(si==cycle_index.end())
				return this->ModOper->zero();
			R c = this->ModOper->ringOper->minus(x.dataArray[i].coeficient);
			return this->ModOper->scalor_mult(c,entries[si->second].full_cycle); 
		};
		if(x.size()!=0){
			auto sbt = Paralell::sum(smd,0,x.size(),adder);
			//    std::cout << this->ModOper->output(sbt) << "sss\n";
			//    std::cout << this->ModOper->output(x) << "se\n";
			x = this->ModOper->add(std::move(x),std::move(sbt));
			//    std::cout << this->ModOper->output(x) << "sdf\n";
		}
		auto nst = this->ModOper->first_invertible(x);
		if(nst>=x.size())
			return curtis_table<R>::Boundary;
		else
			return x.dataArray[nst].ind;
	}
	
	void insert(matrix_index cyc, matrix_index tg, vect const &f_cyc){
		static vect nul = this->ModOper->zero();
		auto const nf_cyc = this->ModOper->unify(cyc,f_cyc);
		#pragma omp parallel for
		for(unsigned i=0;i<datas.data.size();++i){
			this->ModOper->take_away(cyc,datas.data[i].full_cycle,nf_cyc);
		}
		cycle_index.emplace(cyc,datas.data.size());
		typename curtis_table<R>::entry new_entry = {cyc,tg,nf_cyc,nul};
		datas.data.push_back(new_entry);
	}
};

}

//change rational notation to integral notation
BP BP_Op::Q2Z(const BPQ& x, BPQ_Op& BPQoper) {
	std::function<Zp(const Qp&)> ct = [&BPQoper](const Qp &r){ 
		return BPQoper.Qp_oper->int_part(r); };
		return BPQoper.termwise_operation(ct,x);
}

//change rational notation to integral notation
BPBP BP_Op::Q2Z(const BPBPQ& x, BPQ_Op& BPQoper){
	std::function<BP(const BPQ&)> ct = [this,&BPQoper](const BPQ &r){ 
		return Q2Z(r,BPQoper); };
		return BPQoper.BPBPQ_oper.termwise_operation(ct,x);
}

//change rational notation to integral notation
BPBPBP BP_Op::Q2Z(const BPBPBPQ& x, BPQ_Op& BPQoper){
	std::function<BPBP(const BPBPQ&)> ct = [this,&BPQoper](const BPBPQ &r){ 
		return Q2Z(r,BPQoper); };
		return BPQoper.BPBPBPQ_oper.termwise_operation(ct,x);
}

//make the structure tables for the generators
void make_tables_gens(int maxVar, string R2Lfilename, string deltafilename, BPQ_Op BPQoper){
	fstream R2L_file(R2Lfilename, std::ios::out | std::ios::binary);
	for(int i=1; i<=maxVar; ++i){
		
	}
}

// #define DECLARE_ARRAY_MATRIX(ring,rank,module_oper,name) Matrices::matrix_array_data<ring> name##000(rank); Matrices::matrix<ring> name(module_oper, module_oper, &(name##000));

#include"filedmatrices.h"


//direct sum of two matrices
void direct_sum(matrix<R> const &Y, unsigned target_rank);
//direct sum of two matrices
template<typename R>
void matrix<R>::direct_sum(matrix<R> const &Y, unsigned tar_rank){
	
	std::function<void(vect&,matrix_index)> row_summer = [&Y,tar_rank] (vect& ro, matrix_index i) { ro.direct_sum(Y.find(i), tar_rank); };
	this->dataBase->update_all(row_summer);
}


//delete the coloms with index in a set
void del_cols(std::set<matrix_index> to_del){
	std::function<bool(matrix_index)> tule = [&to_del] (matrix_index n) { return to_del.count(n) == 0; };
	filter(tule);
}
