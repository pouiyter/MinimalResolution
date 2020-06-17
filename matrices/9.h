//curtis table
template<typename ring>
class curtis_table{
public:
	//module operations
	static ModuleOp<matrix_index,ring> *ModOper;
	//the entries of a curtis table
	typedef struct{ 
		matrix_index cycle, tag;
		vectors<matrix_index,ring> full_cycle, full_tag;
	} entry;
	
	//the last entry
//	matrix_index last;
	
	//clear the table
	virtual void clear()=0;
	//check if an element is in the table
	virtual bool is_member(matrix_index)=0;
	//find the entry
	virtual entry search(matrix_index)=0;
	//find the reference of entry
	virtual entry* search_ref(matrix_index)=0;
	//the number of entries
	virtual int num_entries()=0;
	
	//flush the cached entries
	virtual void flush(){}
	
	//the index indicating for a boundary
	static constexpr matrix_index Boundary = -1;
	
	//simplify a vector using the table, preserving the homotopy
	matrix_index symplify_to_led(std::vector<ring> &x, matrix_index start, std::vector<ring> *homotopy=NULL){
		//zero is a boundary
		if(start >= x.size()) return Boundary;
		
		//find the current leading term
 		if(ModOper->ringOper->isZero(x[start]))
			return symplify_to_led(x, start+1, homotopy);
		
		//when the leading term is in the table
		if(is_member(start)){
			//fint the leading coeficient
			ring cp = ModOper->ringOper->minus(x[start]);
			//find the entry
			entry *etr;
			etr = search_ref(start);
		//	entry ks = search(start);
		//	etr = &ks;
			//take away the term in the table from x
			etr->full_cycle.add2Dense(cp, x, ModOper->ringOper);
			
			if(homotopy!=NULL)
				etr->full_tag.add2Dense(cp, *homotopy, ModOper->ringOper);

			//take care of the next leading term
			return symplify_to_led(x,start+1,homotopy);
		}
		
		//the leading term is invertible, gives a non trivial cycle
		if(ModOper->ringOper->invertible(x[start])) return start;
		else return symplify_to_led(x,start+1,homotopy);
	}
	
private:
	//simplify a vector using the table, preserving the homotopy
	matrix_index symplify_to_led(vectors<matrix_index,ring> &x,matrix_index start, vectors<matrix_index,ring> *homotopy=NULL){
		//zero is a boundary
		if(ModOper->isZero(x)) return Boundary;
		
		//find the current leading term
		auto nst = ModOper->next_non_trvial_entry(start,x);
		
		//all the entries are nil
		if(nst>ModOper->last_entry(x))
			return Boundary;
		
		//simplify the lleading term
		if(nst>start)
			return symplify_to_led(x,nst,homotopy);
		//when the leading term is in the table
		if(is_member(start)){
			//fint the leading coeficient
			ring cp = ModOper->ringOper->minus(ModOper->component(start,x));
			//take away the term in the table from x
			x = ModOper->add(std::move(x), std::move(ModOper->scalor_mult(cp,search(start).full_cycle)));
			if(homotopy!=NULL)
				*homotopy = ModOper->add(std::move(*homotopy), std::move(ModOper->scalor_mult(cp,search(start).full_tag)));
			//take care of the next leading term
			return symplify_to_led(x,start+1,homotopy);
		}
		
		//the leading term is invertible, gives a non trivial cycle
		if(ModOper->ringOper->invertible(ModOper->component(start,x))) return start;
		else return symplify_to_led(x,start+1,homotopy);
	}
	
public:
	//simplify a vector
	virtual matrix_index symplify_to_led(vectors<matrix_index,ring> &x) { 
		return symplify_to_led(x,0); }
		
	virtual matrix_index symplify_to_led(vectors<matrix_index,ring> &x, vectors<matrix_index,ring> &homotopy) { 
		return symplify_to_led(x,0,&homotopy); }
	
	matrix_index simplify_to_led(int target_rank, vectors<matrix_index,ring> &x){
		//zero is a boundary
		if(ModOper->isZero(x)) return Boundary;
		//convert to dense vector
		matrix_index start = x.dataArray[0].ind;
		std::vector<ring> xx = x.toDense(target_rank);
		//simplify
		auto res = symplify_to_led(xx, start, NULL);
		//change back
		x.deDense(xx,ModOper->ringOper,start);
		return res;
	}
	
	matrix_index simplify_to_led(int target_rank, int source_rank, vectors<matrix_index,ring> &x, vectors<matrix_index,ring> &homotopy){
		//zero is a boundary
		if(ModOper->isZero(x)) return Boundary;
		//convert to dense vector
		matrix_index start = x.dataArray[0].ind;
		std::vector<ring> xx = x.toDense(target_rank);
		std::vector<ring> htpy = homotopy.toDense(source_rank);
		//simplify
		std::cout << "t" << std::flush;
		auto res = symplify_to_led(xx, start, &htpy);
		//change back
		std::cout << "b" << std::flush;
		x.deDense(xx,ModOper->ringOper,start);
		homotopy.deDense(htpy,ModOper->ringOper,0);
//		return symplify_to_led(x,homotopy);
		return res;
	}
	
	//modify all the entries with a function
	virtual void modify_all(std::function<void(entry&)> action){
		std::cerr << "not supported" << std::flush;
	}
	
	//modify all the entries
	virtual void modify_all(std::function<bool(std::function<entry()> original_entry, entry& returned_entry)> action){
		std::cerr << "not supported" << std::flush;
	}
	
	//run through the entries
	virtual void run_through(std::function<void(const entry&)>)=0;
	
	//insert a new entry
	virtual void insert(matrix_index cycle, matrix_index tag, vectors<matrix_index,ring> const &full_cyc, vectors<matrix_index,ring> const &full_tag)=0;
	
	//construct a matrix using the entries of the table
	std::vector<matrix_index> cycle_matrix(matrix<ring> *cyc_mat, matrix<ring> *tag_mat, unsigned source_rank){
		std::vector<matrix_index> leading_term;
		leading_term.resize(source_rank);
		if(cyc_mat!=NULL)
			cyc_mat->clear();
		if(tag_mat!=NULL)
			tag_mat->clear();
		std::function<void(const entry&)> action = [this, &leading_term, &cyc_mat, &tag_mat] (const entry &itm){
			//leading terms of the entries
			leading_term[itm.tag] = itm.cycle;
			//matrix formed by the full cycles
			if(cyc_mat!=NULL)
				cyc_mat->insert(itm.tag, itm.full_cycle);
			//matrix formed by the full tags
			if(tag_mat!=NULL)
				tag_mat->insert(itm.tag, itm.full_tag);
		};
		run_through(action);
		return leading_term;
	}
	
	//construct a matrix using the entries of the table
	std::vector<matrix_index> cycle_matrix(matrix<ring> &result, unsigned source_rank){
		return cycle_matrix(&result, NULL, source_rank);
	}
	
	std::vector<matrix_index> cycle_matrix(unsigned source_rank){
		return cycle_matrix(NULL,NULL,source_rank);
	}
	
	//transforming the cycles
	template<typename ring2>
	std::vector<matrix_index> cycle_matrix(matrix<ring2> &result, unsigned source_rank, std::function<vectors<matrix_index,ring2>(const vectors<matrix_index,ring>&)> transformer, matrix<ring2> *Map) {
		std::vector<matrix_index> leading_term;
		leading_term.resize(source_rank);
		std::function<void(entry&)> action = [this, &leading_term, &result, &transformer,Map] (entry &itm) {
			leading_term[itm.tag] = itm.cycle;
			//transform the full tag
			auto new_tag = transformer(itm.full_tag);
			//the new cycle is the image of the new tag under the Map
			auto new_cyc = Map->maps_to(new_tag);
			result.insert(itm.tag, new_cyc);
		};
		modify_all(action);
		return leading_term;
	}
	
	//save entries
	void save_entry(std::iostream &writer, const entry &et){
		writer.write((char*)&et.cycle, 4);
		writer.write((char*)&et.tag, 4);
		ModOper->save(et.full_cycle, writer);
		ModOper->save(et.full_tag, writer);
	}
	
	//save the table, 4 bytes for the size, followed by entries
	void save(std::iostream &writer){
		int32_t sz = num_entries();
		writer.write((char*)&sz, 4);
		//4 bytes for cycle, 4 bytes for tag, then the full cycle, the full tag
		std::function<void(const entry&)> action = [&writer, this](const entry &et){
			save_entry(writer, et);
		};
		run_through(action);
	}
	
	//load the entry
	entry load_entry(std::iostream &reader){
		entry et;
		reader.read((char*)&et.cycle, 4);
		reader.read((char*)&et.tag, 4);
		et.full_cycle = ModOper->load(reader);
		et.full_tag = ModOper->load(reader);
		return et;
	}
	
	//load the table
	void load(std::iostream &reader){
		clear();
		
		int32_t sz;
		reader.read((char*)&sz, 4);
		for(int i=0; i<sz; ++i){
			entry et = load_entry(reader);
			insert(et.cycle, et.tag, et.full_cycle, et.full_tag);
			flush();
		}
	}
	
	//output the table
	string output(){
		string res;
		std::function<void(entry&)> action = [this, &res](entry &et){
			res += "cycle: " + std::to_string(et.cycle) + ": " + ModOper->output(et.full_cycle) + "\n";
			res += "tag: " + std::to_string(et.tag) + ": " + ModOper->output(et.full_tag) + "\n";
		};
		modify_all(action);
		return res;
	}
};

//the module operations for the curtis tables
template<typename ring>
ModuleOp<matrix_index,ring> *curtis_table<ring>::ModOper;
