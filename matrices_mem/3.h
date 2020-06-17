//Curtis tables stored in memories
template<typename ring>
class curtis_table_mem : public curtis_table<ring>{
	//the index of the entries
	std::map<matrix_index,typename curtis_table<ring>::entry> data;
public:
	//clear the table
	void clear(){
		data.clear(); 
//		last = -1;
	}
	
	//check is an entry exists
	bool is_member(matrix_index n){
		return data.count(n)==1; }
		
	//find an entry
	typename curtis_table<ring>::entry search(matrix_index n){
		return data.at(n); }
		
	//find the reference to entry
	typename curtis_table<ring>::entry* search_ref(matrix_index n){
		return &(data.find(n)->second); }
		
	//modify (or act on entries) all the entries with a function
	void modify_all(std::function<void(typename curtis_table<ring>::entry&)> action){
		for(auto &tm: data)
			action(tm.second);
	}
	
	//modify (or act on entries) all the entries with a function
	void run_through(std::function<void(const typename curtis_table<ring>::entry&)> action){
		for(auto &tm: data)
			action(tm.second);
	}
	
	void modify_all(std::function<bool(std::function<typename curtis_table<ring>::entry()> original_entry, typename curtis_table<ring>::entry& returned_entry)> action){
		for(auto &tm: data){
			std::function<typename curtis_table<ring>::entry()> original_entry = [&tm]{
				return tm.second; };
			typename curtis_table<ring>::entry new_entry;
			if(action(original_entry,new_entry))
				tm.second = new_entry;
		}
	}
	
	//insert a new entry
	void insert(matrix_index cyc, matrix_index tg, vectors<matrix_index,ring> const &f_cyc, vectors<matrix_index,ring> const &f_tg){
		typename curtis_table<ring>::entry new_entry = {cyc, tg, f_cyc, f_tg};
		data[cyc] = new_entry;
//		if(cyc > last) last = cyc;
	}
	
	//the number of entries
	int num_entries(){
		return data.size();
	}
};
