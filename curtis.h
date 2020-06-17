//curtis.h
//make Curtis tables
#pragma once
#include"matrices.h"
#include"matrices_mem.h"
#include"Fp.h"
#include"streams.h"

//curtis table stored in streams
template<typename ring>
class curtisTable_stream : public curtis_table<ring>{
	//the index of the entries
	std::map<matrix_index,std::ios::streampos> data_index;
	//the cached table
	curtis_table_mem<ring> cached;
	//stream of data
	con_streams stream_data;
public:
	//constructor
	curtisTable_stream(std::iostream* datas) : stream_data(datas){}
	
	//clear the table
	void clear(){
		data_index.clear(); 
		cached.clear();
		stream_data.fclear();
	}
	
	//check is an entry exists
	bool is_member(matrix_index n){
		if(cached.is_member(n))
			return true;
		return data_index.count(n)==1;
	}
	
	//find an entry
	typename curtis_table<ring>::entry search(matrix_index n){
		if(cached.is_member(n))
			return cached.search(n);
		std::ios::streampos pos = data_index[n];
		std::function<typename curtis_table<ring>::entry(std::iostream&)> reader = [this](std::iostream &ins){
			return this->load_entry(ins); };
		return stream_data.read(reader, pos);
	}
	
	typename curtis_table<ring>::entry* search_ref(matrix_index n){
		return cached.search_ref(n); }
	
	//insert a new entry
	void insert(matrix_index cyc, matrix_index tg, vectors<matrix_index,ring> const &f_cyc, vectors<matrix_index,ring> const &f_tg){
		cached.insert(cyc, tg, f_cyc, f_tg);
// 		if(cyc > last) last = cyc;
	}
	
	//flush the cached data
	void flush(){
		std::function<void(typename curtis_table<ring>::entry&)> action = [this](typename curtis_table<ring>::entry& et){
			std::function<void(std::iostream&)> writer = [&et,this](std::iostream &outs){
				this->save_entry(outs, et); };
			auto pos = stream_data.write(writer);
			data_index[et.cycle] = pos;
		};
		cached.modify_all(action);
		cached.clear();
	}
	
	//the number of entries
	int num_entries(){
		flush();
		return data_index.size();
	}
	
	//modify (or act on entries) all the entries with a function
	void run_through(std::function<void(const typename curtis_table<ring>::entry&)> action){
		flush();
		for(auto tm: data_index){
			std::ios::streampos pos = tm.second;
			std::function<typename curtis_table<ring>::entry(std::iostream&)> reader = [this](std::iostream &ins){
				return this->load_entry(ins); };
			auto et = stream_data.read(reader, pos);
			action(et);
		}
	}
};

