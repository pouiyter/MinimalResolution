//matrices_stream.h
#pragma once
#include"matrices.h"
#include"matrices_mem.h"
#include"streams.h"

//matrix stored in streams
template<typename ring>
class matrix_stream : public matrix<ring>{
	//the data for a matrix, consisting of a stream and an array of positions
	con_streams *datas;
	std::vector<std::ios::streampos> datapos;
public:
	//the constructor
	matrix_stream(con_streams *ds){
		datas = ds;
		this->rank = 0;
	}
	
	//clear the contents
	void clear(){
		datas->fclear();
		datapos.clear();
		this->rank = 0;
	}
	
	//find the n-th row
	vectors<matrix_index,ring> find(matrix_index n) const{
		auto pos = datapos[n];
		std::function<vectors<matrix_index,ring>(std::iostream&)> reader = [this](std::iostream &is){
			return this->moduleOper->load(is); };
		return datas->read(reader, pos);
	}
	
	//insert a new row
	void insert(matrix_index i, vectors<matrix_index,ring> const& x){
		if(i>=datapos.size()){
			datapos.resize(i+1);
			this->rank = i+1;
		}
		std::function<void(std::iostream&)> writer = [this,&x](std::iostream &os){
			this->moduleOper->save(x, os); };
		datapos[i] = datas->write(writer);
	}
	
	//set the rank, note that the state is invalid after this operation
	void set_rank(unsigned n){
		this->rank = n;
		datapos.resize(n);
	}
	
	//update all rows
	void update_all(std::function<void(vectors<matrix_index,ring>&,matrix_index)> action){
		std::cerr << "action not supported!";
	}
	
	//Gaussian ellimination, with given rows and colums
	void gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols){
		matrix_mem<ring> M;
		M.construct(this);
		M.gaussian(row_cols);
		this->construct(&M);
	}
	
	//delete some columns and then do Gaussion
	void del_and_gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols, std::set<int> const &to_del){
		matrix_mem<ring> M;
		M.construct(this);
		M.del_and_gaussian(row_cols,to_del);
		this->construct(&M);
	}
};

//matrix stored in files
template<typename ring>
class matrix_file : public matrix_stream<ring>{
	con_fstreams files;
public:
	matrix_file(string filename) : matrix_stream<ring>(&files), files(filename){}
};
