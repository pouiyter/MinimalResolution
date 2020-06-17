//output the matrix
template<typename R>
string matrix<R>::output(){
	string res = "rank=";
	res += std::to_string(rank) + "\n";
	for(unsigned i=0; i<rank; ++i)
		res += moduleOper->output(find(i)) + "\n";
	return res;
}

//output with specified format
template<typename R>
string matrix<R>::output(string name0, string name1){
	string res = "";
	for(unsigned i=0; i<rank; ++i)
		res += "[" + name0 + "-" + std::to_string(i) + "]\t->\t" + moduleOper->output(find(i),name1) + "\n";
	return res;
}

//save a matrix. 4 bytes of rank, followed by the rows
template<typename R>
void matrix<R>::save(std::iostream &writer){
	writer.write((char*)&rank, 4);
	for(unsigned i=0; i<rank; ++i)
		moduleOper->save(find(i),writer);
}

//load a matrix with known rank
template<typename R>
void matrix<R>::load(std::iostream &reader, int rk) {
	clear();
	set_rank(rk);
	for(unsigned i=0; i<rank; ++i){
		auto newrow = moduleOper->load(reader);
		insert(i,newrow);
	}
}

//load a matrix
template<typename R>
void matrix<R>::load(std::iostream &reader) {
	reader.read((char*)&rank, 4);
	load(reader, rank);
}

//load and modify a matrix
template<typename R>
void matrix<R>::load_modify(std::iostream &reader, std::function<vectors<matrix_index,R>(const vectors<matrix_index,R>&)> rule){
	int32_t rk;
	reader.read((char*)&rk, 4);
	clear();
	set_rank(rk);
	for(unsigned i=0; i<rank; ++i){
		auto newrow = moduleOper->load(reader);
		insert(i,rule(newrow));
	}
}
