using std::to_string;

//output the vector
template<typename index, typename R>
string ModuleOp<index,R>::output(vectors<index,R> x){
	auto it1 = x.dataArray.begin();
	string outs = "o";
	while(it1!=x.dataArray.end()) {
		outs += "+" + output_term(it1->ind,it1->coeficient);
		it1++;
	}
	return outs;
}

//output vectors with the format a[n-i], where n is a specified name
template<typename index, typename R>
string ModuleOp<index,R>::output(vectors<index,R> x, string name){
	auto it1 = x.dataArray.begin();
	string outs = "o";
	while(it1!=x.dataArray.end()) {
		outs += "+" + ringOper->output(it1->coeficient) + "[" + name + "-" + to_string(it1->ind) + "]";
		it1++;	
	}
	return outs;
}

//save the vector, with format: 4 bytes for the number of terms, then the array of terms, each term the index followed by coeficient
template<typename index, typename R>
void ModuleOp<index,R>::save(vectors<index,R> const &x, std::iostream& writer){
	int32_t length = x.size();
	writer.write((char*)&length, 4);
	
	for(int i=0; i<length; ++i){
		writer.write((char*)&x.dataArray[i].ind, sizeof(index));
		ringOper->save(x.dataArray[i].coeficient,writer);
	}
}

//load the vector
template<typename index, typename R>
vectors<index,R> ModuleOp<index,R>::load(std::iostream& reader){
	int32_t length;
	reader.read((char*)&length,4);
	
	vectors<index,R> result;
	for(int i=0; i<length; i++) {
		typename vectors<index,R>::term newT;
		reader.read((char*)&newT.ind, sizeof(index));
		newT.coeficient = ringOper->load(reader);
		result.push(std::move(newT));
	}
	return result;
}
