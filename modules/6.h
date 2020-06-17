//reindex a vector
template<typename index, typename R>
template<typename index2>
vectors<index2,R> ModuleOp<index,R>::re_index(std::function<index2(index)> rule, vectors<index,R> const &x){
	vectors<index2,R> res;
	res.dataArray.resize(x.size());
	
	for(unsigned i=0; i<x.size(); ++i){
		res.dataArray[i].ind = rule(x.dataArray[i].ind);
		res.dataArray[i].coeficient = x.dataArray[i].coeficient;
	}
	res.sort();
	return res;
}

//desctructive reindex
template<typename index, typename R>
template<typename index2>
vectors<index2,R> ModuleOp<index,R>::re_index(std::function<index2(index)> rule, vectors<index,R> &&x){
	vectors<index2,R> res;
	res.dataArray.resize(x.size());
	for(unsigned i=0;i<x.size();++i){
		res.dataArray[i].ind = rule(x.dataArray[i].ind);
		res.dataArray[i].coeficient = std::move(x.dataArray[i].coeficient);
	}
	res.sort();
	return res;
}

//reindex a vector, filtering out those unwanted terms, the rule produces the invalid index to indicate an unwanted term
template<typename index, typename R> 
template<typename index2>
vectors<index2,R> ModuleOp<index,R>::filtered_reindex(std::function<index2(index)> rule, vectors<index,R> const &x, index2 invalid){
	vectors<index2,R> res;
	res.dataArray.reserve(x.size());
	
	for(unsigned i=0; i<x.size(); ++i){
		auto k = rule(x.dataArray[i].ind);
		if(k==invalid) continue;
		res.dataArray.push_back({ k, x.dataArray[i].coeficient});
	}
	res.sort();
	return res;
}

//filter a vector
template<typename index, typename R>
vectors<index,R> ModuleOp<index,R>::filter(std::function<bool(index)> rule, vectors<index,R> const &x){
	vectors<index,R> res;
	res.dataArray.reserve(x.size());
	
	for(unsigned i=0; i<x.size(); ++i)
		if(rule(x.dataArray[i].ind))
			res.dataArray.push_back(x.dataArray[i]);
	return res;
}

//filter when the rule involves the coeficient of a term
template<typename index, typename R>
vectors<index,R> ModuleOp<index,R>::filter(std::function<bool(index,const R&)> rule, vectors<index,R> const &x){
	vectors<index,R> res;
	res.dataArray.reserve(x.size());
	
	for(unsigned i=0;i<x.size();++i)
		if(rule(x.dataArray[i].ind,x.dataArray[i].coeficient))
			res.dataArray.push_back(x.dataArray[i]);
	return res;
}

//apply an operation termwise
template<typename index, typename R>
template<typename index2>
vectors<index2,R> ModuleOp<index,R>::termwise_operation(std::function<std::pair<index2,R>(index,const R&)> rule, vectors<index,R> const &x){
	vectors<index2,R> res;
	res.dataArray.resize(x.size());
	for(unsigned i=0; i<x.size(); ++i){
		auto newT = rule(x.dataArray[i].ind, x.dataArray[i].coeficient);
		res.dataArray[i].ind = newT.first;
		res.dataArray[i].coeficient = std::move(newT.second);
	}
	res.sort();
	return res;
}

//apply an operation on each coeficient
template<typename index, typename R>
template<typename R2>
vectors<index,R2> ModuleOp<index,R>::termwise_operation(std::function<R2(const R&)> rule, vectors<index,R> const &x){
	vectors<index,R2> res;
	res.dataArray.resize(x.size());
	for(unsigned i=0;i<x.size();++i){
		auto newT = rule(x.dataArray[i].coeficient);
		res.dataArray[i].ind = x.dataArray[i].ind;
		res.dataArray[i].coeficient = std::move(newT);
	}
	return res;
}
