//swap the inside and out side index
template<typename index1, typename index2, typename R> 
vectors<index1,vectors<index2,R>> swapping(vectors<index2,vectors<index1,R>> const &x, std::function<vectors<index1,vectors<index2,R>>(vectors<index1,vectors<index2,R>> &&, vectors<index1,vectors<index2,R>> &&)> merging){
	//swap each entry of the outside vector
	std::function<vectors<index1,vectors<index2,R>>(int)> entries = [merging, &x] (int i) {
		vectors<index1,vectors<index2,R>> result;
		result.dataArray.resize(x.dataArray[i].coeficient.dataArray.size());
		for(unsigned k=0; k<x.dataArray[i].coeficient.dataArray.size(); ++k){
			//construct the inside term of the swapped vector
			typename vectors<index2,R>::term newT(x.dataArray[i].ind, x.dataArray[i].coeficient.dataArray[k].coeficient);
			//set the new outside index
			result.dataArray[k].ind = x.dataArray[i].coeficient.dataArray[k].ind;
			//set the new inside vector
			result.dataArray[k].coeficient.push(std::move(newT));
		}
		
		return result;
	};
	//sum up the swaped terms
	return sum(entries,0,x.dataArray.size(),merging);
}

template<typename index1, typename index2, typename R> 
vectors<index1,vectors<index2,R>> swapping(vectors<index2,vectors<index1,R>> const &x, AbGroupOp<vectors<index1,vectors<index2,R>>> *Abop){
	std::function<vectors<index1,vectors<index2,R>>(vectors<index1,vectors<index2,R>> &&, vectors<index1,vectors<index2,R>> &&)> merging = [Abop](vectors<index1,vectors<index2,R>> &&x, vectors<index1,vectors<index2,R>> &&y){
		return Abop->add(std::move(x), std::move(y)); };
	return swapping(x, merging);
}
