//search the coeficient with index n in a range
template<typename index, typename R>
R ModuleOp<index,R>::component(index n, vectors<index,R> const &x, unsigned startPos, unsigned endPos){
	//not found
	if(endPos<=startPos)
		return ringOper->zero();
	if(endPos-startPos==1){
		//found
		if(x.dataArray[startPos].ind==n)
			return x.dataArray[startPos].coeficient;
		//not found
		else return ringOper->zero();
	}
	
	//recursion
	unsigned midPos = (startPos+endPos)/2;
	if(x.dataArray[midPos].ind==n)
		//found
		return x.dataArray[midPos].coeficient;
	//n is in the left half
	else if(x.dataArray[midPos].ind>n)
		return component(n,x,startPos,midPos);
	//n in the right half
	else return component(n,x,midPos+1,endPos);
}

//search for term with index n
template<typename index, typename R>
R ModuleOp<index,R>::component(index n, vectors<index,R> const &x){ 
	return component(n,x,0,x.dataArray.size()); }
	
//the term next to the term with index n, in a range
template<typename index, typename R>
unsigned ModuleOp<index,R>::next_non_trvial_entry(index n, vectors<index,R> const &x, unsigned startPos, unsigned endPos){
	if(endPos<=startPos) return endPos;
	if(endPos-startPos==1){
		if(x.dataArray[startPos].ind>=n) 	
			return startPos;
		else return endPos;
	}
	unsigned midPos = (startPos+endPos)/2;
	if(x.dataArray[midPos].ind==n) 
		return midPos;
	else if(x.dataArray[midPos].ind<n) 
		return next_non_trvial_entry(n,x,midPos+1,endPos);
	else 
		return next_non_trvial_entry(n,x,startPos,midPos);
}
	
//the index next to the term with index n
template<typename index, typename R>
index ModuleOp<index,R>::next_non_trvial_entry(index n, vectors<index,R> const &x){
	unsigned p = next_non_trvial_entry(n,x,0,x.dataArray.size());
	if(p >= x.dataArray.size())
		return last_entry(x)+1;
	return x.dataArray[p].ind;
}
	
//find the first invertible term
template<typename index, typename R>
unsigned ModuleOp<index,R>::first_invertible(vectors<index,R> const &x){
	unsigned i;
	for(i=0; i<x.size(); ++i)
		if(ringOper->invertible(x.dataArray[i].coeficient))
			break;
	return i;
}
	
//find the first invertible index
#define INVALIDLY_LARGE ((1<<30) - 1)
template<typename index, typename R>
index ModuleOp<index,R>::first_invertible_index(vectors<index,R> const &x){
	unsigned k = first_invertible(x);
	if(k>=x.size()) 
		return INVALIDLY_LARGE;
	return x.dataArray[k].ind;
}
	
//the index of the last entry
template<typename index, typename R>
index ModuleOp<index,R>::last_entry(vectors<index,R> const &x){
	if(x.dataArray.size()==0) return -100;
	else return x.dataArray.rbegin()->ind;
}
