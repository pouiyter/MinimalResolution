//compute the sum of a array of elements
template<typename return_type>
return_type sum(std::function<return_type(int)> summands, int start, int end, std::function<return_type(return_type&&, return_type&&)> add){
	if(end-start==1) return summands(start);
	return_type s1,s2;
	int mid = (start+end)/2;
	#pragma omp parallel sections
	{
		#pragma omp section
		{ s1 = sum(summands, start,mid, add); }
		#pragma omp section
		{ s2 = sum(summands, mid,end, add); }
	}
	return add(std::move(s1), std::move(s2));
}

//the sum of an array of elements in an abelian group
template<typename A>
A sum(std::function<A(int)> summands, int start, int end, AbGroupOp<A> *abOper){
	std::function<A(A&&, A&&)> adding = [abOper] (A &&x, A &&y){ 
		return abOper->add(std::move(x),std::move(y)); };
	if(end<=start) return abOper->zero();
	return sum(summands, start, end, adding);
}

//the sum of an array of vectors
template<typename index, typename R>
vectors<index,R>  ModuleOp<index,R>::sum(std::function<vectors<index,R>(int)> summands, int start, int end){
	return ::sum(summands, start, end, this);
}
