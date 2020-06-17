//multiplication with a monomial
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp<exponent_type,base_ring>::mon_multiply(const poly<exponent_type,base_ring> &x, exponent_type e, const base_ring &r){
	poly<exponent_type,base_ring> result;
	result.dataArray.resize(x.size());
	for(unsigned i=0;i<x.size();++i){
		result.dataArray[i].ind = x.dataArray[i].ind + e;
		result.dataArray[i].coeficient = this->ringOper->multiply(x.dataArray[i].coeficient, r);
	}
	return result;
}

//multiplication algorithm
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp<exponent_type,base_ring>::mult(const poly<exponent_type,base_ring> &x, const poly<exponent_type,base_ring> &y, int startPos, int endPos){
	if(endPos <= startPos) return this->zero();
	//if there is only one term, we take the scalor multiplication
	if(endPos-startPos==1) return mon_multiply(x,y.dataArray[startPos].ind,y.dataArray[startPos].coeficient);
	int midPos = (startPos + endPos) / 2;
	//compute the multiplication of two parts
	poly<exponent_type,base_ring> s1,s2;
	s1 = mult(x,y,startPos,midPos);
	s2 = mult(x,y,midPos,endPos);
	//add the two parts
	return this->add(std::move(s1), std::move(s2));
}

//multiplication
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp<exponent_type,base_ring>::multiply(const poly<exponent_type,base_ring> &x, const poly<exponent_type,base_ring> &y){ 
	//choose y to have smaller size
	if(x.size()>=y.size())
		return mult(x,y,0,y.size()); 
	else return mult(y,x,0,x.size());
}
