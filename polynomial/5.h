//constructor
template<typename exponent_type, typename base_ring>
PolyOp_Para<exponent_type,base_ring>::PolyOp_Para(RingOp<base_ring> *ringop) : ModuleOp<exponent_type,base_ring>(ringop), PolyOp<exponent_type,base_ring>(ringop){}

//multiply with a monomial
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp_Para<exponent_type,base_ring>::mon_multiply(const poly<exponent_type,base_ring> &x, exponent_type e, const base_ring &r) {
	poly<exponent_type,base_ring> result;
	result.dataArray.resize(x.size());
	
#pragma omp parallel for schedule(dynamic)
	for(unsigned i=0; i<x.size(); ++i){
		result.dataArray[i].ind = x.dataArray[i].ind + e;
		result.dataArray[i].coeficient = this->ringOper->multiply(x.dataArray[i].coeficient, r);
	}
	return result;
}

//multiply two polys
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp_Para<exponent_type,base_ring>::multiply(const poly<exponent_type,base_ring> &x, const poly<exponent_type,base_ring> &y){
	if(this->isZero(x)) return this->zero();
	if(this->isZero(y)) return this->zero();
	
	//make y smaller
	const poly<exponent_type,base_ring> *x1,*x2;
	if(x.size()<y.size()){
		x1=&y;
		x2=&x;
	}
	else {
		x1=&x;
		x2=&y;
	}
	
	//take the sum in a paralell way
	std::function<poly<exponent_type,base_ring>(int)> sumds = [x1,x2,this](int i){ 
		return mon_multiply(*x1,x2->dataArray[i].ind,x2->dataArray[i].coeficient); };
	return this->sum(sumds, 0, x2->size());
}
