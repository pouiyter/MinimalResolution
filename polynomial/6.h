using std::to_string;
//the connstructor
template<typename exponent_type, typename base_ring>
PolyOp<exponent_type,base_ring>::PolyOp(RingOp<base_ring> *ringop) : ModuleOp<exponent_type,base_ring>(ringop){
	 this->ringOper = ringop;
	 this->output_term = [ringop](exponent_type x, base_ring y){ 
		 return "(" + ringop->output(y) + ")" + "x" + "^" + to_string(x); };
}

//construct constant poly
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp<exponent_type,base_ring>::constant(const base_ring &r){ 
	return this->singleton(0,r); }

//check if a poly is constant
template<typename exponent_type, typename base_ring>
bool PolyOp<exponent_type,base_ring>::isConst(const poly<exponent_type,base_ring> &x){ 
	return x.size()==1 && x.dataArray[0].ind==0; }
	
//the unit map from integers
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp<exponent_type,base_ring>::unit(int i){ 
	return constant(this->ringOper->unit(i)); }
	
//check if an element is invertible, in this case means a poly of degree 0 with invertible value
template<typename exponent_type, typename base_ring>
bool PolyOp<exponent_type,base_ring>::invertible(const poly<exponent_type,base_ring> &x){ 
	return isConst(x) && this->ringOper->invertible(x.dataArray[0].coeficient); }
	
//compute the inverse of an invertible element
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp<exponent_type,base_ring>::inverse(poly<exponent_type,base_ring> const &x){ 
	return constant(this->ringOper->inverse(x.dataArray[0].coeficient)); }
	
//construction a monomial
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp<exponent_type,base_ring>::monomial(exponent_type e, const base_ring &r){ 
	return this->singleton(e,r); }
	
//construct a monomial with unit coeficient
template<typename exponent_type, typename base_ring>
poly<exponent_type,base_ring> PolyOp<exponent_type,base_ring>::monomial(exponent_type e){ 
	return this->singleton(e,this->ringOper->unit(1)); }
