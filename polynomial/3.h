//operations with paralell algorithms
template<typename exponent_type, typename base_ring>
class PolyOp_Para : public PolyOp<exponent_type,base_ring>{
public:
	//multiply with monomials
	virtual poly<exponent_type,base_ring> mon_multiply(const poly<exponent_type,base_ring> &x, exponent_type e, const base_ring &r);

	//constructor
	PolyOp_Para(RingOp<base_ring> *ringop);
	
	//multiplication
	poly<exponent_type,base_ring> multiply(const poly<exponent_type,base_ring> &x, const poly<exponent_type,base_ring> &y);
};

template<typename base_ring>
using PolynomialOp_Para = PolyOp_Para<exponent,base_ring>;
