//a poly is stored as a vector
template<typename exponent_type, typename base_ring>
using poly = vectors<exponent_type,base_ring>;

//operations on polys
template<typename exponent_type, typename base_ring>
class PolyOp :
virtual public ModuleOp<exponent_type,base_ring>, 
virtual public RingOp<poly<exponent_type,base_ring>>{
public:
	//multiply with a monomial
	virtual poly<exponent_type,base_ring> mon_multiply(const poly<exponent_type,base_ring> &x, exponent_type e, const base_ring &r);
	
	//construction of constant poly
	poly<exponent_type,base_ring> constant(const base_ring &r);
	
	//check if a poly is constant
	bool isConst(const poly<exponent_type,base_ring> &x);
	
private:
	//multiplication
	poly<exponent_type,base_ring> mult(const poly<exponent_type,base_ring> &x, const poly<exponent_type,base_ring> &y, int startPos, int endPos);
	
public:
	//the constructor
	PolyOp(RingOp<base_ring> *ringop);
	
	//multiplication
	poly<exponent_type,base_ring> multiply(const poly<exponent_type,base_ring> &x, const poly<exponent_type,base_ring> &y);
	
	//unit map from integers
	poly<exponent_type,base_ring> unit(int i);
	
	//check if an element is invertible
	bool invertible(const poly<exponent_type,base_ring> &x);
	//the inverse
	poly<exponent_type,base_ring> inverse(poly<exponent_type,base_ring> const&);
	
	//construct a monomial
	poly<exponent_type,base_ring> monomial(exponent_type e, const base_ring &r);
	
	//monomial with unit coeficient
	poly<exponent_type,base_ring> monomial(exponent_type e);
};

//type of the exponent
typedef uint32_t exponent;

//type of ordinary poly
template<typename base_ring>
using PolynomialOp = PolyOp<exponent, base_ring>;

template<typename base_ring>
using polynomial = poly<exponent,base_ring>;
