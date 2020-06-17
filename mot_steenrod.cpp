#include "mot_steenrod.h"

//if x is zero
inline bool tauOper::isZero(tauPoly const &x) { return x==internal_zero; }

//add two elements
inline tauPoly tauOper::add(tauPoly const &x, tauPoly const &y) { 
	if(isZero(x)) return y;
	if(isZero(y)) return x;
	//if x and y are both nontrivial, they are addible only when they have the same degree, and then x=y in F2[tau]
	return internal_zero;
}

//addition
inline tauPoly tauOper::add(tauPoly &&x, tauPoly &&y){ 
	return add(x,y); }
	
//zero element
inline tauPoly tauOper::zero(){ return internal_zero; }

//minus operation, which is trivial modolu two
inline tauPoly tauOper::minus(tauPoly const &x){ return x; }

//multiplication
inline tauPoly tauOper::multiply(tauPoly const &x, tauPoly const &y){
	if(isZero(x) || isZero(y)) return internal_zero;
	return x+y;
}

//the unit map, 0 means the internal degree 0, which is the unit
inline tauPoly tauOper::unit(int n){
	if(n%2==0) return internal_zero;
	else return 0;
}

//only degree 0 elements are invertible
inline bool tauOper::invertible(tauPoly const &x){ return x==0; }

//the inverse of invertibles
inline tauPoly tauOper::inverse(tauPoly const &x){ return -x; }
    
//output in the form t^a
inline string tauOper::output(tauPoly x){
	if(isZero(x)) return "0";
	else if(x==0) return "1";
	else return "t^" + std::to_string(x);
}
    
//save, with 2 bytes
inline void tauOper::save(const tauPoly &x, std::iostream &writer){
    writer.write((char*)&x, sizeof(tauPoly));
}    

//load the 2 bytes
inline tauPoly tauOper::load(std::iostream &reader){
	tauPoly reault;  reader.read((char*)&reault, sizeof(tauPoly));
	return reault;
}

//the tau-valuation
int tauOper::tau_valuation(tauPoly x){ return x;}

//shift the exponent of tau
void tauOper::shift(tauPoly& x, int n){ x+=n;}

//the power of tau
tauPoly tauOper::power_tau(int i){ return i;}
    
//operations of F2[tau]
tauOper tau_oper;
//operations on Amot
MotSteenrodRingOp motSteenrod_oper(&tau_oper);
//operations on A\otimes A
PolynomialOp<motSteenrod> mStmSt_oper(&motSteenrod_oper);

//constructor
MotSteenrodRingOp::MotSteenrodRingOp(tauOper *tauop) : ModuleOp<exponent,tauPoly>(tauop), PolynomialOp<tauPoly>(tauop){
	this->output_term = [tauop](exponent e, tauPoly r){ 
		  return tauop->output(r) + ::output(e); };
}
    
//add the degrees
MotDegree MotDegree::add(MotDegree x, MotDegree y){
	return {x.deg+y.deg, x.weight+y.weight};
}

//add the degrees
template<>
MotDegree cofree_comodule<motSteenrod,MotDegree>::add_degree(MotDegree x,MotDegree y){
	return MotDegree::add(x,y);
}

//the underlying degree
template<>
int cofree_comodule<motSteenrod,MotDegree>::underlyingDeg(MotDegree x){
	return x.deg;
}

//compare
bool MotDegree::operator<(const MotDegree y){
	if(deg<y.deg)
		return true;
	if(deg>y.deg)
		return false;
	return weight < y.weight;
}

//the degree of the n-th variable
int MotDegree::varDeg(int n){ 
      return xnDegs[n]; }

//the weight of n-th variable
int MotDegree::varWeight(int n){ return varDeg(n)/2; }

//the degree of e
int MotDegree::expoDeg(exponent e){ 
	return total_deg(e,varDeg); }

//weight of e
int MotDegree::expoWeight(exponent e){ 
	return total_deg(e,varWeight); }

//output the bidegree
string MotDegree::output(){ 
	return "(" + std::to_string(deg) + "," + std::to_string(weight) + ")"; 
}

//modules over A
ModuleOp<matrix_index,motSteenrod> motSteenrod_module_oper(&motSteenrod_oper);  

//modules over Ctau
ModuleOp<matrix_index,tauPoly> tau_module_oper(&tau_oper);

//constructor
MotSteenrodOp::MotSteenrodOp(matrix<motSteenrod> *cofa, int max_degree){
	//set the matrix for the comultiplication
	cofree_coaction = cofa;

	//initialize operations
	this->ringOper = &tau_oper;
	this->algebroidRingOper = &motSteenrod_oper;
	this->moduleOper = &tau_module_oper;
	this->algebroidModuleOper = &motSteenrod_module_oper;
    
	//set the maximal degree
	this->maxDeg = max_degree;
}
       
// The tau valuation of the monomial e. x_i^2 can be divided by tau, so has valuation 1
int MotSteenrodOp::tauVal(exponent e) {  
	int res=0;
	for(int i=1;i<=maxVar;i++)
		res += xnVal(e,i)/2;
	return res;
}
    
//the degree of the generators in a cofree comodule
MotDegree MotSteenrodOp::cofree_degs(matrix_index n) {
	exponent e = mon_array[n];
	auto degg=MotDegree::expoDeg(e);
	auto wett=MotDegree::expoWeight(e)+tauVal(e);
	return {degg,wett};
}

//transform a element in A into a vector
vectors<matrix_index,tauPoly> MotSteenrodOp::algebroid2vector(motSteenrod const &x,int shift){
	std::function<std::pair<matrix_index,tauPoly>(exponent, const tauPoly&)> rl = [this,shift](exponent e, const tauPoly &r){ 
	//transform x^e to val(e)[index(e)]
		return std::pair<matrix_index,tauPoly>(mon_index[e]+shift, this->ringOper->multiply(tauPoly(tauVal(e)),r)); };
	vectors<matrix_index,tauPoly> result = motSteenrod_oper.termwise_operation(rl,x);
	return result;
}
    
//transform a polynomial in A into a vector in A
vectors<matrix_index,motSteenrod> MotSteenrodOp::algebroid2vectorSt(mStmSt const &x){
	std::function<std::pair<matrix_index,motSteenrod>(exponent, const motSteenrod&)> rl = [this](exponent e, const motSteenrod &r){ 
		return std::pair<matrix_index,motSteenrod>(mon_index[e], this->algebroidRingOper->multiply(this->etaL(tauPoly(tauVal(e))),r)); };
	return mStmSt_oper.termwise_operation(rl,x);
}
    
//transform a vactor back
motSteenrod MotSteenrodOp::vector2algebroid(vectors<matrix_index,tauPoly> const &x){
	std::function<std::pair<exponent,tauPoly>(matrix_index,const tauPoly&)> rl = [this](matrix_index n, const tauPoly &r){ 
		exponent newE = mon_array[n];
		tauPoly val = tauPoly(-tauVal(newE));
		return std::pair<exponent,tauPoly>(newE, this->ringOper->multiply(val,r)); };
	return tau_module_oper.termwise_operation(rl,x);
}
    
//generate the cofree coaction
void MotSteenrodOp::generate_cofree_coaction(string coaction_filename, string expo_filename){
	//transform the data of coactions
	std::fstream coaction_file(coaction_filename, std::ios::in | std::ios::binary);
	std::fstream expo_file(expo_filename, std::ios::in | std::ios::binary);
	int32_t nums;
	expo_file.read((char*)&nums, 4);
	if((unsigned)nums != mon_array.size())
		std::cerr << "datas not compatible!" << std::flush;
	for(int j=0; j<nums; j++){
		mStmSt coactor = mStmSt_oper.load(coaction_file);
		exponentArry eps = load_expArry(expo_file);
		exponent e = pack(eps);
		auto n = mon_index[e];
		//compute the coaction on tau^{-val(e)}x^e
		auto so = motSteenrod_module_oper.scalor_mult(etaL(tauPoly(-tauVal(mon_array[n]))), algebroid2vectorSt(coactor));
		this->cofree_coaction->insert(n, so);
	}
      
	//initialize the degrees
	static std::function<MotDegree(matrix_index)> degs =  [this](matrix_index n){ 
		return this->cofree_degs(n); };
	this->init_cofree_data<MotDegree>(degs);
}

//the left and right units
motSteenrod MotSteenrodOp::etaR(const tauPoly &x){ 
	return motSteenrod_oper.monomial(0,x);}
motSteenrod MotSteenrodOp::etaL(const tauPoly &x){ 
	return etaR(x);}

//number of monomials
unsigned MotSteenrodOp::ranksBelowDeg(unsigned n){ 
	return ranksbelow[n]; }
	
//the comultiplication
vectors<matrix_index, motSteenrod> MotSteenrodOp::delta(matrix_index n){ 
	return cofree_coaction->find(n); }
   
//the Hopf classes
motSteenrod MotSteenrodOp::hi(int i){ 
	return motSteenrod_oper.monomial(1 << (i), - (1 << (i-1))); }

    
//initialize the list of monoials
void MotSteenrodOp::init_mon_array(string filename){
	std::fstream mons_file(filename, std::ios::in | std::ios::binary);
	//get the number of monoials
	int32_t nums;
	mons_file.read((char*)&nums, 4);
	ranksbelow.clear();
	ranksbelow.resize(maxDeg+1);
	mon_array.resize(nums);
	for(int i=0; i<nums; ++i){
		//load the i-th exponent
		exponentArry eps = load_expArry(mons_file);
		exponent e = pack(eps);
		//update the list of monomials
		mon_array[i] = e;
		//update the index
		mon_index.emplace(e, i);
		//update the number of monpmials
		for(int d=0; d<=(int)maxDeg; ++d)
			if(MotDegree::expoDeg(e) <= d)
				ranksbelow[d]++;
	}
}

//output the array of monomials
string MotSteenrodOp::output_monomials(){
	string res = "monomials:\n";
	for(exponent e: mon_array)
		res += ::output(e) + "\n";
	res += "ranks:\n";
	for(auto n:ranksbelow)
		res += std::to_string(n) + "\n";
	return res;
}

//lift elements in F2 to  F2[tau]
tauPoly MotSteenrodOp::lift(F2 x){
	if(x%2 != 0) return tau_oper.unit(1);
	return tau_oper.zero();
}

//lift vectors over F2 to vectors over F2[tau]
vectors<matrix_index,tauPoly> MotSteenrodOp::lift(const vectors<matrix_index,F2>&v){
	vectors<matrix_index, tauPoly> result;
	for(auto tm : v.dataArray)
		result.push({tm.ind, lift(tm.coeficient)});
	return result;
}

