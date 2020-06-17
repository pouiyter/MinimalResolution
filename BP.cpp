//BP.cpp
#include"BP.h"
#include"mon_index.h"
#include<ios>

//constructor
BP_Op::BP_Op(int maxdeg, Z2_Op *Z2_op, matrix<BP> *etaL_mat, matrix<BPBP> *delta_mat, matrix<BP> *R2L_mat) :
ModuleOp<exponent,Z2>(Z2_op), PolynomialOp_Para<Z2>(Z2_op), BPBP_opers(this), BPBPBP_opers(&BPBP_opers), BPMod_opers(this), BPBPMod_opers(&BPBP_opers), Z2Mod_opers(Z2_op), mon_index(maxdeg){
	mon_index.init_mon_array();
	maxDeg = maxdeg;
	Hopf_Algebroid<BP,BPBP>::ringOper = this;
	algebroidRingOper = &BPBP_opers;
    
	moduleOper = &BPMod_opers;
	algebroidModuleOper = &BPBPMod_opers;
    
	Z2_oper = Z2_op;
	etaL_table = etaL_mat;
	delta_table = delta_mat;
	R2L_table = R2L_mat;
}

//constructor
BPBP_Op::BPBP_Op(BP_Op *BP_op) : ModuleOp<exponent,BP>(BP_op), PolynomialOp_Para<BP>(BP_op){}

//constrctor
BPBPBP_Op::BPBPBP_Op(BPBP_Op *BPBP_op) : ModuleOp<exponent,BPBP>(BPBP_op), PolynomialOp_Para<BPBP>(BPBP_op){}

//load etaL table
void BP_Op::load_etaL(string filename){  
	std::cout << "loading etaL table from " << filename << "...\n" << std::flush;
	std::fstream reader(filename, std::ios::in | std::ios::binary);
	etaL_table->load(reader, mon_index.number_of_all_mons()); 
	std::cout << "etaL data loaded\n" << std::flush;
}

//load R2L talbe
void BP_Op::load_R2L(string filename){   
	std::cout << "loading R2L table...\n" << std::flush;
	std::fstream reader(filename, std::ios::in | std::ios::binary);
	R2L_table->load(reader, mon_index.number_of_all_mons());
	std::cout << "R2L data loaded\n" << std::flush;
}

//load delta table
void BP_Op::load_delta(string filename){    
	std::cout << "loading delta table...\n" << std::flush;
	std::fstream reader(filename, std::ios::in | std::ios::binary);
	delta_table->load(reader, mon_index.number_of_all_mons());
	std::cout << "delta data loaded\n" << std::flush;
}

//the left unit
BPBP BP_Op::etaL(const BP &x){
	std::function<BPBP(const Z2&, const BPBP&)> scalor_mult = [this] (const Z2 &r, const BPBP &fm){
		return BPBP_opers.scalor_mult(this->monomial(0,r),fm); };
	auto xvec = mon_index.poly2vec(x, this);
	return etaL_table->maps_to(xvec, scalor_mult, &BPBP_opers);
}

//change the left notation to the right notation
BPBP BP_Op::etaL(const BPBP &x){
	std::function<BPBP(const BP&, const BPBP&)> scalor_mult = [this](const BP &r, const BPBP &fm){
		return BPBP_opers.scalor_mult(r,fm); };
	auto xvec = mon_index.poly2vec(x,&BPBP_opers);
	return etaL_table->maps_to(xvec,scalor_mult, &BPBP_opers);
}

//the right unit, vn is in the outer
BPBP BP_Op::etaR(const BP &x){
	std::function<BP(const Z2&)> tfm = [this](const Z2& r){
		return this->monomial(0,r); };
	return this->termwise_operation(tfm,x);
}

//change algebroid to a vector
vectors<matrix_index, BP> BP_Op::algebroid2vector(const BPBP& x, int shift){
	std::function<matrix_index(exponent)> rd = [this,shift](exponent e){
		return mon_index.mon_index[e]+shift; };
	return BPBP_opers.re_index(rd,std::move(R2L(x)));
}

//change right notation to the left notation and switch ti to the outer
BPBP BP_Op::R2L(const BPBP &x){
	std::function<BPBP(const BP&, const BPBP&)> scalor_mult = [this](const BP &r, const BPBP &fm){
		std::function<BP(const Z2&)> tfm = [this](const Z2& r0){
			return this->monomial(0,r0); };
		auto r1 = termwise_operation(tfm,r);
		return BPBP_opers.multiply(r1,fm); };
	auto xvec = mon_index.poly2vec(x,&BPBP_opers);
	return R2L_table->maps_to(xvec,scalor_mult, &BPBP_opers);
}

//change back
BPBP BP_Op::vector2algebroid(const vectors<matrix_index, BP>&x){
    std::cerr << "not implemented!" ;
    return BPBP_opers.zero();
}
    
//the coaction
vectors<matrix_index, BPBP> BP_Op::delta(matrix_index n){
    return delta_table->find(n);
}

//number of generators
unsigned BP_Op::ranksBelowDeg(unsigned n){ 
	return mon_index.ranksBelow[n]; }

//initialize
void BP_Op::initialize(string etaL_filename, string R2L_filename, string delta_filename){
	load_etaL(etaL_filename);
	load_R2L(R2L_filename);
	load_delta(delta_filename);
    
	//initialize the degree of generators
	std::function<int(matrix_index)>  cofree_degs = [this](matrix_index n){
		exponent e = mon_index.mon_array[n];
		return total_deg(e,xnDeg);
	};
        
	//initialize the data for a cofree coalgebra
	this->init_cofree_data(cofree_degs);
}

//make structure tables
void BP_Op::make_tables(int maxVar, string R2Lfilename, string deltafilename, string etaL_filename, string R2L_filename, string delta_filename, std::ostream &outputfile){
	//the values on the generators
	std::vector<BPBP> R2L_gen(maxVar+1);
	std::vector<BPBPBP> delta_gen(maxVar+1);
	std::vector<BPBP> etaL_gen(maxVar+1);
	
	//load the values on the generators. 
	//for etaR, At this point we use the left unit expression, and the vi are on the outer, ti on the innner
	std::fstream R2L_file(R2Lfilename, std::ios::in | std::ios::binary);
	R2L_gen[0] = BPBP_opers.unit(1);
	for(int i=1; i<=maxVar; ++i)
		R2L_gen[i] = BPBP_opers.load(R2L_file);
	R2L_file.close();
	std::cout << "etaR data loaded\n";
	for(int i=1; i<=maxVar; ++i)
		outputfile << "etaR(v" << i << ") = " << BPBP_opers.output(R2L_gen[i]) << "\n";
	
	//for delta, at this point, we use left unit expression, the vi in the middle, right ti in the outer, the left ti in the inner
	std::fstream delta_file(deltafilename, std::ios::in | std::ios::binary);
	delta_gen[0] = BPBPBP_opers.unit(1);
	for(int i=1; i<=maxVar; ++i)
		delta_gen[i] = BPBPBP_opers.load(delta_file);
	delta_file.close();
	std::cout << "delta data loaded\n";
	for(int i=1; i<=maxVar; ++i)
		outputfile << "delta(t" << i << ") = " << BPBPBP_opers.output(delta_gen[i]) << "\n";
	
	//compute the generators for etaL, in right unit expression,
	etaL_gen[0] = BPBP_opers.unit(1);
	for(int i=1; i<=maxVar; ++i){
		PolynomialOp_Para<BPBP> PBPBPoper(&BPBP_opers);
		polynomial<BPBP> etaL_rule = PBPBPoper.constant(BPBP_opers.monomial(vars(i)));
		//etR = -etaR(vn) + etaL(vn) in left expression , with vi in the outer
		auto etR = R2L_gen[i];
		etR = BPBP_opers.minus(etR);
		etR = BPBP_opers.add(BPBP_opers.monomial(vars(i)), etR);
		
		//lift etR to a polynomial
		std::function<BPBP(const BP&)> ct = [this](const BP& x){ 
			return BPBP_opers.constant(x); };
		auto tR = BPBP_opers.termwise_operation(ct,etR);
		
		//etaL(v_n) = -etaR(vn) + etaL(vn) + etaR(vn)
		etaL_rule = PBPBPoper.add(etaL_rule, tR);
		
		etaL_gen[i] = substitute(etaL_gen, etaL_rule, &BPBP_opers);
	}
	std::cout << "etaL data computed\n";
	for(int i=1; i<=maxVar; ++i)
		outputfile << "etaL(v" << i << ") = " << BPBP_opers.output(etaL_gen[i]) << "\n";
	
	//compute the etaL table
	std::function<BPBP(int)> etaL_gens = [&etaL_gen](int i){
		return etaL_gen[i]; };
	std::fstream etaL_tablefile(etaL_filename, std::ios::out | std::ios::binary);
	mon_index.substitution_table(etaL_gens, etaL_tablefile, &BPBP_opers);
	etaL_tablefile.close();
	
	//compute the R2L table
	std::function<BPBP(int)> R2L_gens = [&R2L_gen, this](int i){
		//switch vi to the outer
		std::function<BPBP(BPBP&&,BPBP&&)> merger = [this](BPBP &&x, BPBP &&y){
			return BPBP_opers.add(std::move(x),std::move(y)); };
		return swapping(R2L_gen[i], merger); };
	std::fstream R2L_tablefile(R2L_filename, std::ios::out | std::ios::binary);
	mon_index.substitution_table(R2L_gens, R2L_tablefile, &BPBP_opers);
	R2L_tablefile.close();
	
	//load the etaL table
	load_etaL(etaL_filename);
	//compute the delta table
	std::function<vectors<matrix_index, BPBP>(const BPBPBP&)> rule = [this](const BPBPBP& srf){
		//change from left unit notation to right unit notation
		std::function< std::pair<matrix_index,BPBP>(exponent,const BPBP&)> cf = [this](exponent e, const BPBP &w){
			return std::make_pair(mon_index.mon_index[e],std::move(etaL(w))); };
		 return BPBPBP_opers.termwise_operation(cf,srf); };
	std::function<void(const vectors<matrix_index,BPBP>&, std::iostream&)> outputer = [this](const vectors<matrix_index,BPBP>& x, std::iostream& writer){
		BPBPMod_opers.save(x, writer); };
	std::function<BPBPBP(int)> delta_gens = [&delta_gen](int i){
		return delta_gen[i]; };
	std::fstream delta_tablefile(delta_filename, std::ios::out | std::ios::binary);
	mon_index.substitution_table(delta_gens, delta_tablefile, &BPBPBP_opers, rule, outputer);
	delta_tablefile.close();
}

//lift elements in F2 to BP
BP BP_Op::lift(F2 x){
	return monomial(0,Z2_oper->lift(x));
}

//lift vectors over F2 to vectors over BP
vectors<matrix_index,BP> BP_Op::lift(const vectors<matrix_index,Fp>&v){
	vectors<matrix_index, BP> result;
	for(auto tm : v.dataArray)
		result.push({tm.ind, lift(tm.coeficient)});
	return result;
}

//the degree of comodules
template<>
int FreeBPCoMod::underlyingDeg(int i){
	return i; }
	
//add degrees
template<>
int FreeBPCoMod::add_degree(int a, int b){
	return a+b; }

//the element v1
BP BP_Op::v1(){
	BP v1 = this->singleton(1);
	return v1;
}
	
//return the element h1 = (etaR(v1) - etaL(v1))/p
BPBP BP_Op::h1(){
	if(mon_index.max_degree<=1){
		std::cerr << "out of range for h1";
		return BPBP_opers.zero();
	}
	auto dv1 = BPBP_opers.add(BPBP_opers.minus(etaL(v1())), etaR(v1()));
	BPBP h1 = divide_power_p(dv1,1);
	std::cout << "h1=" << BPBP_opers.output(h1) << "\n";
	return h1;
}

//return the element h2
BPBP BP_Op::h2(){
	if(mon_index.max_degree<=2){
		std::cerr << "out of range for h2";
		return BPBP_opers.zero();
	}
	//v1^2
	BP v12 = multiply(v1(),v1());
	//etaR(v1^2)-etaL(v1^2)
	auto dv12 = BPBP_opers.add(BPBP_opers.minus(etaL(v12)), etaR(v12));
	//(etaR(v1^2)-etaL(v1^2))/4
	BPBP h2 = divide_power_p(dv12,2);
	std::cout << "h2=" << BPBP_opers.output(h2) << "\n";
	return h2;
}

//return the element h3
BPBP BP_Op::h3(){
	if(mon_index.max_degree<=4){
		std::cerr << "out of range for h3";
		return BPBP_opers.zero();
	}
	//v1^4
	BP v14 = monomial(vars(1)*4);
	//8v2
	BP vv2 = monomial(vars(2),Z2_oper->unit(8));
	//v1^4 + 8v1v2
	v14 = add(v14, multiply(vv2,v1()));
	//etaR - etaL
	auto dv14 = BPBP_opers.add(BPBP_opers.minus(etaL(v14)), etaR(v14));
	//h3 = d(v1^4 + 8v1v2)/16
	BPBP h3 = divide_power_p(dv14,4);
	std::cout << "h3=" << BPBP_opers.output(h3) << "\n";
	return h3;
}

//return the top thetas on the Moore spectrum
std::vector<BPBP> BP_Op::thetas(){
	std::vector<BPBP> theta(10);
	if(mon_index.max_degree<=6){
		std::cerr << "out of range for theta2";
		return theta;
	}
	//v2
	BP v2 = monomial(vars(2));
	std::cout << "v2=" << output(v2) << "\n";
	//v2^2
	BP v22 = multiply(v2,v2);
	//d(v2^2)
	BPBP dv22 = BPBP_opers.add(BPBP_opers.minus(etaL(v22)), etaR(v22));
	//theta2 = d(v2^2)/v1^2
	BPBP theta2 = divide_v1(dv22,2);
	std::cout << "theta2=" << BPBP_opers.output(theta2) << "\n";
	theta[2] = theta2;
	
	if(mon_index.max_degree<=12){
		std::cerr << "out of range for theta3";
		return theta;
	}
	//v2^4
	BP v24 = multiply(v22,v22);
	//d(v2^4 )
	BPBP dv24 = BPBP_opers.add(BPBP_opers.minus(etaL(v24)), etaR(v24));
	//theta3 = d(v2^4 )/v1^4
	BPBP theta3 = divide_v1(dv24,4);
	std::cout << "theta3=" << BPBP_opers.output(theta3) << "\n";
	theta[3] = theta3;
	
	if(mon_index.max_degree<=24){
		std::cerr << "out of range for theta4";
		return theta;
	}
	//v2^8
	BP v28 = multiply(v24,v24);
	//d(v2^8)
	auto dv28 = BPBP_opers.add(BPBP_opers.minus(etaL(v28)), etaR(v28));
	//theta4 = d(v2^8 )/v1^8
	auto theta4 = divide_v1(dv28,8);
	std::cout << "theta4=" << BPBP_opers.output(theta4) << "\n";
	theta[4] = theta4;
	
	if(mon_index.max_degree<=48){
		std::cerr << "out of range for theta5";
		return theta;
	}
	//v2^16
	BP v216 = multiply(v28,v28);
	//d(v2^16)
	BPBP dv216 = BPBP_opers.add(BPBP_opers.minus(etaL(v216)), etaR(v216));
	//theta5 = d(v2^16)/v1^16
	auto theta5 = divide_v1(dv216,16);
	std::cout << "theta5=" << BPBP_opers.output(theta5) << "\n";
	theta[5] = theta5;
	return theta;
}

//divide by p^n
BP BP_Op::divide_power_p(const BP& x, int n){
	std::function<Z2(const Z2&)> rl = [this,n](const Z2 &r){
		return Z2_oper->divide(r,n); };
	return this->termwise_operation(rl,x);
}

//divide by p^n
BPBP BP_Op::divide_power_p(const BPBP& a, int n) {
	std::function<BP(const BP&)> rl = [this,n](const BP& r){ 
		return divide_power_p(r,n);};
	return BPBP_opers.termwise_operation(rl,a);
}

//divide by v1^n
BPBP BP_Op::divide_v1(const BPBP& x,int n){
	//mod 2 reduction
	BPBP x1;
	for(auto tm : x.dataArray){
		BP a1;
		for(auto am : tm.coeficient.dataArray) 
			if(am.coeficient%2!=0){
				a1.push({am.ind, 1});
			}
		if(!isZero(a1))
			x1.push({tm.ind,a1});
	}
	
	//divide by v1^n
	std::function<exponent(exponent)> rl = [this,n](exponent e) { 
		auto nt = unpack(e);
		if(nt[0]<n) std::cerr << "not v1-divisible";
		nt[0] -= n;
		return pack(nt.data());
	};
	return BPBP_opers.re_index(rl,x1);
}
