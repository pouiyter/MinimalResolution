//BPQ.cpp
#include "BPQ.h"

BPQ_Op::BPQ_Op(Qp_Op *Qpoper, Qp_int *Qpintop, int Mvar) : ModuleOp<exponent, Qp>::ModuleOp(Qpoper), PolynomialOp_Para<Qp>(Qpoper), BPBPQ_oper(this), BPBPBPQ_oper(&BPBPQ_oper), BPint_oper(Qpintop), BPBPint_oper(&BPint_oper), BPBPBPint_oper(&BPBPint_oper){
	maxVar = Mvar;
	
	//output function
	this->output_term = [Qpoper](exponent e, Qp r){ 
		return Qpoper->output(r) + ::output(e,"v"); };
	//initialize Qp operations
	Qp_oper = Qpoper;
	Qpint_oper = Qpintop;
    
	PolynomialOp_Para<BPQ> P_BPQ_oper(this);
	PolynomialOp_Para<BPBPBPQ> P_BPBPBPQ_oper(&BPBPBPQ_oper);

	//initialize li table
	std::cout << "constructing li tables\n" << std::flush;
	li_table.resize(maxVar+1);
	li_table[0] = this->unit(1);
	for(int n=1; n<=maxVar; ++n){
		//compute l_n = (1/p)v_n + \sum_k (1/p)v_{n-k}^{p^k} l_k
		polynomial<BPQ> li_rule = P_BPQ_oper.monomial(0, this->monomial(singleVar(n,1),Qp_oper->construct(1,-1)));
		for(int k=1; k<n; ++k)
			li_rule = P_BPQ_oper.add(std::move(li_rule), std::move(P_BPQ_oper.monomial(singleVar(k,1), this->monomial(singleVar(n-k,Qp_oper->power_p_int(k)), Qp_oper->construct(1,-1)))));
		
		li_table[n] = substitute(li_table, li_rule, this);
	}
	
	//initialize vn table
	std::cout << "constructing vn tables\n" << std::flush;
	vn_table.resize(maxVar+1);
	vn_table[0] = this->unit(1);
	for(int n=1; n<=maxVar; ++n){
		//compute vn = pl_n + \sum (-1)l_k v_{n-k}^{p^k}
		polynomial<BPQ> vn_rule = P_BPQ_oper.monomial(0,this->monomial(singleVar(n,1), Qp_oper->construct(1,1)));
		for(int k=1; k<n; ++k)
			vn_rule = P_BPQ_oper.add(std::move(vn_rule),  std::move(P_BPQ_oper.monomial(singleVar(n-k,Qp_oper->power_p_int(k)), this->monomial(singleVar(k,1), Qp_oper->construct(-1,0)))));
		
		vn_table[n] = substitute(vn_table, vn_rule, this);
	}
	

	//initialize etaR table (in terms of l_i)
	etaR_table.resize(maxVar+1);
	etaR_table[0] = BPBPQ_oper.unit(1);
	for(int n=1; n<=maxVar; ++n){
		//compute etaR(l_n) = [t_n] + l_n + \sum l_k[t_{n-k}^{p^k}}
		BPBPQ etaR_rule = BPBPQ_oper.add(BPBPQ_oper.construct(0,0,n,1,Qp_oper->unit(1)), BPBPQ_oper.construct(n,1,0,0,Qp_oper->unit(1)));
		for(int k=1; k<n; ++k)
			etaR_rule = BPBPQ_oper.add(std::move(etaR_rule), BPBPQ_oper.construct(k,1,n-k,Qp_oper->power_p_int(k),Qp_oper->unit(1)));
		
		etaR_table[n] = etaR_rule;
	}
	
	//initialize delta table (in terms of l_i)
	delta_table.resize(maxVar+1);
	delta_table[0] = BPBPBPQ_oper.unit(1);
	for(int n=1; n<=maxVar; ++n){
		//compute D_{s,j} = l_s ([1|t_{j-s}^{p^s}] + [t_{j-s}^{p^s}|1] + \sum [t_k^{p^s}|t_{j-s-k}^{p^{k+s}}])
		auto del_term = [this](int s, int j){
			auto result = BPBPBPQ_oper.add(BPBPBPQ_oper.construct(s,1,0,0,j-s,Qp_oper->power_p_int(s), Qp_oper->unit(1)),  BPBPBPQ_oper.construct(s,1,j-s,Qp_oper->power_p_int(s),0,0, Qp_oper->unit(1)));
			for(int k=1; k<j-s; ++k)
				result = BPBPBPQ_oper.add(std::move(result), BPBPBPQ_oper.construct(s,1, k,Qp_oper->power_p_int(s),j-s-k,Qp_oper->power_p_int(k+s), Qp_oper->unit(1)));
			return result;
		};
		//compute C_j = [1|t_j] + [t_j|1] + \sum [t_k|t_{j-k}^{p^k}] + \sum D_{s,j}
		auto del_cons = [this,del_term](int j){
			auto result = BPBPBPQ_oper.add(BPBPBPQ_oper.construct(0,0,0,0,j,1,Qp_oper->unit(1)), BPBPBPQ_oper.construct(0,0,j,1,0,0,Qp_oper->unit(1)));
			for(int k=1;k<j;++k)
				result = BPBPBPQ_oper.add(std::move(result), BPBPBPQ_oper.construct(0,0,k,1,j-k,Qp_oper->power_p_int(k),Qp_oper->unit(1)));
			for(int s=1;s<j;++s)
				result = BPBPBPQ_oper.add(std::move(result), del_term(s,j));
			return result;
		};
		
		//compute delta(t_n) =  = C_n + \sum (-1) l_k delta(t_{n-k})^{p^k} 
		polynomial<BPBPBPQ> delta_rule = P_BPBPBPQ_oper.monomial(0,del_cons(n));
		for(int k=1; k<n; ++k)
			delta_rule = P_BPBPBPQ_oper.add(std::move(delta_rule),  P_BPBPBPQ_oper.monomial(singleVar(n-k,Qp_oper->power_p_int(k)), BPBPBPQ_oper.construct(k,1,0,0,0,0,Qp_oper->unit(-1))));
		
		delta_table[n] = substitute(delta_table, delta_rule, &BPBPBPQ_oper);
	}
}

//constructor
BPBPQ_Op::BPBPQ_Op(BPQ_Op* op) : ModuleOp<exponent, BPQ>::ModuleOp(op), PolynomialOp_Para<BPQ>(op){
	this->output_term = [this](exponent e, BPQ r){
		string result = BPQ_opers->Qp_oper->output(r.dataArray[0].coeficient) + ::output(e,"v") + "[" + ::output(r.dataArray[0].ind,"t") + "]"; 
		for(unsigned i=1; i<r.dataArray.size(); ++i)
			result += "+" + BPQ_opers->Qp_oper->output(r.dataArray[i].coeficient) + ::output(e,"v") + "[" + ::output(r.dataArray[i].ind,"t") + "]";
		return result;
	};
	
	BPQ_opers = op; 
}

//constructor
BPBPBPQ_Op::BPBPBPQ_Op(BPBPQ_Op* op) : ModuleOp<exponent, BPBPQ>::ModuleOp(op), PolynomialOp_Para<BPBPQ>(op){
	BPBPQ_opers = op; }

//v_n^s[t_k^i] ,vn goes in the outer level
BPBPQ BPBPQ_Op::construct(int n, int s, int k, int i, Qp coef){
	return this->monomial(singleVar(n,s), BPQ_opers->monomial(singleVar(k,i),coef));
}

//v_n^s[t_{k1}^{i1}|t_{k2}^{i2}] vn goes in the outer, left ti goes in the inner, right ti goes in the middle
BPBPBPQ BPBPBPQ_Op::construct(int n, int s, int k1, int i1, int k2, int i2, Qp coef){
	return this->monomial(singleVar(n,s), BPBPQ_opers->construct(k2,i2,k1,i1,coef));
}

//transform li expression to vn expression
BPQ BPQ_Op::l2v(const BPQ &x){
	std::function<BPQ(int n)> values = [this](int n){
		return this->li_table[n]; };
	std::function<BPQ(Qp)> crule = [this](Qp y){
		return this->constant(y); };
	
    return substitute(values,crule,x,this);
}

//transform vn expression to li expression
BPQ BPQ_Op::v2l(const BPQ &x){
	std::function<BPQ(int n)> values = [this](int n){
		return vn_table[n]; };
	std::function<BPQ(Qp)> crule = [this](Qp y){
		return this->constant(y); };
		
    return substitute(values,crule,x,this);
}

//left unit
BPBPQ BPQ_Op::etaL(const BPQ &x){
	std::function<BPQ(const Qp&)> rl = [this](const Qp &a){
		return this->constant(a); };
		return this->termwise_operation(rl,x);
}

//transformation on BPBPQ
BPBPQ BPQ_Op::l2v(const BPBPQ &x){
	std::function<BPBPQ(int n)> values = [this](int n){
		auto va = li_table[n]; 
		return etaL(va);
	};
	std::function<BPBPQ(BPQ)> crule = [this](BPQ y){
		return BPBPQ_oper.constant(y); };
		
    return substitute(values,crule,x,&BPBPQ_oper);
}

//transformation on BPBPQ
BPBPQ BPQ_Op::v2l(const BPBPQ &x){
	std::function<BPBPQ(int n)> values = [this](int n){
		auto va = vn_table[n]; 
		return etaL(va);
	};
	std::function<BPBPQ(BPQ)> crule = [this](BPQ y){
		return BPBPQ_oper.constant(y); };
		
	return substitute(values,crule,x,&BPBPQ_oper);
}

//transform on BPBPBPQ
BPBPBPQ BPQ_Op::l2v(const BPBPBPQ &x){
	std::function<BPBPQ(const Qp&)> rl = [this](const Qp &a){
		return BPBPQ_oper.constant(this->constant(a)); };
	std::function<BPBPBPQ(int n)> values = [this,rl](int n){
		auto va = li_table[n]; 
		return this->termwise_operation(rl,va);
	};
	std::function<BPBPBPQ(BPBPQ)> crule = [this](BPBPQ y){
		return BPBPBPQ_oper.constant(y); };
    
	return substitute(values,crule,x,&BPBPBPQ_oper);
}

//right unit in terms of li
BPBPQ BPQ_Op::etaR_l(const BPQ &x){
	std::function<BPBPQ(int n)> values = [this](int n){
		return etaR_table[n]; };
	std::function<BPBPQ(Qp)> crule = [this](Qp y){
			return BPBPQ_oper.constant(this->constant(y)); };
	
	return substitute(values,crule,x,&BPBPQ_oper);
}
	
//right unit in terms of vn
BPBPQ BPQ_Op::etaR_v(const BPQ &x){
	auto ls = v2l(x);
	auto etaRl = etaR_l(ls);
	return l2v(etaRl);
}

//delta in terms of li
BPBPBPQ BPQ_Op::delta_l(exponent e){
	std::function<BPBPBPQ(int n)> values = [this](int n){
		return delta_table[n]; };
	return substitute(values,e,&BPBPBPQ_oper);
}

//delta in terms of vn
BPBPBPQ BPQ_Op::delta_v(exponent e){
    auto ds = delta_l(e);
    return l2v(ds);
}

//output the formula for li
string BPQ_Op::show_li(){
	string result;
	for(int i=1; i<=maxVar; ++i){
		result += "l" + std::to_string(i) + " = ";
		result += this->output(li_table[i]);
		result += "\n";
	}
	
	return result;
}

//check the consistancy of the formula for vn
string BPQ_Op::show_vn(){
	string result;
	for(int i=1; i<=maxVar; ++i){
		BPQ x = this->monomial(singleVar(i,1));
		auto vn = v2l(x);
		vn = l2v(vn);
		result += "v" + std::to_string(i) + " = ";
		result += this->output(vn);
		result += "\n";
	}
	
	return result;
}

//output the formula for etaR
string BPQ_Op::show_etaR(){
	string result;
	for(int i=1; i<=maxVar; ++i){
		BPQ x = this->monomial(singleVar(i,1));
		result += "etaR(v" + std::to_string(i) + ") = ";
		result += BPBPQ_oper.output(etaR_v(x));
		result += "\n";
	}
	
	return result;
}

//output the etaR formular, in intergral form, vn in the outer, ti in the inner
void BPQ_Op::output_R2L(std::iostream &writer){
	for(int i=1; i<=maxVar; ++i){
		BPQ x = this->monomial(singleVar(i,1));
		BPBPQ y = etaR_v(x);
		BPBPint_oper.save(y, writer);
	}
}

//output the delta formular, in intergral form, vn in the middle, right ti in the outer, left ti in the inner
void BPQ_Op::output_delta(std::iostream &writer){
	for(int i=1; i<=maxVar; ++i){
		BPBPBPQ x = delta_v(singleVar(i,1));
		BPBPBPQ y = swapping(x, &BPBPBPQ_oper);
		BPBPBPint_oper.save(y, writer);
	}
}
