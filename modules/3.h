using std::to_string;
//template for operations on a R-module
template<typename index, typename R>
class ModuleOp : virtual public AbGroupOp<vectors<index,R>>{
public:
	//operations on the base ring
	RingOp<R> *ringOper;
	//function to output a term
	std::function<string(index,R)> output_term;
	
	//the constructor, initialize ringOper and output_term
	ModuleOp(RingOp<R> *ring_op) {
		ringOper = ring_op;
		output_term = [ring_op](index e, R r){
			return ring_op->output(r) + "[" + to_string(e) + "]"; };
	}
	
	//the zero element
	virtual vectors<index,R> zero(){ 
		return vectors<index,R>(); }
		
	//tell if a vector is zero
	virtual bool isZero(const vectors<index,R> &x){ 
		return x.dataArray.empty(); }
			
	//addition
	virtual vectors<index,R> add(const vectors<index,R> &x, const vectors<index,R> &y);
			
	//destructive addition
	virtual vectors<index,R> add(vectors<index,R> &&x, vectors<index,R> &&y);
			
	//the negation
	virtual vectors<index,R> minus(const vectors<index,R> &x) { 
		vectors<index,R> result;
		result.dataArray.resize(x.size());
		for(unsigned i=0; i<x.size(); ++i){
			result.dataArray[i].ind = x.dataArray[i].ind;
			result.dataArray[i].coeficient = ringOper->minus(x.dataArray[i].coeficient);
		}
		return result;
	}
			
			
	//scalor multiplication with R
	virtual vectors<index,R> scalor_mult(R const &r, const vectors<index,R> &x){
		if(ringOper->isZero(r))
			return zero();
				
		vectors<index,R> res;
		res.dataArray.resize(x.size());
				
		#pragma omp parallel for schedule(dynamic)
		for(unsigned i=0; i<x.size(); ++i){
			res.dataArray[i].ind = x.dataArray[i].ind;
			res.dataArray[i].coeficient = ringOper->multiply(x.dataArray[i].coeficient, r);
		}
		return res;
	}	
			
	//output the vector
	virtual string output(vectors<index,R> x);
			
	//output vectors with the format a[n-i], where n is a specified name
	virtual string output(vectors<index,R> x, string name);
			
	//save the vector
	virtual void save(vectors<index,R> const &x, std::iostream& writer);
			
	//load the vector
	virtual vectors<index,R> load(std::iostream& reader);
			
	//construct a vector with a single entry
	virtual vectors<index,R> singleton(index e, const R &r){
		vectors<index,R> result;
		result.push(typename vectors<index,R>::term(e,r));
		return result;
	}
			
	//singleton with coeficient the unit in R
	virtual vectors<index,R> singleton(index e) {
		return singleton(e,ringOper->unit(1));
	}
			
	//reindex a vector
	template<typename index2>
	vectors<index2,R> re_index(std::function<index2(index)> rule, vectors<index,R> const &x);
			
	//desctructive reindex
	template<typename index2>
	vectors<index2,R> re_index(std::function<index2(index)> rule, vectors<index,R> &&x);
			
	//reindex a vector, filtering out those unwanted terms, the rule produces the invalid index to indicate an unwanted term
	template<typename index2>
	vectors<index2,R> 	filtered_reindex(std::function<index2(index)> rule, vectors<index,R> const &x, index2 invalid = (index2)-1);
			
	//filter a vector
	virtual vectors<index,R> filter(std::function<bool(index)> rule, vectors<index,R> const &x);
			
	//filter when the rule involves the coeficient of a term
	virtual vectors<index,R> filter(std::function<bool(index,const R&)> rule, vectors<index,R> const &x);
			
	//apply an operation termwise
	template<typename index2>
	vectors<index2,R> termwise_operation(std::function<std::pair<index2,R>(index,const R&)> rule, vectors<index,R> const &x);
			
	//apply an operation on each coeficient
	template<typename ring2>
	vectors<index,ring2> termwise_operation(std::function<ring2(const R&)> rule, vectors<index,R> const &x);
			
	//search the coeficient with index n in a range
	virtual R component(index n, vectors<index,R> const &x, unsigned startPos, unsigned endPos);
			
	//search for term with index n
	virtual R component(index n, vectors<index,R> const &x);
			
	//the term next to the term with index n, in a range
	virtual unsigned next_non_trvial_entry(index n, vectors<index,R> const &x, unsigned startPos, unsigned endPos);
			
	//the index next to the term with index n
	virtual index next_non_trvial_entry(index n, vectors<index,R> const &x);
			
	//find the first invertible term
	virtual unsigned first_invertible(vectors<index,R> const &x);
			
	//find the first invertible index
	index first_invertible_index(vectors<index,R> const &x);
			
	//the index of the last entry
	virtual index last_entry(vectors<index,R> const &x);
	
	//the sum of an array of vectors
	vectors<index,R> sum(std::function<vectors<index,R>(int)> summands, int start, int end);
};
