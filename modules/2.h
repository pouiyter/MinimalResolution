//data structure for vectors with entries in a ring R
template<typename index, typename R>
class vectors{
public:
	//the terms in a vector, including an index and a coeficient
	class term { 
	public:
		//the index of the term
		index ind;
		//the coeficient of the term
		R coeficient; 
		
		//the default constructor
		term(){}
		
		//constructor with given index and coeficient
		term(index i, const R &c) {
			ind = i; coeficient = c; }
			
		//constructor which is destructive for the inputs
		term(index i, R &&c) {
			ind=i; coeficient=std::move(c); }
				
		//comparison operator compares the index of the terms
		bool operator<(const term &y) { return ind<y.ind; }	
	};
	
	//default constructor
	vectors(){}
	
	//the data for the vector, which is an array of terms
	std::vector<term> dataArray; 
	
	//clear the contents in a vector
	inline void clear(){
		dataArray.clear(); }
		
	//the size of the vector
	inline unsigned size() const{
		return dataArray.size();}
			
	//add a new term in the back
	inline void push(const term &tm){
		dataArray.push_back(tm); } 
				
	//add a new term, destructive in the argument
	inline void push(term &&tm){
		dataArray.push_back(std::move(tm)); }
					
	//sort the terms in an unordered array
	virtual void sort(){
		std::sort(dataArray.begin(), dataArray.end()); }
						
	//move the data from an array of terms
	void get_data(vectors<index,R> &&y){ 
		dataArray = std::move(y.dataArray); }
							
	//direct sum of two vectors
	virtual void direct_sum(const vectors<index,R> &y, unsigned rank){
		for(auto &tm:y.dataArray){
			term newTerm(tm.ind + rank, tm.coeficient);
			push(std::move(newTerm));
		}
	}
							
	//the inverse operation of direct sum: split a vector into two vectors, with given splitting point
	vectors<index,R> un_direct_sum(unsigned rank){
		unsigned i;
		//find the term at which to split
		for(i=0; i<size(); ++i)
			if(dataArray[i].ind>=rank) break;
									
		//construct the new vector with the tail
		vectors<index,R> result;
		for(unsigned k=i; k<size(); ++k){
			term newTerm(dataArray[k].ind - rank, std::move(dataArray[k].coeficient));
			result.push(std::move(newTerm));
		}
		//truncate the tail from the old vector
		dataArray.resize(i);
		return result;
	}
	
	//change into dense vector, 0 must transform to zeo
	std::vector<R> toDense(unsigned rank){
		std::vector<R> result(rank);
		if(size()==0) return result;
		if(dataArray[size()-1].ind >= rank)
			std::cerr << "tout:" << dataArray[size()-1].ind << "/" << rank << std::flush;
#pragma omp parallel for schedule(dynamic)
		for(unsigned i=0; i<size(); ++i)
			result[(unsigned) dataArray[i].ind] = dataArray[i].coeficient;
		return result;
	}
	
	//change back from dense vector
	void deDense(std::vector<R> const& y, RingOp<R> *ringop, unsigned start){
		clear();
		for(unsigned i = start; i<y.size(); ++i)
			if(!ringop->isZero(y[i])){
				term newTerm((index) i, y[i]);
				push(newTerm);
			}
	}
	
	//add to a dense vector
	void add2Dense(R c, std::vector<R> &y, RingOp<R> *ringop){
		if(size()==0) return;
		if(dataArray[size()-1].ind >= y.size())
			std::cerr << "aout:" << dataArray[size()-1].ind << "/" << y.size() << std::flush;
#pragma omp parallel for schedule(dynamic)
		for(unsigned i=0; i<size(); ++i){
			R nf = ringop->multiply(c, dataArray[i].coeficient);
			y[(unsigned) dataArray[i].ind] = ringop->add(y[(unsigned) dataArray[i].ind], nf);
		}
	}
};
