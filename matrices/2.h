//template class for matrices. matrices are regarded as an array of row veectors
template<typename R>
class matrix{
public:
	//the number of ranks
	unsigned rank;

	//operations on R-modules
	static ModuleOp<matrix_index,R> *moduleOper;
	
	//clear the contents of a matrix
	virtual void clear()=0;
	
	//find the n-th row 
	virtual vectors<matrix_index,R> find(matrix_index n) const=0;
	
	//insert a new row to the matrix
	virtual void insert(matrix_index i, vectors<matrix_index,R> const& x)=0;
	
	//set the rank of the matrix
	virtual void set_rank(unsigned)=0;
	
	//update all rows
	virtual void update_all(std::function<void(vectors<matrix_index,R>&,matrix_index)>)=0;
	
	//filter the columns
	void filter(std::function<bool(matrix_index)> rule){
		std::function<void(vectors<matrix_index,R>&,matrix_index)> action = [this, &rule] (vectors<matrix_index,R>& x, matrix_index){ 
			x = moduleOper->filter(rule,x); };
		update_all(action);
	}
	
	//delete some columns
	void del_cols(std::set<int> const &to_del){
		std::function<bool(matrix_index)> rule = [&to_del](matrix_index n){ 
			return to_del.count(n) == 0; };
		filter(rule);
	}
	
	//construct a matrix with specified rows
	void construct(int rk, std::function<vectors<matrix_index,R>(int)> rows){
		clear();
		set_rank(rk);
		for(unsigned i=0; i<rank; ++i)
			insert(i,rows(i));
	}
	
	//construct with parallel algorithm
	virtual void construct_parallel(int rk, std::function<vectors<matrix_index,R>(int)> rows){
		clear();
		set_rank(rk);
		#pragma omp parallel for schedule(dynamic)
		for(unsigned i=0; i<rank; ++i)
			insert(i,rows(i));
	}
	
	//copy from another matrix
	void construct(int rk, matrix<R> *M){
		set_rank(rk);
		std::function<vectors<matrix_index,R>(int)> rows = [M](int i){
			return M->find(i); };
			construct(rank, rows);
	}
	
	//copy from another matrix
	void construct(matrix<R> *M){
		construct(M->rank, M);
	}
	
	//set a matrix to the unit matrix
	void set2unit(int rk);
	
	//set a matrix to zero
	void set2zero(int rk);
	
	bool equal(matrix<R>&);
	
	//output a matrix
	string output();
	
	//output with format [a-i]->[b-j]+...
	string output(string name0, string name1);
	
	//save a matrix. 4 bytes of rank, followed by the rows
	void save(std::iostream &writer);
	
	//load a matrix
	void load(std::iostream &reader);
	
	//load a matrix with known rank
	void load(std::iostream &reader, int);
	
	//load and modify
	void load_modify(std::iostream&, std::function<vectors<matrix_index,R>(const vectors<matrix_index,R>&)>);
	
	//right multiplication on a row vector v
	vectors<matrix_index,R> maps_to(const vectors<matrix_index,R> &v);
	
	//compose this with rows, with rank of rows equal to res_rank
	void compose(std::function<vectors<matrix_index,R>(int)> rows, int result_rank, matrix<R> *result);
	
	//compose with Y
	void compose(const matrix<R> *Y, matrix<R> *result);
	
	//compose with another matrix with different type of coeficients
	template<typename coef_type,typename module2>
	module2 maps_to(vectors<matrix_index,coef_type> const &x, std::function<module2(const coef_type&, const vectors<matrix_index,R>&)> scalor_mult, AbGroupOp<module2> *modop2);
	template<typename coef_type,typename module2>
	module2 maps_to_p(vectors<matrix_index,coef_type> const &x, std::function<module2(const coef_type&, const vectors<matrix_index,R>&)> scalor_mult, AbGroupOp<module2> *modop2);
	
	//make a quotient matrix. This should be allready in an echelon form
	void make_quotient(std::vector<matrix_index> const &inverse_ind, std::map<matrix_index,matrix_index> const &quot_ind, std::vector<matrix_index> const &inj_index, matrix<R> *result);

	//operations on each rows
	void row_operation(std::function<vectors<matrix_index,R>(vectors<matrix_index,R> const&,int)> row_rule, matrix<R> *result);
	
	//filter each row vector
	void filter(std::function<bool(matrix_index)> column_rule, matrix<R> *result);
	
	//filter each row vector
	void filter(std::function<bool(matrix_index,const R&)> column_rule, matrix<R> *result);
	
	//Gaussian ellimination
	virtual void gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols)=0;
	
	//delete some columns and then do Gaussion
	virtual void del_and_gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols, std::set<int> const &to_del)=0;
	
	//merge two matrices
	void merge(matrix<R> const *Y){
		std::function<void(vectors<matrix_index,R>&,matrix_index)> row_summer = [Y,this](vectors<matrix_index,R>& ro, matrix_index i){ 
			ro.direct_sum(Y->find(i),0); };
		this->update_all(row_summer);
	}

	//return the leading term index for a quotient map. inj_index lists the leading indices of an injective map. Return the indices left out in the quotient space and the leading term correspondence for the quotient map
	static std::pair<std::vector<matrix_index>,std::map<matrix_index,matrix_index>> quot_index(std::vector<matrix_index> const  &inj_index, unsigned rank);
	
	
	//take the direct sum with another matrix
	void direct_sum(matrix<R> const &Y, unsigned tar_rank) {
		std::function<void(vectors<matrix_index,R>&,matrix_index)> row_summer = [&Y,tar_rank](vectors<matrix_index,R>& ro, matrix_index i){ 
			ro.direct_sum(Y.find(i), tar_rank);  };
		update_all(row_summer);
	}
};

template<typename R>
ModuleOp<matrix_index,R> *matrix<R>::moduleOper;
