//matrix stored in memeries
template<typename ring>
class matrix_mem : public matrix<ring>{
	//the data for a matrix, consisting of an array of vectors
	std::vector<vectors<matrix_index,ring>> data;
public:
	//default constructor
	matrix_mem(){
		this->rank = 0; }
	
	matrix_mem(string){
		this->rank = 0; }
	
	//clear the contents
	void clear(){
		std::vector<vectors<matrix_index,ring>> emd;
		data.swap(emd);
		data.clear();
		this->rank = 0;
	}
	
	//find the n-th row
	vectors<matrix_index,ring> find(matrix_index n) const{
		return data[n];
	}
	
	//insert a new row
	void insert(matrix_index i, vectors<matrix_index,ring> const& x){
		if(i>=data.size()){
			data.resize(i+1);
			this->rank = i+1;
		}
		data[i] = x;
	}
	
	//set the rank
	void set_rank(unsigned n){
		this->rank = n;
		data.resize(n);
	}
	
	//update all rows
	void update_all(std::function<void(vectors<matrix_index,ring>&,matrix_index)> action){
		for(unsigned i=0; i<data.size(); ++i)
			action(data[i],i);
	}
	
	//row transfermations, subtract the given row to cancel the given column	
	void row_reduction(matrix_index row, matrix_index col);
	
	//Gaussian ellimination, with given rows and colums
	void gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols);
	
	//delete some columns and then do Gaussion
	void del_and_gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols, std::set<int> const &to_del){
		this->del_cols(to_del);
		gaussian(row_cols);
	}
	
	//unify a vector to make a given entry 1
	vectors<matrix_index,ring> unify(matrix_index n, vectors<matrix_index,ring> const &x);
	
	//ellimination a given entry
	void take_away(matrix_index n, vectors<matrix_index,ring> &x, vectors<matrix_index,ring> const &y);
};
