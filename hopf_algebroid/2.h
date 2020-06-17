//template for comodules over a Hopf algebroid, assumed to be free over the base ring with a specified base
template<typename algebroid, typename degree_type>
class CoModule{
public:
	//the rank over the base ring
	virtual int rank() const=0;
	
	//the co-action on the i-th element of the base
	virtual vectors<matrix_index, algebroid> coaction(int i) const=0;
	
	//the degree of the i-th element of the base
	virtual degree_type degree(int i) const=0;
};

//a generic comodule
template<typename algebroid, typename degree_type>
class comodule_generic : virtual public CoModule<algebroid, degree_type>{
public:
	//the coaction is stored in a matrix
	matrix<algebroid> *coaction_matrix;
	//the underlying module
	modules<degree_type> base_module;
	
	comodule_generic(matrix<algebroid> *coactor) { coaction_matrix = coactor; }
	
	//the rank function 
	int rank() const{ 
		return base_module.rank; }
	
	//the coaction is determined by the coaction matrix
	vectors<matrix_index, algebroid> coaction(int i) const{ 
		return coaction_matrix->find(i); }
		
	//the degree function
	degree_type degree(int i) const{ 
		return base_module.degree[i]; }
};
