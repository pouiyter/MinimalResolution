template<typename ring, typename algebroid>
class Hopf_Algebroid{
public:
	//operations on base ring
	RingOp<ring> *ringOper;
	//ring operations on algebroid
	RingOp<algebroid> *algebroidRingOper;
	
	//module operations over base ring
	ModuleOp<matrix_index, ring> *moduleOper;
	//module operations over algebroid
	ModuleOp<matrix_index, algebroid> *algebroidModuleOper;
	
	//left and right units
	virtual algebroid etaL(const ring&)=0;
	virtual algebroid etaR(const ring&)=0;
	 
	//convert an element in the algebroid into a vector over the base ring
	virtual vectors<matrix_index, ring> algebroid2vector(const algebroid&,int shift)=0;
	//the inverse of the previous function
	virtual algebroid vector2algebroid(const vectors<matrix_index, ring>&)=0;
	 
	//the co-multiplication on the i-th base element of the algebroid 
	virtual vectors<matrix_index, algebroid> delta(matrix_index)=0;
	
	//the range of the degrees to be computed
	unsigned maxDeg;
	//how many generators for the algbroid is below a certain degree
	virtual unsigned ranksBelowDeg(unsigned)=0; 
	
	//the file name for the data
	string maps_file_name;
	 
	//convert a vector V into a vector over algebroid using the right unit and then scalor multiply with A
	vectors<matrix_index, algebroid> right_scalor_mult (algebroid const &A, vectors<matrix_index, ring> const &V);
	
	//initialize the structure functions of the algebroid
	template<typename degree_type>
	void init_cofree_data(std::function<degree_type(matrix_index)> cofree_degree);
	 
	//adjoint to the projection to the n-th coordinate, returning the resulting map
	template<typename degree_type>
	cofree_comodule<algebroid,degree_type> adjoint(const CoModule<algebroid,degree_type> *X, int n, matrix<ring> *adjoint_map,int shift);
	 
	//adjoint to projectin on many coordinates, the i-th component
	template<typename degree_type>
	vectors<matrix_index,ring> adjoint(const CoModule<algebroid,degree_type> *X, std::vector<int> const& gens, std::vector<uint32_t> const& pos, int i, int shift = 0);
	
	//the trivial comodule of rank one over the base ring 
	template<typename degree_type>
	void set_to_trivial(comodule_generic<algebroid,degree_type> &X, degree_type deg);
	
	//the quotient comodule
	template<typename degree_type>
	void quotient(const CoModule<algebroid,degree_type> *X, matrix<ring> *quot, std::vector<matrix_index> inverse_ind, comodule_generic<algebroid,degree_type> &result);
	template<typename degree_type>
	void quotient_p(const CoModule<algebroid,degree_type> *X, matrix<ring> *quot, std::vector<matrix_index> inverse_ind, comodule_generic<algebroid,degree_type> &result);
	
	//embed a comodule into a cofree comodule
	template<typename degree_type>
	cofree_comodule<algebroid,degree_type> embed2cofree(const CoModule<algebroid,degree_type> *X, matrix<ring> *inj, curtis_table<ring> *table, std::vector<int> *gens, matrix<ring>*, std::vector<int> *basis_order, std::iostream&);
	
	//one step in constructing a relative injective resolution
	template<typename degree_type>
	cofree_comodule<algebroid,degree_type> resolvor(comodule_generic<algebroid,degree_type> &X, matrix<ring> *inj, matrix<ring> *indj, matrix<ring> *quot, curtis_table<ring> *table, matrix<ring>*, std::vector<int> *gens, std::vector<int> *basis_order, std::iostream&);
	 
	//construct a pre-resolution together with the data for the co-generators
	template<typename degree_type>
	void pre_resolution_tab(comodule_generic<algebroid,degree_type>& resolved, string filename_maps, string filename_generators, int resolution_length, std::vector<curtis_table<ring>*>& result, std::vector<std::vector<int>> &gens, matrix<ring>*, matrix<ring>*, matrix<ring>*, matrix<ring>*, std::string tablename = "table.tmp", std::vector<std::vector<int>> *basis_orders = NULL);

	//embed into a cofree comodule using a model
	template<typename degree_type>
	cofree_comodule<algebroid,degree_type> embed2cofree_modeled(const CoModule<algebroid,degree_type> *X, matrix<ring> *inj, std::vector<int> *gens, matrix<ring>* new_maps);
	template<typename degree_type>
	cofree_comodule<algebroid,degree_type> embed2cofree_modeled(const CoModule<algebroid,degree_type> *X, matrix<ring> *inj, std::vector<int> *gens);
	 
	//resolution with a model 
	template<typename degree_type, typename table_type>
	cofree_comodule<algebroid,degree_type> resolvor_modeled(comodule_generic<algebroid,degree_type> &X, matrix<ring> *inj, matrix<ring> *quot, matrix<ring> *indj, matrix<ring> *new_map, curtis_table<table_type> *table, std::vector<int> *gens, std::vector<int> *nex_gens, std::function<vectors<matrix_index,ring>(vectors<matrix_index,table_type> const&)> transformer, std::iostream &maps_file);
	
	//pre-resolution with a model
	template<typename degree_type, typename table_type>
	void pre_resolution_modeled( comodule_generic<algebroid,degree_type>& resolved, string filename_maps, string filename_generators, int resolution_length, string filename_table, curtis_table<table_type>* table, std::vector<std::vector<int>> &gens,std::function<vectors<matrix_index,ring>(vectors<matrix_index,table_type> const&)> transformer, matrix<ring> *inj, matrix<ring> *qut, matrix<ring> *indj, matrix<ring> *new_map, string back_up_file_name, int res_start=-1);
	
	//combine the generators data into one sigle file
	template<typename degree_type>
	void gens_file_combiner(string filename_generators, int resolution_length, comodule_generic<algebroid,degree_type>&);
	
	//composing short exact sequences into a relolution
	template<typename degree_type>
	void resolution(string filename_maps, string filename_generators, string filename_resolution, int resolution_length, comodule_generic<algebroid,degree_type>& tmp, matrix<ring> *inj, matrix<ring>* qut, matrix<ring>* composed, std::ostream *outfile=NULL);
	
	//construct a multiplication table
	void make_multiplication_table(algebroid const &x, int deg_x, matrix<ring> *result);
	
	//load the data for generators
	static void load_gens(std::vector<std::vector<int>>&,string);
};
