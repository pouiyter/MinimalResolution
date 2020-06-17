//generate an index for the positions of the co-generators
template<typename algebroid, typename degree_type>
std::map<matrix_index,matrix_index> cofree_comodule<algebroid,degree_type>::gen_pos_index(){
	std::map<matrix_index,matrix_index> result;
	for(unsigned i=0; i<position_of_gens.size(); ++i)
		result.emplace(position_of_gens[i],i);
	return result;
}

//save the data of generators
template<typename algebroid, typename degree_type>
void cofree_comodule<algebroid,degree_type>::save(std::iostream &writer){
	int32_t rk = total_rank;
	writer.write((char*)&rk,4);
	generators.save(writer);
	uint32_t n_gens = position_of_gens.size();
	writer.write((char*)&n_gens, 4);
	writer.write((char*)position_of_gens.data(), n_gens*4);
}

//load a cofree comodule
template<typename algebroid, typename degree_type>
void cofree_comodule<algebroid,degree_type>::load(std::iostream &reader){
	int32_t rk; reader.read((char*)&rk, 4);
	total_rank = rk;
	generators.load(reader);
	uint32_t n_gens; reader.read((char*)&n_gens, 4);
	position_of_gens.resize(n_gens);
	reader.read((char*)position_of_gens.data(), n_gens*4);
}

//default generator
template<typename algebroid, typename degree_type>
cofree_comodule<algebroid,degree_type>::cofree_comodule(){ 
	total_rank = 0; }
	

//find the position of a generator
template<typename algebroid, typename degree_type>
unsigned cofree_comodule<algebroid,degree_type>::findPos(unsigned startPos, unsigned endPos, unsigned n) const{
	if(endPos<=startPos) return startPos-1;
	if(endPos-startPos==1){
		if(position_of_gens[startPos]<=n) return startPos;
		else return startPos-1;
	}
	unsigned midPos = (startPos+endPos)/2;
	if(position_of_gens[midPos]==n) return midPos;
	else if(position_of_gens[midPos]>n) return findPos(startPos,midPos,n);
	else return findPos(midPos+1,endPos,n);
}

template<typename algebroid, typename degree_type>
unsigned cofree_comodule<algebroid,degree_type>::findPos(unsigned n) const{ 
	return findPos(0,position_of_gens.size(),n); }

//the rank of the comodule
template<typename algebroid, typename degree_type>
int cofree_comodule<algebroid,degree_type>::rank() const{ 
	return total_rank; }

//the coaction function
template<typename algebroid, typename degree_type>
vectors<matrix_index, algebroid> cofree_comodule<algebroid,degree_type>::coaction(int i) const{
	unsigned pos = findPos(i);
	std::function<matrix_index(matrix_index)> shifting = [pos, this] (matrix_index k){
		return k + this->position_of_gens[pos]; };
	return modOpers->re_index(shifting, cofree_coaction(i-position_of_gens[pos]));
}

//the degree of the generators
template<typename algebroid, typename degree_type>
degree_type cofree_comodule<algebroid,degree_type>::degree(int i) const{
	unsigned pos = findPos(i);
	return add_degree(cofree_degree(i-position_of_gens[pos]), generators.degree[pos]);
}

//direct sum with another cofree comodule
template<typename algebroid, typename degree_type>
void cofree_comodule<algebroid,degree_type>::direct_sum( cofree_comodule<algebroid,degree_type> const &Y){
	generators.direct_sum(Y.generators);
	for(auto p: Y.position_of_gens)
		position_of_gens.push_back(p+total_rank);
	total_rank+=Y.total_rank;
}
		
// find if n is a cogenerators
template<typename algebroid, typename degree_type>
int cofree_comodule<algebroid,degree_type>::find_index(int n){
	int h = findPos(n);
	if((int)position_of_gens[h] == n) return h;
		else  return invalid_pos;
}
