//BPcomplex.cpp
#include"BPcomplex.h"

//BP operations
BP_Op *primitive_data::BPoper;
//Zp module operations
ModuleOp<matrix_index,Z2> *primitive_data::Z2Mod_oper;

//set operators
void primitive_data::set_oper(BP_Op* BPop){
	BPoper = BPop;
	Z2Mod_oper = &BPop->Z2Mod_opers;
}

//compute the complex on primitives, the formula being etaR(v^e)[shift]
void BPComplex::make_prim_map(const primitive_data& source, const primitive_data& target, matrix<BP>* X, matrix<BP> *result, matrix<Z2> *map){
	std::function<vectors<matrix_index,BP>(int)> rows = [&source,X](int i){
		//get the i-th primitive
		auto v = primitive_data::BPoper->monomial(source[i].coeficient);
		auto v1 = primitive_data::BPoper->etaR(v);
		auto w = primitive_data::BPoper->algebroid2vector(v1,source.gen_shift[i]);

		//compute the image of primitives
		return X->maps_to(w);
	};
	result->construct(source.size(),rows);
	
	//convert the map to a map of Zp modules
	std::function<vectors<matrix_index,Z2>(int)> rs = [&result, &target](int k){
		return target.expand(result->find(k)); };
		
	map->construct(source.size(),rs);
}

//add a new entry for the primitives
void primitive_data::add(prim_entry itm, int shift, int deg){
	auto k = size();
	prim_index.emplace(itm, k);
	push_back(itm);
	gen_shift.push_back(shift);
	gen_deg.push_back(deg);
}

//make the primitive data using the list of monimials
primitive_data primitive_data::make_primitives(int deg, int pos, int shift, monomial_index* mon_ind){
	primitive_data result;
	for(auto e: mon_ind->mon_array){
		if (total_deg(e) + deg > mon_ind->max_degree) break;
		result.add({(matrix_index)pos,e},shift,total_deg(e) + deg );
	}
	return result;
}

//combining two datas
void primitive_data::direct_summed(const primitive_data& y){
	for(unsigned i=0; i<y.size(); ++i)
		add(y[i],y.gen_shift[i],y.gen_deg[i]);
}

//make data of primitive from a cofree comodule
primitive_data primitive_data::primitive_maker(cofree_comodule<BPBP,int> &F, monomial_index* mon_ind){
	primitive_data result;
	for(int i=0; i<(int)F.generators.rank; ++i)
		result.direct_summed(make_primitives(F.generators.degree[i], i, F.position_of_gens[i], mon_ind));
	return result;
}

//comparison of entries
bool prim_entry::operator<(const prim_entry y) const{
	if(gen_pos==y.gen_pos) return coeficient<y.coeficient;
	else return gen_pos < y.gen_pos;
}

//expand a verctor over BP into a vector over Zp
vectors<matrix_index, Z2> primitive_data::expand(const vectors<matrix_index,BP>& v) const{
	vectors<matrix_index, Z2> result;

	for(unsigned i=0; i<v.size(); ++i)
		for(unsigned j=0; j<v.dataArray[i].coeficient.size(); ++j){
			//the corresponding element is v^(e)[g], set ns={g,e}
			prim_entry ns = {v.dataArray[i].ind,v.dataArray[i].coeficient.dataArray[j].ind};
			
			//find the place of this entry
			auto n = prim_index.at(ns);
			
			//add this new term to result
			result = Z2Mod_oper->add(result, Z2Mod_oper->singleton(n,v.dataArray[i].coeficient.dataArray[j].coeficient));
		}
	return result;
}

//output primitive entris
string prim_entry::output() const {
	auto cf = coeficient;
	return ::output(cf,"v") + "[" + std::to_string(gen_pos) + "]";
}

//output the primitives
string primitive_data::output(){
	string res;
	for(auto sf : *this)
		res += std::to_string(prim_index[sf]) + ":" + sf.output() + "\n";
	return res;
}

//load the primitives
std::vector<int> BPComplex::load_prim(int resolution_length, std::fstream &gens_file){
	int32_t M_rank;
	std::vector<int> F_rank(resolution_length+1);
	for(int i=0; i<=resolution_length; ++i){
		FreeBPCoMod F;
		gens_file.read((char*)&M_rank,4);
		F.load(gens_file);
		Prims[i] = primitive_data::primitive_maker(F,&primitive_data::BPoper->mon_index);
		F_rank[i]=F.rank();
	}
	return F_rank;
}

//comstruct a complex from a resolution by taking the primitives
void BPComplex::load(int resolution_length, string generator_filename, string maps_filename, matrix<BP>* mp, matrix<BP>* pm, std::function<matrix<Z2>*(int)> map_constr){
	//open the files for the generators and maps of a resolution
	std::fstream gens_file(generator_filename, std::ios::in | std::ios::binary);
	std::fstream maps_file(maps_filename, std::ios::in | std::ios::binary);
	
	//intialize the containers
	Prims.resize(resolution_length+1);
	Maps.resize(resolution_length+1);
	
	//initialize the maps
	for(int i=0;i<=resolution_length;++i){
		Maps[i] = map_constr(i);
		Maps[i]->clear();
	}
		
	std::cout << "constructing primitve data...\n" << std::flush;
	
	//load the data of the generators
	auto F_rank = load_prim(resolution_length,gens_file);
	
	//construct the maps
	for(int i=1;i<=resolution_length; ++i){
		std::cout << i << " " << std::flush;
		//load the matrix in the resolution
		mp->load(maps_file);
		
		//construct the map of primitives
		make_prim_map(Prims[i-1],Prims[i],mp,pm,Maps[i]);
	}
}

//save to file
void BPComplex::save_matrix(string filename){
	std::fstream witer(filename, std::ios::out | std::ios::binary);
	for(unsigned i=1;i<Maps.size();++i)
		Maps[i]->save(witer);
}

//load from file
void BPComplex::load_matrix(int resolution_length, string generator_filename, string matrix_filename){
	//open the files
	std::fstream gens_file(generator_filename, std::ios::in | std::ios::binary);
	std::fstream matrix_file(matrix_filename, std::ios::in | std::ios::binary);
	
	//initialize
	Prims.resize(resolution_length+1);
	for(int i=0;i<=resolution_length;++i)
		Maps[i]->clear();
	
	//load the generators
	load_prim(resolution_length,gens_file);
	
	for(int i=1;i<=resolution_length; ++i)
		Maps[i]->load(matrix_file);
}

unsigned BPComplex::size(){
	return Prims.size(); }

//get the data for the generators
std::vector<FreeBPCoMod> BPComplex::get_generator(int resolution_length, string generator_filename){
	std::fstream gens_file(generator_filename, std::ios::in | std::ios::binary);
	std::vector<FreeBPCoMod> result(resolution_length+1);
	int32_t M_rank;
	for(int i=0; i<=resolution_length; ++i){
		gens_file.read((char*)&M_rank,4);
		result[i].load(gens_file);
	}
	return result;
}
