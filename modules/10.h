//class for free graded modules over a ring
template<typename degree_type>
class modules{
public:
	//the rank of the module
	unsigned rank;
	//array of the degrees of the generators
	std::vector<degree_type> degree;
	
	//the default constructor
	modules(){ rank=0; }
	
	//save the module. 4 bytes for the rank, then the data for the degrees.
	void save(std::iostream &writer){
		uint32_t rk = rank;
		writer.write((char*)&rk, 4);
		writer.write((char*)degree.data(), sizeof(degree_type)*rk);
	}
	
	//load a module
	void load(std::iostream &reader){
		uint32_t rk;  
		reader.read((char*)&rk, 4);
		rank = rk;
		degree.resize(rk);
		reader.read((char*)degree.data(), sizeof(degree_type)*rk);
	}
	
	//direct sum with another module
	void direct_sum(const modules<degree_type> &Y){
		rank+=Y.rank;
		degree.insert(degree.end(), Y.degree.begin(), Y.degree.end());
	}
	
	//output the module, with spoecified way to output the degrees
	string output(std::function<string(degree_type)> deg2string){
		string res = "";
		for(auto d: degree) 
			res += deg2string(d) + " ";
		return res;
	}
	
	//default output function
	string output(){ 
		return output([] (degree_type d) {return std::to_string(d); }); 
	}
};
