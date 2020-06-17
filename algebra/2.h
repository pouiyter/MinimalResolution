// template for the operations on abelian groups
template<typename A>
class AbGroupOp{
public:
	//addition
	virtual A add(const A&, const A&)=0;
	
	//the zero element
	virtual A zero()=0;
	
	//if an element is zero
	virtual bool isZero(const A&)=0;
	
	//negation
	virtual A minus(const A&)=0;
	
	//addition which is destructive for its arguments
	virtual A add(A&&,A&&)=0;
	
	//output as a string
	virtual string output(A)=0;
	
	//save the data to a stream
	virtual void save(const A&, std::iostream&)=0;
	
	//load from a stream
	virtual A load(std::iostream&)=0;
};
