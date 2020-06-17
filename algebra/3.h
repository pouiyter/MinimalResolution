template<typename R>
class RingOp : virtual public AbGroupOp<R>{
public:
	//multliplication
	virtual R multiply(const R&, const R&)=0;
	  
	//the unit map from the integers
	virtual R unit(int)=0;
	  
	//if an element is invertible
	virtual bool invertible(const R&)=0;
	  
	//the inverse of an invertible element
	virtual R inverse(const R&)=0;
	  
	//power of an elememnt
	R power(R x, unsigned n){
		if(n==0) return unit(1);
		if(n==1) return x;
		auto y = power(x,n/2);
		y = multiply(y,y);
		if(n%2 == 1) return multiply(y,x);
		else return y;
	}
};
