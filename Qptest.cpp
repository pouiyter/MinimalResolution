#include"Qp.h"

int main(){
	Qp w;
	Q2_Op q2;
	Q2_int qi;
	w = q2.unit(-8);
	
	std::cout << q2.output(w) << "\n";
	std::cout << qi.output(w);
}
