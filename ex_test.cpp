#include"ex_index.h"
#include"ctau_steenrod.h"

main(){
	monomial_index mi(142);
 	mi.init_mon_array();
	for(auto e:mi.mon_array)
		std::cout << total_deg(e) << "\t" << output(e) << ":\t" << "\n";
}
