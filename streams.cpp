//streams.cpp
#include"streams.h"

//constructor
con_streams::con_streams(std::iostream* stm){
	steam = stm;
}

//clear the contents to free resouces
void con_streams::fclear(){
	steam->seekp(0, std::ios::beg);
}

//the destructor
con_streams::~con_streams(){
//	fclear();
}

//write and return the position
std::ios::streampos con_streams::write(std::function<void(std::iostream&)> writer){
	std::lock_guard<std::mutex> lk(locker);
	auto pos = steam->tellp();
	writer(*steam);
	return pos;
}

//constructor
con_fstreams::con_fstreams(string fname) : con_streams(&file), file(fname, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc){
	filename = fname;
	if(!file.is_open())
		std::cerr << "fail to open" << filename << std::flush;
}

//clear the contents
void con_fstreams::fclear(){
	std::lock_guard<std::mutex> lk(locker);
	file.close();
	file.open(filename,  std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
	if(!file.is_open())
		std::cerr << "fail to open" << filename << std::flush;
}
