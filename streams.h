//streams.h
#pragma once
#include<functional>
#include<mutex>
#include<iostream>
#include<fstream>
#include<string>

using std::string;

//streams with locks for concurrent accessing
class con_streams{
protected:
	std::mutex locker;
	std::iostream *steam;
	
public:
	//constructor
	con_streams(std::iostream*);
	//clear the contents to free resouces
	virtual void fclear();
	//the destructor
	~con_streams();
	
	//access the stream
	template<typename return_type>
	return_type access(std::function<return_type(std::iostream&)> actor) {
		std::lock_guard<std::mutex> lk(locker);
		return actor(steam);
	}
	
	//read from a given position
	template<typename return_type>
	return_type read(std::function<return_type(std::iostream&)> reader, std::ios::streampos pos){
		std::lock_guard<std::mutex> lk(locker);
		steam->seekg(pos);
		return reader(*steam);
	}
	
	//write and return the position
	std::ios::streampos write(std::function<void(std::iostream&)> writer);
};

//file streams
class con_fstreams : public con_streams{
	//the file for the data
	std::fstream file;
	//the filename
	string filename;
public:
	//constructor
	con_fstreams(string);
	
	//clear the contents to free resouces
	void fclear();
};
