#ifndef SYNTH_HEADER
#define SYNTH_HEADER
#include <functional>
#include <algorithm>
#include <vector>
#include <assert.h>
#include <iostream>
#include <string>
#include <sstream>
// Define P Data structure
template <typename T>
using Comparator = std::function<bool (std::vector<T>&,std::vector<T>&)>;

template <typename T>
class Permutation{
private:
    std::vector<std::vector<T>> d;
    Comparator<T> sortConstraint;
public:
    Permutation(Comparator<T> sortConstraint): sortConstraint(sortConstraint){}
    Permutation(){
       this->sortConstraint = NULL;
    }
    void insert(std::vector<T> tup){
        d.push_back(tup);
	if (sortConstraint != NULL){
	    std::sort(d.begin(),d.end(),sortConstraint);
	}
    }
    int get(std::vector<T> tup){
        auto it = std::find(d.begin(),d.end(),tup);
	if (it == d.end()) {
	    std::stringstream ss;
	    ss << "Permutation::get: Tuple {";

	    for(int j = 0; j  < tup.size(); j++){
	        ss << tup[j] << ",";
	    }
	    ss << "} not found";
	    std::cerr << ss.str();
	    assert(0 && ss.str().c_str());
	}
	return it - d.begin();
    }
    std::string toString(){
	std::stringstream ss;
	for(int i = 0; i < d.size(); i++){
	    ss<< "[" << i << "] => {";
	    for(int j = 0; j  < d[i].size(); j++){
	        ss << d[i][j] << ",";
	    }
	    ss << "}";
	}
	return ss.str();
    }
};

#endif SYNTH_HEADER
