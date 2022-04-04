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
    int tupleSplit = 0;
    Comparator<T> sortConstraint;
public:
    Permutation(Comparator<T> sortConstraint): tupleSplit(tupleSplit),
	sortConstraint(sortConstraint){}
    Permutation(){
       this->sortConstraint = NULL;
    }
    Permutation(int tupleSplit): tupleSplit(tupleSplit) {}
    void insert(std::vector<T> tup){
        d.push_back(tup);
	if (sortConstraint != NULL){
	    std::sort(d.begin(),d.end(),sortConstraint);
	}
    }
    int get(std::vector<T> tup){
        typename std::vector<std::vector<T>>::iterator it;
    	if (tupleSplit == 0){
	    it = std::find(d.begin(),d.end(),tup);
	}else{
	    it = std::find_if(d.begin(),d.end(),[this,&tup](std::vector<T> &a){
			        for(int i=0; i < tupleSplit; i++){
				    if(a[i] != tup[i]) return false;
				}
				return true;
			    });
	}
	if (it == d.end()) {
	    std::stringstream ss;
	    assert(0 && "Tuple does not exist");
	}
	if (tupleSplit == 0) return it - d.begin();
	else return (*it)[tupleSplit];
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
// This abstraction is helpful for DIA and ELL
// When a tuple is inserted,say {-1}
// it gets ordered to 0
// then after {-2} is inserted,
// {-1} becomees 1
// {-2} becomes 0
// This is all dependent on the comparator.
// If the comparator is ascending, we have 
// the example shown above.
template <typename T>
class POS {
private:
    Comparator<T> constraint;
    std::vector<std::vector<T>> d;
public:
    POS() {}
    void insert (std::vector<T> tup){
        typename std::vector<std::vector<T>>::iterator it;
	it = std::find(d.begin(),d.end(),tup);
	if (it == d.end()){
	    d.push_back(tup);
	    if (constraint!= NULL)
	        std::sort(d.begin(),d.end(),constraint);
	    else 
	        std::sort(d.begin(),d.end());
	}
    }
    int get(std::vector<T> tup){
        typename std::vector<std::vector<T>>::iterator it;
	it = std::find(d.begin(),d.end(),tup);
	if (it == d.end()) {
	    assert(0 &&  "tuple does not exist");
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

#endif
