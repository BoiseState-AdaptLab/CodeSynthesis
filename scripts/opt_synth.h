#ifndef SYNTH_HEADER
#define SYNTH_HEADER
#include <functional>
#include <algorithm>
#include <vector>
#include <set>
#include <assert.h>
#include <iostream>
#include <string>
#include <sstream>
// Define P Data structure
template <typename T>
using Comparator = std::function<bool (std::vector<T>&,std::vector<T>&)>;
struct P1_Comp{
bool operator() (const std::vector<int>& a,const std::vector<int>& b){
if (a[0] < b[0] )  return true;
return a < b;
}
};


template <typename T,typename C = std::less<std::vector<T> > >
class Permutation{
private: 
    std::set<std::vector<T>,C> d;
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
        d.insert(tup);
    }
    int get(std::vector<T> tup){
        typename std::set<std::vector<T>>::iterator it;
    	if (tupleSplit == 0){
	    it = d.find(tup);
	}else{
	    it = std::find_if(d.begin(),d.end(),[this,&tup](const std::vector<T> &a){
			        for(int i=0; i < tupleSplit; i++){
				    if(a[i] != tup[i]) return false;
				}
				return true;
			    });
	}
	if (it == d.end()) {
	    assert(0 && "Tuple get tuple does not exist");
	}
	if (tupleSplit == 0) return std::distance(d.begin(),it);
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

#endif //SYNTH_HEADER
