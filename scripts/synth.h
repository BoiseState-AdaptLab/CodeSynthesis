#ifndef SYNTH_HEADER
#define SYNTH_HEADER
#include <functional>
#include <algorithm>
#include <set>
#include <vector>
#include <assert.h>
#include <iostream>
#include <string>
#include <sstream>
#define MAX_LEVEL 5
template <typename T, typename C, typename E>class SkipList;
template <typename T> class SkipNode;
template <typename T,typename C = std::less<std::vector<T> > >   
class Permutation{   
private:    
    std::set<std::vector<T>,C> d;   
    int tupleSplit = 0;   
public:   
    Permutation(C comp):d(comp), tupleSplit(0){}   
    Permutation(int tupleSplit): tupleSplit(tupleSplit) {}   
    void insert(std::vector<T> tup){   
        d.insert(tup);   
    }   
    int get(std::vector<T> tup){   
        typename std::set<std::vector<T>>::iterator it;   
    	if (tupleSplit == 0){   
	    it = d.find(tup);   
	}else{   
	    it = std::find_if(d.begin(),d.end(),
[this,&tup](const std::vector<T> &a){   
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

#endif   
