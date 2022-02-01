// Compile this file manually 
// g++ coo_csr_gen_code.cc -g -std=c++11 -o coo_csr_gen_code

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
	    ss << "}\n";
	}
	return ss.str();
    }
};

// This code turns off SSA
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
#undef s0
#undef s_0
#undef s1
#undef s_1
#undef s2
#undef s_2
#undef s3
#undef s_3
#undef s4
#undef s_4
#undef s5
#undef s_5
// Modification here.
#define s_0(n, tv1, tv2)   P->insert({row1(n), col1(n)}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6)   s_0(a1, a3, a5);
#define s_1(n)   rowptr(row1(n)) = min(rowptr(row1(n)),P(row1(n), col1(n))) 
#define s1(__x0, a1, tv2, __x2, __x3, __x4, __x5, __x6, _x7)   s_1(a1);
#define s_2(n)   rowptr(row1(n) + 1) = max(rowptr(row1(n) + 1),P(row1(n), col1(n)) + 1) 
#define s2(__x0, a1, tv2, __x2, __x3, __x4, __x5, __x6, _x7)   s_2(a1);
#define s_3(n, k)   col2(k)=col1(n) 

// There is a bug here... We had to fix mmanually by replacing a3 by __x4 in s_3 
#define s3(__x0, a1, tv2, __x2, a3, __x4, __x5, __x6, _x7)   s_3(a1, __x4);
#define s_4(e1, e2)   if ( not (rowptr(e1) <= rowptr(e2))){rowptr(e2) = rowptr(e1);} 
#define s4(__x0, a1, __x2, a3, __x4, __x5, __x6)   s_4(a1, a3);
// Manual Modification here
#define s_5(n, k)   ACSR[k] = ACOO[n] 
// There is a bug here... We had to fix mmanually by replacing a3 by __x4 in s_5.
// Figure out how to avoid generating this buggy code
// It might have something to do with what happens after nested ufs are replaced by 
// tuple variable assignments.
#define s5(__x0, a1, tv2, __x2, a3, __x4, __x5, __x6, _x7)   s_5(a1, __x4);

// Deleted undef macros, causing errors
// Data access change here
#define P(t0,t1) P->get({t0,t1})
// Possible bug here  something is going on with tv2 been 
// swapped for tv3i, temporary fix has been made to fix bug
// This bug has to do with nested ufs in a multi variate
// uninterpreted function. TODO: Test cases in iegenlib 
// must cover this corner case and pass correctly!
#define P_2(__tv0, __tv1, __tv2, __tv3) P(__tv2 + 1, __tv3)
#define P_3(__tv0, __tv1, __tv2, __tv3) P(__tv2, __tv3)
// Same swap happening here as well
#define P_4(__tv0, __tv1, __tv2, __tv3) P(__tv2, __tv3)
// This is flipped too
#define P_5(__tv0, __tv1, __tv2, __tv3) P(__tv2 + 1, __tv3)
#define col1(t0) col1[t0]
#define col1_1(__tv0, __tv1) col1(__tv1)
#define col2(t0) col2[t0]
#define col2_6(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) col2(__tv5)
#define row1(t0) row1[t0]
#define row1_0(__tv0, __tv1) row1(__tv1)
#define rowptr(t0) rowptr[t0]
#define rowptr_7(__tv0, __tv1, __tv2) rowptr(__tv2)
#define rowptr_8(__tv0, __tv1, __tv2) rowptr(__tv2 + 1)

bool compare_array(int *ar1, int *ar2, int length) {
  for (int i = 0; i < length; i++) {
    if (ar1[i] != ar2[i]) {
      return false;
    }
  }
  return true;
}

bool compare_array_d(double *ar1, double *ar2, int length) {
  for (int i = 0; i < length; i++) {
    if (ar1[i] != ar2[i]) {
      return false;
    }
  }
  return true;
}
int main (){
// Variable type definition begin
int t1 = 5; 
int t2 = 0; 
int t3 = 0; 
int t4 = 0; 
int t5 = 0; 
int t6 = 0; 
int t7 = 0; 
int t8 = 0; 
int t9 = 0; 

Permutation<int> * P = new Permutation <int>([](std::vector<int>& a,std::vector<int>& b){
		    for(int i = 0; i < a.size(); i++){
		       if (a[i]  < b[i] ) return true;
		       else if (a[i] > b[i]) return false;
		    }
		    return false;
		    });
// Example extracted dissertation proposal
double ACOO[] = {4,1,7,6,1,3,5};
int row1[] = {0,0,1,1,3,3,3};
int col1[] = {1,0,1,0,1,2,3};

int NR = 4;
int NC = 4;
int NNZ = 7;

int *col2 = (int*) calloc(NNZ,sizeof(int));
int *rowptr = (int*) calloc(NR+1,sizeof(int));
double *ACSR = (double*) calloc(NNZ,sizeof(double));

if (NR >= 1 && NC >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s0(0,t2,0,t4,0,t6,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      s1(1,t2,t3,t4,0,0,0,0,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      s2(2,t2,t3,t4,0,0,0,0,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      t6=P_4(t1,t2,t3,t4);
      s3(3,t2,t3,t4,0,t6,0,0,0);
    }
  }
}
for(t2 = 0; t2 <= NR-1; t2++) {
  for(t4 = t2+1; t4 <= NR; t4++) {
    s4(4,t2,0,t4,0,0,0);
  }
}
if (NC >= 1 && NR >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (col1_1(t1,t2) >= 0 && row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1) {
      t3=row1_0(t1,t2);
      if (rowptr_8(t1,t2,t3) >= rowptr_7(t1,t2,t3)+1) {
        t4=col1_1(t1,t2);
        if (rowptr_8(t1,t2,t3) >= P_4(t1,t2,t3,t4)+1 && P_4(t1,t2,t3,t4) >= rowptr_7(t1,t2,t3)) {
          t6=P_4(t1,t2,t3,t4);
          if (col2_6(t1,t2,t3,t4,t5,t6) == col1_1(t1,t2)) {
            s5(5,t2,t3,t4,0,t6,0,0,0);
          }
        }
      }
    }
  }
}


int col2_res[] = {0,1,0,1,1,2,3};
int rowptr_res[] = { 0,2,4,4,7};
double ACSR_res[] = {1,4,6,7,1,3,5};
// Expected result
assert(compare_array(rowptr_res,rowptr,NR+1) && "Invalid rowptr conversion");
assert(compare_array_d(ACSR_res,ACSR,NNZ) && "Invalid Value array");
assert(compare_array(col2_res,col2,NNZ) && "Invalid Col2 conversion");

std::cout << "Successfully performed Synthesis\n";


return 0;
}
// Removed defs
