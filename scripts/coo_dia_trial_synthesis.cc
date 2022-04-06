#include <m_synth.h> 
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
      std::cerr << "[" << i << "] => "<< ar1[i] << " , "<< ar2[i] << "\n";
      return false;
    }
  }
  return true;
}
int main () {
int t1 = 7; 
int t2 = 0; 
int t3 = 0; 
int t4 = 0; 
int t5 = 0; 
int t6 = 0; 
int t7 = 0; 
int t8 = 0; 
int t9 = 0; 
int t10 = 0; 
int t11 = 0; 

// Example extracted dissertation proposal
double ACOO[] = {4,1,7,6,1,3,5};
int row1[] = {0,0,1,1,2,2,3};
int col1[] = {1,0,1,0,1,3,3};

int NR = 4;
int NC = 4;
int NNZ = 7;
// This should be generated
int ND = 3;


int *off = (int*) calloc(ND,sizeof(int));
double *ADIA = (double*) calloc(12,sizeof(double));

#include <coo_dia.h>


int off_res[] = {-1,0,1};
double ADIA_res[] = {0,1,4,6,7,0,1,0,3,0,5,0};
// Expected result
assert(compare_array(off_res,off,ND) && "Invalid offset conversion");
assert(compare_array_d(ADIA_res,ADIA,12) && "Invalid Value array");

std::cout << "Successfully performed Synthesis\n";

}
