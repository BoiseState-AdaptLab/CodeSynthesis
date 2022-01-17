#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
Permutation<int> * P0 = new Permutation<int>();
int t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
#define ACSR(i,k) ACSR[k]
#define ACOO(n) ACOO[n]
#undef s0
#undef s_0
#undef s1
#undef s_1
#undef s2
#undef s_2
#undef s3
#undef s_3
#define s_0(i, k, tv2, tv3)   P0->insert({i, col2(k)}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8)   s_0(a1, a3, a5, a7);
#define s_1(i, k, n)   col1(n)=col2(k) 
#define s1(__x0, a1, __x2, a3, tv4, __x5, a5, __x7, __x8, __x9)   s_1(a1, a3, a5);
#define s_2(i, k, n)   row1(n)=i 
#define s2(__x0, a1, __x2, a3, tv4, __x5, a5, __x7, __x8, __x9)   s_2(a1, a3, a5);
#define s_3(i, k, n)   ACOO(n) = ACSR(i,k )  
#define s3(__x0, a1, __x2, a3, tv4, __x5, a5, __x7, __x8, __x9)   s_3(a1, a3, a5);
#undef col1_5
#undef col2_0
#undef row1_4
#undef P0_3
#undef rowptr_1
#undef rowptr_2
#define col1(t0) col1[t0]
#define col1_5(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6) col1(__tv6)
#define col2(t0) col2[t0]
#define col2_0(__tv0, __tv1, __tv2, __tv3) col2(__tv3)
#define row1(t0) row1[t0]
#define row1_4(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6) row1(__tv6)
#define P0(t0,t1) P0->get({t0,t1})
#define P0_3(__tv0, __tv1, __tv2, __tv3, __tv4) P0(__tv1, __tv4)
#define rowptr(t0) rowptr[t0]
#define rowptr_1(__tv0, __tv1) rowptr(__tv1)
#define rowptr_2(__tv0, __tv1) rowptr(__tv1 + 1)

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

int main(){
t1 = 3; 
t2 = 0; 
t3 = 0; 
t4 = 0; 
t5 = 0; 
t6 = 0; 
t7 = 0; 
t8 = 0; 
t9 = 0; 
t10 = 0; 


int NR = 4;
int NC = 4;
int NNZ = 7;
int col2[] = {0,1,0,1,1,2,3};
int rowptr[] = { 0,2,4,4,7};
double ACSR[] = {1,4,6,7,1,3,5};

int *col1 = (int*) calloc(NNZ,sizeof(int));
int *row1 = (int*) calloc(NNZ,sizeof(int));
double *ACOO = (double*) calloc(NNZ,sizeof(double));


if (NNZ >= 1 && NC >= 1) {
  for(t2 = 0; t2 <= NR-1; t2++) {
    for(t4 = rowptr_1(t1,t2); t4 <= rowptr_2(t1,t2)-1; t4++) {
      if (col2_0(t1,t2,t3,t4) >= 0 && NC >= col2_0(t1,t2,t3,t4)+1) {
        t8=col2_0(t1,t2,t3,t4);
        s0(0,t2,0,t4,0,t2,0,t8,0);
      }
    }
  }
  for(t2 = 0; t2 <= NR-1; t2++) {
    for(t4 = rowptr_1(t1,t2); t4 <= rowptr_2(t1,t2)-1; t4++) {
      if (NC >= col2_0(t1,t2,t3,t4)+1 && NNZ >= 1 && col2_0(t1,t2,t3,t4) >= 0) {
        t5=col2_0(t1,t2,t3,t4);
        if (P0_3(t1,t2,t3,t4,t5) >= 0 && NNZ >= P0_3(t1,t2,t3,t4,t5)+1) {
          t7=P0_3(t1,t2,t3,t4,t5);
          s1(1,t2,0,t4,t5,0,t7,0,0,0);
        }
      }
    }
  }
  for(t2 = 0; t2 <= NR-1; t2++) {
    for(t4 = rowptr_1(t1,t2); t4 <= rowptr_2(t1,t2)-1; t4++) {
      if (NC >= col2_0(t1,t2,t3,t4)+1 && NNZ >= 1 && col2_0(t1,t2,t3,t4) >= 0) {
        t5=col2_0(t1,t2,t3,t4);
        if (P0_3(t1,t2,t3,t4,t5) >= 0 && NNZ >= P0_3(t1,t2,t3,t4,t5)+1) {
          t7=P0_3(t1,t2,t3,t4,t5);
          s2(2,t2,0,t4,t5,0,t7,0,0,0);
        }
      }
    }
  }
  for(t2 = 0; t2 <= NR-1; t2++) {
    for(t4 = rowptr_1(t1,t2); t4 <= rowptr_2(t1,t2)-1; t4++) {
      if (col2_0(t1,t2,t3,t4) >= 0 && NNZ >= 1 && NC >= col2_0(t1,t2,t3,t4)+1) {
        t5=col2_0(t1,t2,t3,t4);
        if (P0_3(t1,t2,t3,t4,t5) >= 0 && NNZ >= P0_3(t1,t2,t3,t4,t5)+1) {
          t7=P0_3(t1,t2,t3,t4,t5);
          if (row1_4(t1,t2,t3,t4,t5,t6,t7) == t2 && col1_5(t1,t2,t3,t4,t5,t6,t7) == t5) {
            s3(3,t2,0,t4,t5,0,t7,0,0,0);
          }
        }
      }
    }
  }
}


double ACOO_res[] =  {1,4,6,7,1,3,5};
int row1_res[] = {0,0,1,1,3,3,3};
int col1_res[] = {0,1,0,1,1,2,3};
assert(compare_array(row1_res,row1,NNZ) && "Invalid rowptr conversion");
assert(compare_array_d(ACOO_res,ACOO,NNZ) && "Invalid Value array");
assert(compare_array(col1_res,col1,NNZ) && "Invalid Col2 conversion");


std::cout << "Successfully performed Synthesis\n";
return 0;
}
#undef s0
#undef s_0
#undef s1
#undef s_1
#undef s2
#undef s_2
#undef s3
#undef s_3
#undef col1_5
#undef col2_0
#undef row1_4
#undef P0_3
#undef rowptr_1
#undef rowptr_2
