#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
#define ACOO(n) ACOO[n]
#define ACSR(i,k) ACSR[k]
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
#undef s6
#undef s_6
#undef s7
#undef s_7
//This was modified to tupleSplit
#define s_0(n, tv1, tv2)   P0->insert({row1(n), col1(n), row1(n)}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6)   s_0(a1, a3, a5);
#define s_1(n, tv1, tv2)   P1->insert({row1(n), col1(n)}) 
#define s1(__x0, a1, __x2, a3, __x4, a5, __x6)   s_1(a1, a3, a5);
#define s_2(n)   col2(P1(row1(n), col1(n)))=col1(n) 
#define s2(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_2(a1);
#define s_3(n, i, k)   col2(k)=col1(n) 
#define s3(__x0, a1, tv2, tv3, __x4, a3, __x6, a5, __x8)   s_3(a1, a3, a5);
#define s_4(n)   rowptr(row1(n)) = min(rowptr(row1(n)),P1(row1(n), col1(n))) 
#define s4(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_4(a1);
#define s_5(n)   rowptr(row1(n) + 1) = max(rowptr(row1(n) + 1),P1(row1(n), col1(n)) + 1) 
#define s5(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_5(a1);
#define s_6(e1, e2)   if ( not (rowptr(e1) <= rowptr(e2))){rowptr(e2) = rowptr(e1);} 
#define s6(__x0, a1, __x2, a3, __x4, __x5, __x6)   s_6(a1, a3);
#define s_7(n, i, k)   ACSR(i,k) = ACOO(n ) 
#define s7(__x0, a1, tv2, tv3, __x4, a3, __x6, a5, __x8)   s_7(a1, a3, a5);

#undef P1_3
#undef col1_1
#undef col2_4
#undef row1_0
#undef rowptr_5
#undef rowptr_6
#undef P0_2
#define P1(t0,t1) P1->get({t0,t1})
#define P1_3(__tv0, __tv1, __tv2, __tv3) P1(__tv2, __tv3)
#define col1(t0) col1[t0]
#define col1_1(__tv0, __tv1) col1(__tv1)
#define col2(t0) col2[t0]
#define col2_4(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6, __tv7) col2(__tv7)
#define row1(t0) row1[t0]
#define row1_0(__tv0, __tv1) row1(__tv1)
#define rowptr(t0) rowptr[t0]
#define rowptr_5(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) rowptr(__tv5)
#define rowptr_6(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) rowptr(__tv5 + 1)
#define P0(t0,t1) P0->get({t0,t1})
#define P0_2(__tv0, __tv1, __tv2, __tv3) P0(__tv2, __tv3)

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
int t1 = 4; 
int t2 = 0; 
int t3 = 0; 
int t4 = 0; 
int t5 = 0; 
int t6 = 0; 
int t7 = 0; 
int t8 = 0; 
int t9 = 0; 

double ACOO[] = {4,1,7,6,1,3,5};
int row1[] = {3,0,1,1,3,3,0};
int col1[] = {3,0,1,0,1,2,1};

int NR = 4;
int NC = 4;
int NNZ = 7;

int *col2 = (int*) calloc(NNZ,sizeof(int));
int *rowptr = (int*) calloc(NR+1,sizeof(int));
double *ACSR = (double*) calloc(NNZ,sizeof(double));

Permutation<int> * P0 = new Permutation <int>(2);

Permutation<int>* P1 = new Permutation <int>([]( std::vector<int>& a, std::vector<int>& b){
if (a[0] < b[0] )  return true;
return false;
});


if (NC >= 1 && NR >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (col1_1(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s0(0,t2,0,t4,0,t6,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (col1_1(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s1(1,t2,0,t4,0,t6,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (row1_0(t1,t2) == P0_2(t1,t2,t3,t4)) {
        s2(2,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (col1_1(t1,t2) >= 0 && NR >= row1_0(t1,t2)+1 && NC >= col1_1(t1,t2)+1 && row1_0(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (P0_2(t1,t2,t3,t4) == t3) {
        t6=row1_0(t1,t2);
        t8=P1_3(t1,t2,t3,t4);
        s3(3,t2,t3,t4,0,t6,0,t8,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (row1_0(t1,t2) == P0_2(t1,t2,t3,t4)) {
        s4(4,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (row1_0(t1,t2) == P0_2(t1,t2,t3,t4)) {
        s5(5,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
}
for(t2 = 0; t2 <= NR-1; t2++) {
  for(t4 = t2+1; t4 <= NR; t4++) {
    s6(6,t2,0,t4,0,0,0);
  }
}
if (NR >= 1 && NC >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && col1_1(t1,t2) >= 0 && NR >= row1_0(t1,t2)+1) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (P0_2(t1,t2,t3,t4) == t3) {
        t6=P0_2(t1,t2,t3,t4);
        if (P1_3(t1,t2,t3,t4) >= rowptr_5(t1,t2,t3,t4,t5,t6) && rowptr_6(t1,t2,t3,t4,t5,t6) >= P1_3(t1,t2,t3,t4)+1) {
          t8=P1_3(t1,t2,t3,t4);
          if (col2_4(t1,t2,t3,t4,t5,t6,t7,t8) == col1_1(t1,t2)) {
            s7(7,t2,t3,t4,0,t6,0,t8,0);
          }
        }
      }
    }
  }
}


int col2_res[] = {0,1,1,0,3,1,2};
int rowptr_res[] = { 0,2,4,4,7};
double ACSR_res[] = {1,5,7,6,4,1,3};
// Expected result
assert(compare_array(rowptr_res,rowptr,NR+1) && "Invalid rowptr conversion");
assert(compare_array_d(ACSR_res,ACSR,NNZ) && "Invalid Value array");
assert(compare_array(col2_res,col2,NNZ) && "Invalid Col2 conversion");

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
#undef s4
#undef s_4
#undef s5
#undef s_5
#undef s6
#undef s_6
#undef s7
#undef s_7
#undef P1_3
#undef col1_1
#undef col2_4
#undef row1_0
#undef rowptr_5
#undef rowptr_6
#undef P0_2

