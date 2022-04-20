#include <opt_synth.h> 
// Optimization summary.
// Removing Permutes that don't affect reordering. 
// 1. This will be called optimize out Case 1 involving Permutes
// in this case P0 will be removed. We need a mathematical proof
// that permute functions dont need to be checked since they dont
// exist in reality except for reordering
// 2. Redundant Statement Elimination
// 3. Read reduction fusion
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
//This was modified to include another macro for indirection see bench_harness for how this is intended to be used.
#define ACOO(n) EX_ACOO(n)
//This was modified to include another macro for indirection see bench_harness for how this is intended to be used.
#define ACSR(i,k) EX_ACSR(k)
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
#define s_1(n, tv1, tv2)   P1->insert({row1(n), col1(n)}) 
#define s1(__x0, a1, __x2, a3, __x4, a5, __x6)   s_1(a1, a3, a5);
#define s_2(n,k)   col2(k)=col1(n) 
#define s2(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_2(a1,__x4);
#define s_4(n,k)   rowptr(row1(n)) = min(rowptr(row1(n)),k) 
#define s4(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_4(a1,__x4);
#define s_5(n,k)   rowptr(row1(n) + 1) = max(rowptr(row1(n) + 1),k + 1) 
#define s5(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_5(a1,__x4);
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
//This was modified to include another macro for indirection see bench_harness for how this is intended to be used.
#define col1(t0) EX_COL1(t0)
#define col1_1(__tv0, __tv1) col1(__tv1)
//This was modified to include another macro for indirection see bench_harness for how this is intended to be used.
#define col2(t0) EX_COL2(t0)
#define col2_4(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6, __tv7) col2(__tv7)
//This was modified to include another macro for indirection see bench_harness for how this is intended to be used.
#define row1(t0) EX_ROW1(t0)
#define row1_0(__tv0, __tv1) row1(__tv1)
//This was modified to include another macro for indirection see bench_harness for how this is intended to be used.
#define rowptr(t0) EX_ROWPTR(t0)
#define rowptr_5(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) rowptr(__tv5)
#define rowptr_6(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) rowptr(__tv5 + 1)
#define P0(t0,t1) P0->get({t0,t1})
#define P0_2(__tv0, __tv1, __tv2, __tv3) P0(__tv2, __tv3)

int t1 = 4; 
int t2 = 0; 
int t3 = 0; 
int t4 = 0; 
int t5 = 0; 
int t6 = 0; 
int t7 = 0; 
int t8 = 0; 
int t9 = 0; 



Permutation<int,P1_Comp>* P1 = new Permutation <int,P1_Comp>();
//  std::cerr << "S0: Permute Code Ended\n";
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (col1_1(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s1(1,t2,0,t4,0,t6,0);
    }
  }
// Fusing copy, rowptr and col2 synthesis
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      t6 = t3;
      t8=P1_3(t1,t2,t3,t4);
      s2(2,t2,t3,t4,t8,0,0,0,0);
      s4(5,t2,t3,t4,t8,0,0,0,0);
      s5(5,t2,t3,t4,t8,0,0,0,0);
      s7(7,t2,t3,t4,0,t6,0,t8,0);
    }
  }
//  std::cerr << "S5: completed \n";
for(t2 = 0; t2 <= NR-1; t2++) {
  for(t4 = t2+1; t4 <= NR; t4++) {
    s6(6,t2,0,t4,0,0,0);
  }
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
#undef rowptr_6
#undef P0_2

