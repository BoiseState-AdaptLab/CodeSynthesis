#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
Permutation<int> * P0 = new Permutation <int>(2);
Permutation<int>* P1 = new Permutation <int>([]( std::vector<int>& a, std::vector<int>& b){
if (a[1] < b[1] )    return true;
return false;
});
#define ACOO(n) ACOO[n]
#define ACSC(j,k) ACSC[k]
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
#define s_0(n, tv1, tv2)   P0->insert({row1(n), col1(n), col1(n)}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6)   s_0(a1, a3, a5);
#define s_1(n, tv1, tv2)   P1->insert({row1(n), col1(n)}) 
#define s1(__x0, a1, __x2, a3, __x4, a5, __x6)   s_1(a1, a3, a5);
#define s_2(n)   colptr(P0(row1(n), col1(n))) = min(colptr(P0(row1(n), col1(n))),P1(row1(n), col1(n))) 
#define s2(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_2(a1);
#define s_3(n)   colptr(P0(row1(n), col1(n)) + 1) = max(colptr(P0(row1(n), col1(n)) + 1),P1(row1(n), col1(n)) + 1) 
#define s3(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_3(a1);
#define s_4(n)   row4(P1(row1(n), col1(n)))=row1(n) 
#define s4(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_4(a1);
#define s_5(n, j, k)   row4(k)=row1(n) 
#define s5(__x0, a1, tv2, tv3, __x4, a3, __x6, a5, __x8)   s_5(a1, a3, a5);
#define s_6(e1, e2)   if ( not (colptr(e1) <= colptr(e2))){colptr(e2) = colptr(e1);} 
#define s6(__x0, a1, __x2, a3, __x4, __x5, __x6)   s_6(a1, a3);
#define s_7(n, j, k)   ACSC(j,k) = ACOO(n ) 
#define s7(__x0, a1, tv2, tv3, __x4, a3, __x6, a5, __x8)   s_7(a1, a3, a5);

#undef P0_2
#undef P1_3
#undef col1_1
#undef colptr_5
#undef colptr_6
#undef row1_0
#undef row4_4
#define P0(t0,t1) P0->get({t0,t1})
#define P0_2(__tv0, __tv1, __tv2, __tv3) P0(__tv2, __tv3)
#define P1(t0,t1) P1->get({t0,t1})
#define P1_3(__tv0, __tv1, __tv2, __tv3) P1(__tv2, __tv3)
#define col1(t0) col1[t0]
#define col1_1(__tv0, __tv1) col1(__tv1)
#define colptr(t0) colptr[t0]
#define colptr_5(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) colptr(__tv5)
#define colptr_6(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) colptr(__tv5 + 1)
#define row1(t0) row1[t0]
#define row1_0(__tv0, __tv1) row1(__tv1)
#define row4(t0) row4[t0]
#define row4_4(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6, __tv7) row4(__tv7)

t1 = 7; 
t2 = 0; 
t3 = 0; 
t4 = 0; 
t5 = 0; 
t6 = 0; 
t7 = 0; 
t8 = 0; 
t9 = 0; 

if (NC >= 1 && NR >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NR >= row1_0(t1,t2)+1) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s0(0,t2,0,t4,0,t6,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NR >= row1_0(t1,t2)+1) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s1(1,t2,0,t4,0,t6,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (col1_1(t1,t2) == P0_2(t1,t2,t3,t4)) {
        s2(2,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (col1_1(t1,t2) == P0_2(t1,t2,t3,t4)) {
        s3(3,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (col1_1(t1,t2) == P0_2(t1,t2,t3,t4)) {
        s4(4,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NR >= row1_0(t1,t2)+1 && NC >= col1_1(t1,t2)+1 && row1_0(t1,t2) >= 0 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (P0_2(t1,t2,t3,t4) == t4) {
        t6=col1_1(t1,t2);
        t8=P1_3(t1,t2,t3,t4);
        s5(5,t2,t3,t4,0,t6,0,t8,0);
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
    if (col1_1(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && row1_0(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (P0_2(t1,t2,t3,t4) == col1_1(t1,t2)) {
        t6=col1_1(t1,t2);
        if (P1_3(t1,t2,t3,t4) >= colptr_5(t1,t2,t3,t4,t5,t6) && colptr_6(t1,t2,t3,t4,t5,t6) >= P1_3(t1,t2,t3,t4)+1) {
          t8=P1_3(t1,t2,t3,t4);
          if (row4_4(t1,t2,t3,t4,t5,t6,t7,t8) == t3) {
            s7(7,t2,t3,t4,0,t6,0,t8,0);
          }
        }
      }
    }
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
#undef P0_2
#undef P1_3
#undef col1_1
#undef colptr_5
#undef colptr_6
#undef row1_0
#undef row4_4
