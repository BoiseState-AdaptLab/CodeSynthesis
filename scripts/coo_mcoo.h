#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b

#include <user_defs.h>

Permutation<int> * P0 = new Permutation <int>([](std::vector<int>& a,std::vector<int>& b){
if(-MORTON(a[0], a[1]) + MORTON(b[0], b[1]) - 1>= 0)
return true;
if(-MORTON(a[0], a[1]) + MORTON(b[0], b[1]) - 1>= 0)
return true;
    

return false;
});
#define ACOO(n) ACOO[n]
#define AMCOO(n1) AMCOO[n]
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
#define s_0(n, tv1, tv2)   P0->insert({row1(n), col1(n)}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6)   s_0(a1, a3, a5);
#define s_1(n)   col3(P0(row1(n), col1(n)))=col1(n) 
#define s1(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_1(a1);
#define s_2(n, n1)   col3(n1)=col1(n) 
#define s2(__x0, a1, tv2, tv3, __x4, a3, __x6, __x7, __x8)   s_2(a1, a3);
#define s_3(n)   row3(P0(row1(n), col1(n)))=row1(n) 
#define s3(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_3(a1);
#define s_4(n, n1)   row3(n1)=row1(n) 
#define s4(__x0, a1, tv2, tv3, __x4, a3, __x6, __x7, __x8)   s_4(a1, a3);
#define s_5(n, n1)   AMCOO(n1) = ACOO(n ) 
#define s5(__x0, a1, tv2, tv3, __x4, a3, __x6, __x7, __x8)   s_5(a1, a3);

#undef P0_2
#undef col1_1
#undef col3_3
#undef row1_0
#undef row3_4
#define P0(t0,t1) P0->get({t0,t1})
#define P0_2(__tv0, __tv1, __tv2, __tv3) P0(__tv2, __tv3)
#define col1(t0) col1[t0]
#define col1_1(__tv0, __tv1) col1(__tv1)
#define col3(t0) col3[t0]
#define col3_3(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) col3(__tv5)
#define row1(t0) row1[t0]
#define row1_0(__tv0, __tv1) row1(__tv1)
#define row3(t0) row3[t0]
#define row3_4(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) row3(__tv5)

int t1 = 5; 
int t2 = 0; 
int t3 = 0; 
int t4 = 0; 
int t5 = 0; 
int t6 = 0; 
int t7 = 0; 
int t8 = 0; 
int t9 = 0; 

if (NR >= 1 && NC >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s0(0,t2,0,t4,0,t6,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0 && row1_0(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (NNZ >= P0_2(t1,t2,t3,t4)+1 && P0_2(t1,t2,t3,t4) >= 0) {
        s1(1,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (NNZ >= P0_2(t1,t2,t3,t4)+1 && P0_2(t1,t2,t3,t4) >= 0) {
        t6=P0_2(t1,t2,t3,t4);
        s2(2,t2,t3,t4,0,t6,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0 && row1_0(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (NNZ >= P0_2(t1,t2,t3,t4)+1 && P0_2(t1,t2,t3,t4) >= 0) {
        s3(3,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (NNZ >= P0_2(t1,t2,t3,t4)+1 && P0_2(t1,t2,t3,t4) >= 0) {
        t6=P0_2(t1,t2,t3,t4);
        s4(4,t2,t3,t4,0,t6,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0 && row1_0(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (NNZ >= P0_2(t1,t2,t3,t4)+1 && P0_2(t1,t2,t3,t4) >= 0) {
        t6=P0_2(t1,t2,t3,t4);
        if (row3_4(t1,t2,t3,t4,t5,t6) == t3 && col3_3(t1,t2,t3,t4,t5,t6) == col1_1(t1,t2)) {
          s5(5,t2,t3,t4,0,t6,0,0,0);
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
#undef P0_2
#undef col1_1
#undef col3_3
#undef row1_0
#undef row3_4
