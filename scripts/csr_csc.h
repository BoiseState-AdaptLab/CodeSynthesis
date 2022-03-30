#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
Permutation<int> * P0 = new Permutation <int>([](std::vector<int>& a,std::vector<int>& b){
    return false;
});
Permutation<int> * P1 = new Permutation <int>([](std::vector<int>& a,std::vector<int>& b){
    return false;
});
Permutation<int> * P2 = new Permutation <int>([](std::vector<int>& a,std::vector<int>& b){
    return false;
});
#define ACOO(n) ACOO[n]
#define ADIA(q,r,k) ADIA[k]
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
#define s_0(n, tv1, tv2)   P0->insert({row1(n), col1(n)}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6, __x7, __x8)   s_0(a1, a3, a5);
#define s_1(n, tv1, tv2)   P1->insert({row1(n), col1(n)}) 
#define s1(__x0, a1, __x2, a3, __x4, a5, __x6, __x7, __x8)   s_1(a1, a3, a5);
#define s_2(n, tv1, tv2)   P2->insert({row1(n), col1(n)}) 
#define s2(__x0, a1, __x2, a3, __x4, a5, __x6, __x7, __x8)   s_2(a1, a3, a5);
#define s_3(n)   offset(P0(row1(n), col1(n)))=col1(n) - row1(n) 
#define s3(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8, __x9, __x10)   s_3(a1);
#define s_4(n, q)   offset(q)=col1(n) - row1(n) 
#define s4(__x0, a1, tv2, tv3, __x4, a3, __x6, __x7, __x8, __x9, __x10)   s_4(a1, a3);
#define s_5(n)   offset(P0(row1(n), col1(n))) = max(offset(P0(row1(n), col1(n))),-row1(n)) 
#define s5(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8, __x9, __x10)   s_5(a1);
#define s_6(n)   offset(P0(row1(n), col1(n))) = min(offset(P0(row1(n), col1(n))),NC - row1(n) - 1) 
#define s6(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8, __x9, __x10)   s_6(a1);
#define s_7(n, q, r, k)   ADIA(q,r,k) = ACOO(n ) 
#define s7(__x0, a1, tv2, tv3, __x4, a3, __x6, a5, __x8, a7, __x10)   s_7(a1, a3, a5, a7);

#undef P0_3
#undef col1_1
#undef offset_5
#undef row1_0
#undef P1_2
#undef P2_4
#define P0(t0,t1) P0->get({t0,t1})
#define P0_3(__tv0, __tv1, __tv2, __tv3) P0(__tv2, __tv3)
#define col1(t0) col1[t0]
#define col1_1(__tv0, __tv1) col1(__tv1)
#define offset(t0) offset[t0]
#define offset_5(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) offset(__tv5)
#define row1(t0) row1[t0]
#define row1_0(__tv0, __tv1) row1(__tv1)
#define P1(t0,t1) P1->get({t0,t1})
#define P1_2(__tv0, __tv1, __tv2, __tv3) P1(__tv2, __tv3)
#define P2(t0,t1) P2->get({t0,t1})
#define P2_4(__tv0, __tv1, __tv2, __tv3) P2(__tv2, __tv3)

t1 = 7; 
t2 = 0; 
t3 = 0; 
t4 = 0; 
t5 = 0; 
t6 = 0; 
t7 = 0; 
t8 = 0; 
t9 = 0; 
t10 = 0; 
t11 = 0; 

if (NR >= 1 && C >= 1 && NC >= 1 && R >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_1(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s0(0,t2,0,t4,0,t6,0,0,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_1(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s1(1,t2,0,t4,0,t6,0,0,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_1(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s2(2,t2,0,t4,0,t6,0,0,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NR >= row1_0(t1,t2)+1 && NC >= col1_1(t1,t2)+1 && row1_0(t1,t2) >= 0 && R >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (99*P0_3(t1,t2,t3,t4)+P1_2(t1,t2,t3,t4) == P2_4(t1,t2,t3,t4) && P1_2(t1,t2,t3,t4) == row1_0(t1,t2) && C >= P0_3(t1,t2,t3,t4)+1 && P0_3(t1,t2,t3,t4) >= 0) {
        s3(3,t2,t3,t4,0,0,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NR >= row1_0(t1,t2)+1 && NC >= col1_1(t1,t2)+1 && row1_0(t1,t2) >= 0 && R >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (99*P0_3(t1,t2,t3,t4)+P1_2(t1,t2,t3,t4) == P2_4(t1,t2,t3,t4) && P1_2(t1,t2,t3,t4) == row1_0(t1,t2) && C >= P0_3(t1,t2,t3,t4)+1 && P0_3(t1,t2,t3,t4) >= 0) {
        t6=P0_3(t1,t2,t3,t4);
        s4(4,t2,t3,t4,0,t6,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NR >= row1_0(t1,t2)+1 && NC >= col1_1(t1,t2)+1 && row1_0(t1,t2) >= 0 && R >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (99*P0_3(t1,t2,t3,t4)+P1_2(t1,t2,t3,t4) == P2_4(t1,t2,t3,t4) && P1_2(t1,t2,t3,t4) == row1_0(t1,t2) && C >= P0_3(t1,t2,t3,t4)+1 && P0_3(t1,t2,t3,t4) >= 0) {
        s5(5,t2,t3,t4,0,0,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_1(t1,t2)+1 && col1_1(t1,t2) >= 0 && row1_0(t1,t2) >= 0 && NR >= row1_0(t1,t2)+1 && R >= row1_0(t1,t2)+1) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (99*P0_3(t1,t2,t3,t4)+P1_2(t1,t2,t3,t4) == P2_4(t1,t2,t3,t4) && P1_2(t1,t2,t3,t4) == row1_0(t1,t2) && C >= P0_3(t1,t2,t3,t4)+1 && P0_3(t1,t2,t3,t4) >= 0) {
        s6(6,t2,t3,t4,0,0,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_1(t1,t2)+1 && row1_0(t1,t2) >= 0 && col1_1(t1,t2) >= 0 && NR >= row1_0(t1,t2)+1 && R >= row1_0(t1,t2)+1) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (99*P0_3(t1,t2,t3,t4)+P1_2(t1,t2,t3,t4) == P2_4(t1,t2,t3,t4) && P1_2(t1,t2,t3,t4) == row1_0(t1,t2) && P0_3(t1,t2,t3,t4) >= 0 && C >= P0_3(t1,t2,t3,t4)+1) {
        t6=P0_3(t1,t2,t3,t4);
        if (col1_1(t1,t2) == offset_5(t1,t2,t3,t4,t5,t6)+row1_0(t1,t2)) {
          t10=-offset_5(t1,t2,t3,t4,t5,t6)+t4+99*t6;
          s7(7,t2,t3,t4,0,t6,0,t3,0,t10,0);
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
#undef P0_3
#undef col1_1
#undef offset_5
#undef row1_0
#undef P1_2
#undef P2_4
