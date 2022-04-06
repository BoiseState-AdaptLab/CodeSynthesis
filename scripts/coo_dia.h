#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
POS <int> * pos = new POS<int>();
Permutation<int> * P0 = new Permutation<int>(2);
Permutation<int> * P1 = new Permutation<int>(2);
Permutation<int> * P2 = new Permutation<int>(2);
#define ACOO(n) ACOO[n]
#define ADIA(id,dd,kd) ADIA[kd]
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
#define s_0(n, tv1)   pos->insert({col1(n) - row1(n)}) 
#define s0(__x0, a1, __x2, a3, __x4, __x5, __x6, __x7, __x8)   s_0(a1, a3);
#define s_1(n)   P0->insert({row1(n), col1(n),row1(n)}) 
#define s1(__x0, a1, __x2, __x3, __x4, __x5, __x6, __x7, __x8)   s_1(a1);
#define s_2(n)   P1->insert({row1(n), col1(n),pos(col1(n) - row1(n))}) 
#define s2(__x0, a1, __x2, __x3, __x4, __x5, __x6, __x7, __x8)   s_2(a1);
#define s_3(n)   P2->insert({row1(n), col1(n),ND*P0(row1(n), col1(n)) + P1(row1(n), col1(n))}) 
#define s3(__x0, a1, __x2, __x3, __x4, __x5, __x6, __x7, __x8)   s_3(a1);
#define s_4(n)   off(P1(row1(n), col1(n)))=-P0(row1(n), col1(n)) + col1(n) 
#define s4(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8, __x9, __x10)   s_4(a1);
#define s_5(n)   off(P1(row1(n), col1(n))) = max(off(P1(row1(n), col1(n))),-P0(row1(n), col1(n))) 
#define s5(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8, __x9, __x10)   s_5(a1);
#define s_6(n)   off(P1(row1(n), col1(n))) = min(off(P1(row1(n), col1(n))),NC - P0(row1(n), col1(n)) - 1) 
#define s6(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8, __x9, __x10)   s_6(a1);
#define s_7(n, id, dd, kd)   ADIA(id,dd,kd) = ACOO(n ) 
#define s7(__x0, a1, tv2, tv3, __x4, a3, __x6, a5, __x8, a7, __x10)   	s_7(a1, a3, a5, a7);

#undef P0_2
#undef P1_3
#undef col1_0
#undef off_6
#undef pos_4
#undef row1_1
#undef P2_5
#define P0(t0,t1) P0->get({t0,t1})
#define P0_2(__tv0, __tv1, __tv2, __tv3) P0(__tv2, __tv3)
#define P1(t0,t1) P1->get({t0,t1})
#define P1_3(__tv0, __tv1, __tv2, __tv3) P1(__tv2, __tv3)
#define col1(t0) col1[t0]
#define col1_0(__tv0, __tv1) col1(__tv1)
#define off(t0) off[t0]
#define off_6(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6, __tv7) off(__tv7)
#define pos(t0) pos->get({t0})
#define pos_4(__tv0, __tv1, __tv2, __tv3) pos(-__tv2 + __tv3)
#define row1(t0) row1[t0]
#define row1_1(__tv0, __tv1) row1(__tv1)
#define P2(t0,t1) P2->get({t0,t1})
#define P2_5(__tv0, __tv1, __tv2, __tv3) P2(__tv2, __tv3)

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

if (ND >= 1 && NC >= 1 && NR >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_0(t1,t2)+1 && col1_0(t1,t2) >= 0) {
      t4=-row1_1(t1,t2)+col1_0(t1,t2);
      s0(0,t2,0,t4,0,0,0,0,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_0(t1,t2)+1 && col1_0(t1,t2) >= 0) {
      s1(1,t2,0,0,0,0,0,0,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_0(t1,t2)+1 && col1_0(t1,t2) >= 0) {
      s2(2,t2,0,0,0,0,0,0,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col1_0(t1,t2)+1 && col1_0(t1,t2) >= 0) {
      s3(3,t2,0,0,0,0,0,0,0);
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NR >= row1_1(t1,t2)+1 && NC >= col1_0(t1,t2)+1 && col1_0(t1,t2) >= 0 && row1_1(t1,t2) >= 0) {
      t3=row1_1(t1,t2);
      t4=col1_0(t1,t2);
      if (P1_3(t1,t2,t3,t4)+ND*row1_1(t1,t2) == P2_5(t1,t2,t3,t4) && P1_3(t1,t2,t3,t4) == pos_4(t1,t2,t3,t4) && row1_1(t1,t2) == P0_2(t1,t2,t3,t4) && P1_3(t1,t2,t3,t4) >= 0 && ND >= P1_3(t1,t2,t3,t4)+1) {
        s4(4,t2,t3,t4,0,0,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NR >= row1_1(t1,t2)+1 && NC >= col1_0(t1,t2)+1 && col1_0(t1,t2) >= 0 && row1_1(t1,t2) >= 0) {
      t3=row1_1(t1,t2);
      t4=col1_0(t1,t2);
      if (P1_3(t1,t2,t3,t4)+ND*row1_1(t1,t2) == P2_5(t1,t2,t3,t4) && P1_3(t1,t2,t3,t4) == pos_4(t1,t2,t3,t4) && row1_1(t1,t2) == P0_2(t1,t2,t3,t4) && P1_3(t1,t2,t3,t4) >= 0 && ND >= P1_3(t1,t2,t3,t4)+1) {
        s5(5,t2,t3,t4,0,0,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NR >= row1_1(t1,t2)+1 && NC >= col1_0(t1,t2)+1 && col1_0(t1,t2) >= 0 && row1_1(t1,t2) >= 0) {
      t3=row1_1(t1,t2);
      t4=col1_0(t1,t2);
      if (P1_3(t1,t2,t3,t4)+ND*row1_1(t1,t2) == P2_5(t1,t2,t3,t4) && P1_3(t1,t2,t3,t4) == pos_4(t1,t2,t3,t4) && row1_1(t1,t2) == P0_2(t1,t2,t3,t4) && P1_3(t1,t2,t3,t4) >= 0 && ND >= P1_3(t1,t2,t3,t4)+1) {
        s6(6,t2,t3,t4,0,0,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_1(t1,t2) >= 0 && NC >= col1_0(t1,t2)+1 && col1_0(t1,t2) >= 0 && NR >= row1_1(t1,t2)+1) {
      t3=row1_1(t1,t2);
      t4=col1_0(t1,t2);
      if (P1_3(t1,t2,t3,t4)+ND*row1_1(t1,t2) == P2_5(t1,t2,t3,t4) && P1_3(t1,t2,t3,t4) == pos_4(t1,t2,t3,t4) && row1_1(t1,t2) == P0_2(t1,t2,t3,t4) && ND >= P1_3(t1,t2,t3,t4)+1 && P1_3(t1,t2,t3,t4) >= 0) {
        t6=row1_1(t1,t2);
        t8=P2_5(t1,t2,t3,t4)-ND*t3;
        if (col1_0(t1,t2) == off_6(t1,t2,t3,t4,t5,t6,t7,t8)+P0_2(t1,t2,t3,t4)) {
          t10=-ND*off_6(t1,t2,t3,t4,t5,t6,t7,t8)+ND*col1_0(t1,t2)+t8;
	  s7(7,t2,t3,t4,0,t6,0,t8,0,t10,0);
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
#undef col1_0
#undef off_6
#undef pos_4
#undef row1_1
#undef P2_5
