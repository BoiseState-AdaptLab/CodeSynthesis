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
#define s_0(n, tv1, tv2)   P.insert(row1(n), col1(n)) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6)   s_0(a1, a3, a5);
#define s_1(n)   rowptr__w__0(row1(n)) = min(rowptr(row1(n)),P(row1(n), col1(n))) 
#define s1(__x0, a1, tv2, __x2, __x3, __x4, __x5, __x6, _x7)   s_1(a1);
#define s_2(n)   rowptr__w__1(row1(n) + 1) = max(rowptr__w__0(row1(n) + 1),P(row1(n), col1(n)) + 1) 
#define s2(__x0, a1, tv2, __x2, __x3, __x4, __x5, __x6, _x7)   s_2(a1);
#define s_3(n, k)   col2(k)=col1(n) 
#define s3(__x0, a1, tv2, __x2, a3, __x4, __x5, __x6, _x7)   s_3(a1, a3);
#define s_4(e1, e2)   if ( not (rowptr(e1) <= rowptr(e2))){rowptr(e2) = rowptr__w__1(e1);} 
#define s4(__x0, a1, __x2, a3, __x4, __x5, __x6)   s_4(a1, a3);
#define s_5(n, k)   ACSR(n,k) = ACOO(n,k) 
#define s5(__x0, a1, tv2, __x2, a3, __x4, __x5, __x6, _x7)   s_5(a1, a3);

#undef P(t0,t1)
#undef P_2(__tv0, __tv1, __tv2, __tv3)
#undef P_3(__tv0, __tv1, __tv2, __tv3)
#undef P_4(__tv0, __tv1, __tv2, __tv3)
#undef P_5(__tv0, __tv1, __tv2, __tv3)
#undef col1(t0)
#undef col1_1(__tv0, __tv1)
#undef col2(t0)
#undef col2_6(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5)
#undef row1(t0)
#undef row1_0(__tv0, __tv1)
#undef rowptr(t0)
#undef rowptr_7(__tv0, __tv1, __tv2)
#undef rowptr_8(__tv0, __tv1, __tv2)
#define P(t0,t1) P[t0][t1]
#define P_2(__tv0, __tv1, __tv2, __tv3) P(__tv3 + 1, __tv2)
#define P_3(__tv0, __tv1, __tv2, __tv3) P(__tv2, __tv3)
#define P_4(__tv0, __tv1, __tv2, __tv3) P(__tv3, __tv2)
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
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (P_3(t1,t2,t3,t4) >= P_2(t1,t2,t3,t4)+1) {
        s1(1,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (P_3(t1,t2,t3,t4) >= P_2(t1,t2,t3,t4)+1) {
        s2(2,t2,t3,t4,0,0,0,0,0);
      }
    }
  }
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_0(t1,t2) >= 0 && NC >= col1_1(t1,t2)+1 && NR >= row1_0(t1,t2)+1 && col1_1(t1,t2) >= 0) {
      t3=row1_0(t1,t2);
      t4=col1_1(t1,t2);
      if (P_3(t1,t2,t3,t4) >= P_5(t1,t2,t3,t4)+1) {
        t6=P_4(t1,t2,t3,t4);
        s3(3,t2,t3,t4,0,t6,0,0,0);
      }
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
#undef P(t0,t1)
#undef P_2(__tv0, __tv1, __tv2, __tv3)
#undef P_3(__tv0, __tv1, __tv2, __tv3)
#undef P_4(__tv0, __tv1, __tv2, __tv3)
#undef P_5(__tv0, __tv1, __tv2, __tv3)
#undef col1(t0)
#undef col1_1(__tv0, __tv1)
#undef col2(t0)
#undef col2_6(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5)
#undef row1(t0)
#undef row1_0(__tv0, __tv1)
#undef rowptr(t0)
#undef rowptr_7(__tv0, __tv1, __tv2)
#undef rowptr_8(__tv0, __tv1, __tv2)

