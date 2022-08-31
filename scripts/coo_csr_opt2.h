#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
ReorderStream* P1= new ReorderStream(2);
ComparatorInt P1Comp = [&](const int a,const int b){
if (a[0] < b[0] )    return true;
return false;
}; 
P1->setComparator(P1Comp);
#define ACOO(n,ii,jj) ACOO[n]
#define ACSR(ii,k,jj) ACSR[k]
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
#undef s8
#undef s_8
#undef s9
#undef s_9
#undef s10
#undef s_10
#define s_0(n, ii, jj)   P1->insert({ii, jj}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6, __x7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_0(a1, a3, a5);
#define s_1(0)   P1->sort() 
#define s1(__x0, __x1, __x2, __x3, __x4, __x5, __x6, __x7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_1(a1);
#define s_2(_no, _n, n, ii, jj, ii1, k)   rowptr(ii1) = min(rowptr(ii1),k) 
#define s2(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16)   s_2(a1, a3, a5, a7, a9, a11, a13);
#define s_3(_no, _n, n, ii, jj, ii1, k)   rowptr(ii1 + 1) = max(rowptr(ii1 + 1),k + 1) 
#define s3(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16)   s_3(a1, a3, a5, a7, a9, a11, a13);
#define s_4(_no, _n, n, ii, jj, ii1, k)   col2(k)=jj 
#define s4(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16)   s_4(a1, a3, a5, a7, a9, a11, a13);
#define s_5(_no, _n, n, ii, jj, ii1, k, jj1)   col2(k)=jj1 
#define s5(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, a15, __x16)   s_5(a1, a3, a5, a7, a9, a11, a13, a15);
#define s_6(_no, _n, n, ii, jj, ii1, k)   col2(k)=col1(n) 
#define s6(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16)   s_6(a1, a3, a5, a7, a9, a11, a13);
#define s_7(_no, _n, n, ii, jj, ii1, k)   col2(k) = max(col2(k),) 
#define s7(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16)   s_7(a1, a3, a5, a7, a9, a11, a13);
#define s_8(_no, _n, n, ii, jj, ii1, k)   col2(k) = min(col2(k),NC - 1) 
#define s8(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16)   s_8(a1, a3, a5, a7, a9, a11, a13);
#define s_9(e1, e2)   if ( not (rowptr(e1) <= rowptr(e2))){rowptr(e2) = rowptr(e1);} 
#define s9(__x0, a1, __x2, a3, __x4, __x5, __x6, __x7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_9(a1, a3);
#define s_10(_no, _n, n, ii, jj, ii1, k, jj1)   ACSR(ii,k,jj) = ACOO(n,ii,jj ) 
#define s10(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, a15, __x16)   s_10(a1, a3, a5, a7, a9, a11, a13, a15);

#undef P1DIM0_3
#undef P1DIM1_4
#undef P1MAP_2
#undef col1_1
#undef row1_0
#define P1DIM0(idx) P1->getDim(0,idx)
#define P1DIM0_3(__tv0, __tv1) P1DIM0(__tv1)
#define P1DIM1(idx) P1->getDim(1,idx)
#define P1DIM1_4(__tv0, __tv1) P1DIM1(__tv1)
#define P1MAP(idx) P1->getMap(idx)
#define P1MAP_2(__tv0, __tv1) P1MAP(__tv1)
#define col1(t0) col1[t0]
#define col1_1(__tv0, __tv1) col1(__tv1)
#define row1(t0) row1[t0]
#define row1_0(__tv0, __tv1) row1(__tv1)

t1 = 18; 
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
t12 = 0; 
t13 = 0; 
t14 = 0; 
t15 = 0; 
t16 = 0; 
t17 = 0; 

if (NR >= 1 && NC >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (col1_1(t1,t2) >= 0 && NR >= row1_0(t1,t2)+1 && NC >= col1_1(t1,t2)+1 && row1_0(t1,t2) >= 0) {
      t4=row1_0(t1,t2);
      t6=col1_1(t1,t2);
      s0(0,t2,0,t4,0,t6,0,0,0,0,0,0,0,0,0,0,0);
    }
  }
}
s1(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
for(t2 = 0; t2 <= P1SIZE-1; t2++) {
  t4=P1MAP_2(t1,t2);
  t8=P1DIM0_3(t1,t2);
  t10=P1DIM1_4(t1,t2);
  t14=P1MAP_2(t1,t2);
  s2(2,t2,0,t4,0,t2,0,t8,0,t10,0,t8,0,t14,0,0,0);
}
for(t2 = 0; t2 <= P1SIZE-1; t2++) {
  t4=P1MAP_2(t1,t2);
  t8=P1DIM0_3(t1,t2);
  t10=P1DIM1_4(t1,t2);
  t12=P1DIM0_3(t1,t2);
  t14=P1MAP_2(t1,t2);
  s3(5,t2,0,t4,0,t2,0,t8,0,t10,0,t12,0,t14,0,0,0);
}
for(t2 = 0; t2 <= P1SIZE-1; t2++) {
  t4=P1MAP_2(t1,t2);
  t8=P1DIM0_3(t1,t2);
  t10=P1DIM1_4(t1,t2);
  t14=P1MAP_2(t1,t2);
  s4(9,t2,0,t4,0,t2,0,t8,0,t10,0,t8,0,t14,0,0,0);
}
for(t2 = 0; t2 <= P1SIZE-1; t2++) {
  t4=P1MAP_2(t1,t2);
  t8=P1DIM0_3(t1,t2);
  t10=P1DIM1_4(t1,t2);
  t12=P1DIM0_3(t1,t2);
  t14=P1MAP_2(t1,t2);
  t16=P1DIM1_4(t1,t2);
  s5(10,t2,0,t4,0,t2,0,t8,0,t10,0,t12,0,t14,0,t16,0);
}
for(t2 = 0; t2 <= P1SIZE-1; t2++) {
  t4=P1MAP_2(t1,t2);
  t8=P1DIM0_3(t1,t2);
  t10=P1DIM1_4(t1,t2);
  t14=P1MAP_2(t1,t2);
  s6(12,t2,0,t4,0,t2,0,t8,0,t10,0,t8,0,t14,0,0,0);
}
for(t2 = 0; t2 <= P1SIZE-1; t2++) {
  t4=P1MAP_2(t1,t2);
  t8=P1DIM0_3(t1,t2);
  t10=P1DIM1_4(t1,t2);
  t12=P1DIM0_3(t1,t2);
  t14=P1MAP_2(t1,t2);
  s7(14,t2,0,t4,0,t2,0,t8,0,t10,0,t12,0,t14,0,0,0);
}
for(t2 = 0; t2 <= P1SIZE-1; t2++) {
  t4=P1MAP_2(t1,t2);
  t8=P1DIM0_3(t1,t2);
  t10=P1DIM1_4(t1,t2);
  t12=P1DIM0_3(t1,t2);
  t14=P1MAP_2(t1,t2);
  s8(16,t2,0,t4,0,t2,0,t8,0,t10,0,t12,0,t14,0,0,0);
}
for(t2 = 0; t2 <= NR-1; t2++) {
  for(t4 = t2+1; t4 <= NR; t4++) {
    s9(17,t2,0,t4,0,0,0,0,0,0,0,0,0,0,0,0,0);
  }
}
for(t2 = 0; t2 <= P1SIZE-1; t2++) {
  t4=P1MAP_2(t1,t2);
  t8=P1DIM0_3(t1,t2);
  t10=P1DIM1_4(t1,t2);
  t12=P1DIM0_3(t1,t2);
  t14=P1MAP_2(t1,t2);
  t16=P1DIM1_4(t1,t2);
  s10(18,t2,0,t4,0,t2,0,t8,0,t10,0,t12,0,t14,0,t16,0);
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
#undef s8
#undef s_8
#undef s9
#undef s_9
#undef s10
#undef s_10
#undef P1DIM0_3
#undef P1DIM1_4
#undef P1MAP_2
#undef col1_1
#undef row1_0
delete P1; 
