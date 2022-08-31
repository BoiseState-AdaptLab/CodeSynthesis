#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
#define A3DCOO(n,ii,jj,kk) EX_A3DCOO(n)
#define AM3DCOO(n1,i,j,k) EX_A3DMCOO(n1)
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
#undef s11
#undef s_11
#undef s12
#undef s_12
#undef s13
#undef s_13
#undef s14
#undef s_14
#undef s15
#undef s_15
#undef s16
#undef s_16
#undef s17
#undef s_17
#define s_0(n, ii, jj, kk)   P0->insert({ii, jj, kk}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_0(a1, a3, a5, a7);
#define s_1(x)   P0->sort() 
#define s1(__x0, __x1, __x2, __x3, __x4, __x5, __x6, __x7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_1(a1);
#define s_2(_no, _n, n, ii, jj, kk, n1)   mz3d(n1)=kk 
#define s2(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_2(a1, a3, a5, a7, a9, a11, a13);
#define s_3(_no, _n, n, ii, jj, kk, n1, i, j, k)   mz3d(n1)=k 
#define s3(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, a15, __x16, a17, __x18, a19, __x20)   s_3(a1, a3, a5, a7, a9, a11, a13, a15, a17, a19);
#define s_4(_no, _n, n, ii, jj, kk, n1)   mz3d(n1)=z3d(n) 
#define s4(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_4(a1, a3, a5, a7, a9, a11, a13);
#define s_5(_no, _n, n, ii, jj, kk, n1)   mz3d(n1) = max(mz3d(n1),) 
#define s5(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_5(a1, a3, a5, a7, a9, a11, a13);
#define s_6(_no, _n, n, ii, jj, kk, n1)   mz3d(n1) = min(mz3d(n1),NZ - 1) 
#define s6(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_6(a1, a3, a5, a7, a9, a11, a13);
#define s_7(_no, _n, n, ii, jj, kk, n1)   mrow3(n1)=ii 
#define s7(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_7(a1, a3, a5, a7, a9, a11, a13);
#define s_8(_no, _n, n, ii, jj, kk, n1, i)   mrow3(n1)=i 
#define s8(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, a15, __x16, __x17, __x18, __x19, __x20)   s_8(a1, a3, a5, a7, a9, a11, a13, a15);
#define s_9(_no, _n, n, ii, jj, kk, n1)   mrow3(n1)=row3d(n) 
#define s9(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_9(a1, a3, a5, a7, a9, a11, a13);
#define s_10(_no, _n, n, ii, jj, kk, n1)   mrow3(n1) = max(mrow3(n1),) 
#define s10(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_10(a1, a3, a5, a7, a9, a11, a13);
#define s_11(_no, _n, n, ii, jj, kk, n1)   mrow3(n1) = min(mrow3(n1),NR - 1) 
#define s11(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_11(a1, a3, a5, a7, a9, a11, a13);
#define s_12(_no, _n, n, ii, jj, kk, n1)   mcol3(n1)=jj 
#define s12(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_12(a1, a3, a5, a7, a9, a11, a13);
#define s_13(_no, _n, n, ii, jj, kk, n1, i, j)   mcol3(n1)=j 
#define s13(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, a15, __x16, a17, __x18, __x19, __x20)   s_13(a1, a3, a5, a7, a9, a11, a13, a15, a17);
#define s_14(_no, _n, n, ii, jj, kk, n1)   mcol3(n1)=col3d(n) 
#define s14(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_14(a1, a3, a5, a7, a9, a11, a13);
#define s_15(_no, _n, n, ii, jj, kk, n1)   mcol3(n1) = max(mcol3(n1),) 
#define s15(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_15(a1, a3, a5, a7, a9, a11, a13);
#define s_16(_no, _n, n, ii, jj, kk, n1)   mcol3(n1) = min(mcol3(n1),NC - 1) 
#define s16(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16, __x17, __x18, __x19, __x20)   s_16(a1, a3, a5, a7, a9, a11, a13);
#define s_17(_no, _n, n, ii, jj, kk, n1, i, j, k)   AM3DCOO(n1,i,j,k) = A3DCOO(n,ii,jj,kk ) 
#define s17(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, a15, __x16, a17, __x18, a19, __x20)   s_17(a1, a3, a5, a7, a9, a11, a13, a15, a17, a19);

#undef P0DIM0_4
#undef P0DIM1_5
#undef P0DIM2_6
#undef P0MAP_3
#undef col3d_1
#undef row3d_0
#undef z3d_2
#define P0SIZE P0->getSize()
#define P0DIM0(idx) P0->getDim(0,idx)
#define P0DIM0_4(__tv0, __tv1) P0DIM0(__tv1)
#define P0DIM1(idx) P0->getDim(1,idx)
#define P0DIM1_5(__tv0, __tv1) P0DIM1(__tv1)
#define P0DIM2(idx) P0->getDim(2,idx)
#define P0DIM2_6(__tv0, __tv1) P0DIM2(__tv1)
#define P0MAP(idx) P0->getMap(idx)
#define P0MAP_3(__tv0, __tv1) P0MAP(__tv1)
#define col3d_1(__tv0, __tv1) col3d(__tv1)
#define row3d_0(__tv0, __tv1) row3d(__tv1)
#define z3d_2(__tv0, __tv1) z3d(__tv1)

#define col3d(t0) EX_COL3D(t0)
#define mcol3(t0) EX_MCOL3(t0)
#define mrow3(t0) EX_MROW3(t0)
#define row3d(t0) EX_ROW3D(t0)
#define z3d(t0)   EX_Z3D(t0)
#define mz3d(t0)  EX_MZ3D(t0)

ReorderStream* P0= new ReorderStream(3);
ComparatorInt P0Comp = [&](const int a,const int b){
if(MORTON3D(row3d(a), col3d(a) , z3d(a)) <  MORTON3D(row3d(b), col3d(b) , z3d(b)))
return true;
return false;
}; 
P0->setComparator(P0Comp);

int t1 = 29; 
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
int t12 = 0; 
int t13 = 0; 
int t14 = 0; 
int t15 = 0; 
int t16 = 0; 
int t17 = 0; 
int t18 = 0; 
int t19 = 0; 
int t20 = 0; 
int t21 = 0; 

if (NC >= 1 && NZ >= 1 && NR >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row3d_0(t1,t2) >= 0 && NR >= row3d_0(t1,t2)+1 && NZ >= z3d_2(t1,t2)+1 && NC >= col3d_1(t1,t2)+1 && col3d_1(t1,t2) >= 0 && z3d_2(t1,t2) >= 0) {
      t4=row3d_0(t1,t2);
      t6=col3d_1(t1,t2);
      t8=z3d_2(t1,t2);
      s0(0,t2,0,t4,0,t6,0,t8,0,0,0,0,0,0,0,0,0,0,0,0,0);
    }
  }
}
auto start = std::chrono::high_resolution_clock::now();
s1(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
auto stop = std::chrono::high_resolution_clock::now();
for(t2 = 0; t2 <= P0SIZE-1; t2++) {
  t4=P0MAP_3(t1,t2);
  t8=P0DIM0_4(t1,t2);
  t10=P0DIM1_5(t1,t2);
  t12=P0DIM2_6(t1,t2);
  t14=P0MAP_3(t1,t2);
  s2(3,t2,0,t4,0,t2,0,t8,0,t10,0,t12,0,t14,0,0,0,0,0,0,0);
  s7(12,t2,0,t4,0,t2,0,t8,0,t10,0,t12,0,t14,0,0,0,0,0,0,0);
  s12(21,t2,0,t4,0,t2,0,t8,0,t10,0,t12,0,t14,0,0,0,0,0,0,0);
  s17(29,t2,0,t4,0,t2,0,t8,0,t10,0,t12,0,t4,0,t8,0,t10,0,t12,0);
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
#undef s11
#undef s_11
#undef s12
#undef s_12
#undef s13
#undef s_13
#undef s14
#undef s_14
#undef s15
#undef s_15
#undef s16
#undef s_16
#undef s17
#undef s_17
#undef P0DIM0_4
#undef P0DIM1_5
#undef P0DIM2_6
#undef P0MAP_3
#undef col3d_1
#undef row3d_0
#undef z3d_2
delete P0; 
