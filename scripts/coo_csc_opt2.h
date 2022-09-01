#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
#define ACOO(n,ii,jj) EX_ACOO(n)
#define ACSC(jj,k,ii) EX_ACSC(k)
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
#define s_0(n, ii, jj)   P1->insert({ii, jj}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6, __x7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_0(a1, a3, a5);
#define s_1(x)   P1->sort() 
#define s1(__x0, __x1, __x2, __x3, __x4, __x5, __x6, __x7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_1(a1);
#define s_2(_no, _n, n, ii, jj, jj1, k)   row4(k)=ii 
#define s2(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16)   s_2(a1, a3, a5, a7, a9, a11, a13);
#define s_3(_no, _n, n, ii, jj, jj1, k)   colptr(jj1) = min(colptr(jj1),k) 
#define s3(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16)   s_3(a1, a3, a5, a7, a9, a11, a13);
#define s_4(_no, _n, n, ii, jj, jj1, k)   colptr(jj1 + 1) = max(colptr(jj1 + 1),k + 1) 
#define s4(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, __x15, __x16)   s_4(a1, a3, a5, a7, a9, a11, a13);
#define s_5(e1, e2)   if ( not (colptr(e1) <= colptr(e2))){colptr(e2) = colptr(e1);} 
#define s5(__x0, a1, __x2, a3, __x4, __x5, __x6, __x7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_5(a1, a3);
#define s_6(_no, _n, n, ii, jj, jj1, k, ii1)   ACSC(jj1,k,ii1) = ACOO(n,ii,jj ) 
#define s6(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, a15, __x16)   s_6(a1, a3, a5, a7, a9, a11, a13, a15);

#undef P1DIM0_3
#undef P1DIM1_4
#undef P1MAP_2
#undef col1_1
#undef row1_0
#define P1SIZE P1->getSize()
#define P1DIM0(idx) P1->getDim(0,idx)
#define P1DIM0_3(__tv0, __tv1) P1DIM0(__tv1)
#define P1DIM1(idx) P1->getDim(1,idx)
#define P1DIM1_4(__tv0, __tv1) P1DIM1(__tv1)
#define P1MAP(idx) P1->getMap(idx)
#define P1MAP_2(__tv0, __tv1) P1MAP(__tv1)
#define col1_1(__tv0, __tv1) col1(__tv1)
#define row1_0(__tv0, __tv1) row1(__tv1)


#define col1(t0) EX_COL1(t0)
#define colptr(t0) EX_COLPTR(t0)
#define row1(t0) EX_ROW1(t0)
#define row4(t0) EX_ROW(t0)

ReorderStream* P1= new ReorderStream(2);
ComparatorInt P1Comp = [&](const int a,const int b){
if (P1->getDim(1,a) < P1->getDim(1,b) )    return true;
return false;
}; 
P1->setComparator(P1Comp);


int t1 = 10; 
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
auto start = std::chrono::high_resolution_clock::now();
for(t2 = 0; t2 <= P1SIZE-1; t2++) {
  t4=P1MAP_2(t1,t2);
  t8=P1DIM0_3(t1,t2);
  t10=P1DIM1_4(t1,t2);
  s2(2,t2,0,t4,0,t2,0,t8,0,t10,0,t10,0,t4,0,0,0);
  s3(3,t2,0,t4,0,t2,0,t8,0,t10,0,t10,0,t4,0,0,0);
  s4(6,t2,0,t4,0,t2,0,t8,0,t10,0,t10,0,t4,0,0,0);
  s6(10,t2,0,t4,0,t2,0,t8,0,t10,0,t10,0,t4,0,t8,0);
}
for(t2 = 0; t2 <= NR-1; t2++) {
    t4 = t2+1;
    s5(9,t2,0,t4,0,0,0,0,0,0,0,0,0,0,0,0,0);
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
#undef P1DIM0_3
#undef P1DIM1_4
#undef P1MAP_2
#undef col1_1
#undef row1_0
delete P1; 
