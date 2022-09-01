// Optimizations
// 1. Redundant Statement Elimination
// 2. Fusion
// 3. 

#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
#define ACSR(i,k) EX_ACSR(k)
#define ACSC(j,k) EX_ACSC(k)
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
#undef s6
#undef s_6
#undef s7
#undef s_7
#define s_0(i, k)   P0->insert({i, col2(k),col2(k)}) 
#define s0(__x0, a1, __x2, a3, __x4, __x5, __x6, __x7, __x8)   s_0(a1, a3);
#define s_1(i, k, tv2, tv3)   P1->insert({i, col2(k), k}) 
#define s1(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8)   s_1(a1, a3, a5, a7);
#define s_2(i, k, k1)   colptr(col2(k)) = min(colptr(col2(k)),k1) 
#define s2(__x0, a1, __x2, a3, tv4, __x5, __x6, __x7, __x8, __x9)   s_2(a1, a3, __x5);
#define s_3(i, k, k1)   colptr(col2(k) + 1) = max(colptr(col2(k) + 1),k1 + 1) 
#define s3(__x0, a1, __x2, a3, tv4, __x5, __x6, __x7, __x8, __x9)   s_3(a1, a3, __x5);
#define s_4(i, k)   row4(k)=i 
#define s4(__x0, a1, __x2, a3, tv4, __x5, __x6, __x7, __x8, __x9)   s_4(a1, a3);
#define s_5(i, k, j, k1)   row4(k1)=i 
#define s5(__x0, a1, __x2, a3, tv4, __x5, a5, __x7, a7, __x9)   s_5(a1, a3, a5, a7);
#define s_6(e1, e2)   if ( not (colptr(e1) <= colptr(e2))){colptr(e2) = colptr(e1);} 
#define s6(__x0, a1, __x2, a3, __x4, __x5, __x6, __x7, __x8)   s_6(a1, a3);
// TODO: bug in tuple variable names
#define s_7(i, k, j, k1)   ACSC(j,k1) = ACSR(i,k ) 
#define s7(__x0, a1, __x2, a3, tv4, __x5, a5, __x7, a7, __x9)   s_7(a1, a3, a5, a7);

#undef P1_4
#undef col2_2
#undef colptr_6
#undef colptr_7
#undef row4_5
#undef P0_3
#undef rowptr_0
#undef rowptr_1
#define P1_INV(d) P1->getInv(d)
#define P1(t0,t1) P1->get({t0,t1})
#define P1_4(__tv0, __tv1, __tv2, __tv3, __tv4) P1_INV(__tv3)[2]
#define col2(t0) EX_COL(t0)
#define col2_2(__tv0, __tv1, __tv2, __tv3) col2(__tv3)
#define colptr(t0) EX_COLPTR(t0)
#define colptr_6(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6) colptr(__tv6)
#define colptr_7(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6) colptr(__tv6 + 1)
#define row4(t0) EX_ROW(t0)
#define row4_5(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6, __tv7, __tv8) row4(__tv8)
#define P0(t0,t1) P0->get({t0,t1})
#define P0_3(__tv0, __tv1, __tv2, __tv3, __tv4) P0(__tv1, __tv4)
#define rowptr(t0) EX_ROWPTR(t0)
#define rowptr_0(__tv0, __tv1) rowptr(__tv1)
#define rowptr_1(__tv0, __tv1) rowptr(__tv1 + 1)

int t1 = 7; 
int t2 = 7; 
int t3 = 0; 
int t4 = 0; 
int t5 = 0; 
int t6 = 0; 
int t7 = 0; 
int t8 = 0; 
int t9 = 0; 
int t10 = 0; 

  // Assumed to be sorted so not timed
auto P1Comp = [&]( const std::vector<int>& a,const  std::vector<int>& b){
return a[1]*NR + a[0] < NR*b[1] + b[0];
};

auto P1Line = [&]( const std::vector<int>& a){
return a[1]*NR + a[0];
};
PermuteSimp<decltype(P1Comp)>* P1 = new PermuteSimp <decltype(P1Comp)>
(P1Comp);
for( t2 = 0; t2 <= NR-1; t2++) {
   for( t4 = rowptr_0(t1,t2); t4 <= rowptr_1(t1,t2)-1; t4++) {
      t8=col2_2(t1,t2,t3,t4);
      s1(1,t2,0,t4,0,t2,0,t8,0);
   }
}
P1->sort();
auto start = std::chrono::high_resolution_clock::now();
// Fused row, colptr and data copy
  for( t2 = 0; t2 <= NR-1; t2++) {
    for( t4 = rowptr_0(t1,t2); t4 <= rowptr_1(t1,t2)-1; t4++) {
      if (NC >= col2_2(t1,t2,t3,t4)+1 && col2_2(t1,t2,t3,t4) >= 0) {
        t5=col2_2(t1,t2,t3,t4);
        t7=t5;
        t9=P1_4(t1,t2,t3,t4,t5); 
	s3(3,t2,0,t4,t5,t9,0,0,0,0);
        s2(2,t2,0,t4,t5,t9,0,0,0,0);
        s7(7,t2,0,t4,t5,0,t7,0,t9,0);
        s5(4,t2,0,t4,t5,0,t7,0,t9,0);
      }
    }
  }
// Optimize UQ verification.
for(t2 = 0; t2 <= NR-1; t2++) {
   t4 = t2 + 1; 
   s6(6,t2,0,t4,0,0,0,0,0);
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
#undef P1_4
#undef col2_2
#undef colptr_6
#undef colptr_7
#undef row4_5
#undef P0_3
#undef rowptr_0
#undef rowptr_1
delete P1; 
