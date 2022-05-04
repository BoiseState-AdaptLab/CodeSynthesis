#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
	auto P0Comp = [](const std::vector<int>& a, const std::vector<int>& b){
	return MORTON(a[0], a[1]) < MORTON(b[0], b[1]);
	};
	auto P0Lin = [](const std::vector<int>& a){
	   return MORTON(a[0],a[1]);
	};
	Permute<int,decltype(P0Lin),decltype(P0Comp)> * P0 = 
	new Permute<int,decltype(P0Lin),decltype(P0Comp)>(P0Lin,P0Comp,22);
#define ACOO(n) EX_ACOO(n)
//TODO n generated for RHS, fix!!!
#define AMCOO(n1) EX_AMCOO(n1)
#undef s0
#undef s_0
#undef s1
#undef s_1
#undef s2
#undef s_2
#undef s3
#undef s_3
#define s_0(n, tv1, tv2)   P0->insert({row1(n), col1(n),n}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6)   s_0(a1, a3, a5);
#define s_1(n)   col3(n)=P0_INV(n)[1] 
#define s1(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_1(a1);
#define s_2(n)   row3(n)=P0_INV(n)[0] 
#define s2(__x0, a1, tv2, tv3, __x4, __x5, __x6, __x7, __x8)   s_2(a1);
#define s_3(n, n1)   AMCOO(n) = ACOO(n1 ) 
#define s3(__x0, a1, tv2, tv3, __x4, a3, __x6, __x7, __x8)   s_3(a1, a3);


#undef P0_2
#undef col1_1
#undef col3_3
#undef row1_0
#undef row3_4
#define P0_INV(t0) P0->getInv(t0)
#define P0(t0,t1) P0->get({t0,t1})
#define P0_2(__tv0, __tv1, __tv2, __tv3) P0_INV(__tv1)[2]
#define col1(t0) EX_COL1(t0)
#define col1_1(__tv0, __tv1) col1(__tv1)
#define col3(t0) EX_COL3(t0)
#define col3_3(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) col3(__tv5)
#define row1(t0) EX_ROW1(t0)
#define row1_0(__tv0, __tv1) row1(__tv1)
#define row3(t0) EX_ROW3(t0)
#define row3_4(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) row3(__tv5)

int t1 = 3; 
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
      s1(1,t2,t3,t4,0,0,0,0,0);
      s2(2,t2,t3,t4,0,0,0,0,0);
      t6=P0_2(t1,t2,t3,t4);
      s3(3,t2,t3,t4,0,t6,0,0,0);
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
#undef P0_2
#undef col1_1
#undef col3_3
#undef row1_0
#undef row3_4
delete P0; 
