#include <synth.h> 
#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
auto offComp = [&](int a,int b){
// Manually add the comparator
if (a < b) return true;
return false;
}; 
// TODO: Fix this.
GrowthFunc<decltype(offComp)>* off = new GrowthFunc<decltype(offComp)>(offComp);
#define ACOO(n) EX_ACOO(n)
#define ADIA(id,dd,jj,kd) EX_ADIA(kd)
#undef s0
#undef s_0
#undef s1
#undef s_1
#define s_0(n)   off->insert(col1(n) - row1(n)) 
#define s0(__x0, a1, __x2, __x3, __x4, __x5, __x6, __x7, __x8, __x9, __x10)   s_0(a1);
#define s_1(n, id, dd, jj, kd)   ADIA(id,dd,jj,kd) = ACOO(n ) 
#define s1(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10)   s_1(a1, a3, a5, a7, a9);

#undef col1_0
#undef off_2
#undef row1_1
#define col1(t0) EX_COL1(t0)
#define col1_0(__tv0, __tv1) col1(__tv1)
// Added this.
#define off(t0) off->get(t0)
#define off_2(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5) off(__tv5)
#define row1(t0) EX_ROW1(t0)
#define row1_1(__tv0, __tv1) row1(__tv1)

uint64_t t1 = 1; 
uint64_t t2 = 0; 
uint64_t t3 = 0; 
uint64_t t4 = 0; 
uint64_t t5 = 0; 
uint64_t t6 = 0; 
uint64_t t7 = 0; 
uint64_t t8 = 0; 
uint64_t t9 = 0; 
uint64_t t10 = 0; 
uint64_t t11 = 0; 

if (NR >= 1 && NC >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (row1_1(t1,t2) >= 0 && NC >= col1_0(t1,t2)+1 && NR >= row1_1(t1,t2)+1 && col1_0(t1,t2) >= 0) {
      s0(0,t2,0,0,0,0,0,0,0,0,0);
    }
  }
  // Added this.
  // d < ND
  ND = off->getSize(); 
  // Added this
  off->sort();
  values.resize(ND * NR,INT_MIN);
  if (ND >= 1) {
    for(t2 = 0; t2 <= NNZ-1; t2++) {
      if (col1_0(t1,t2) >= 0 && NR >= row1_1(t1,t2)+1 && row1_1(t1,t2) >= 0 && NC >= col1_0(t1,t2)+1) {
        t4=row1_1(t1,t2);
        t8=col1_0(t1,t2);
        // Replace with binary search
	t6=off->search(t8-t4);	
	// Changed 99 to ND
        t10=ND*row1_1(t1,t2)+t6;
	s1(1,t2,0,t4,0,t6,0,t8,0,t10,0);
      }
    }
  }
}

#undef s0
#undef s_0
#undef s1
#undef s_1
#undef col1_0
#undef off_2
#undef row1_1
// Removed deletion of <off>
