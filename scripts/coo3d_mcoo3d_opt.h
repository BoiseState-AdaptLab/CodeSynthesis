#define min(a,b) a < b ? a : b
#define max(a,b) a > b ? a: b
auto P0Comp = [&](std::vector<int>& a,std::vector<int>& b){
if(-MORTON3D(a[0], a[1], a[2]) + MORTON3D(b[0], b[1], b[2]) - 1>= 0)
   return true;
return false;
}; 
PermuteSimp<decltype(P0Comp)>* P0 = new PermuteSimp <decltype(P0Comp)>(P0Comp);
#define A3DCOO(n,ii,jj,kk) EX_A3DCOO(n)
#define A3DMCOO(n1,i,j,k)  EX_A3DMCOO(n1)
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
#define s_0(n, ii, jj, kk)   P0->insert({ii, jj, kk, n}) 
#define s0(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_0(a1, a3, a5, a7);
#define s_1(t1)   P0->sort() 
#define s1(__x0, __x1, __x2, __x3, __x4, __x5, __x6, __x7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_1(a1);
#define s_2(n, ii, jj, kk)   mz3d(P0(ii, jj, kk))=kk 
#define s2(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_2(a1, a3, a5, a7);
#define s_3(n, ii, jj, kk)   mz3d(n)=z3d(P0_INV(n)[3]) 
#define s3(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_3(a1, a3, a5, a7);
#define s_4(n, ii, jj, kk)   mrow3(P0(ii, jj, kk))=ii 
#define s4(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_4(a1, a3, a5, a7);
#define s_5(n, ii, jj, kk)   mrow3(n)=row3d(P0_INV(n)[3]) 
#define s5(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_5(a1, a3, a5, a7);
#define s_6(n, ii, jj, kk)   mcol3(P0(ii, jj, kk))=jj 
#define s6(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_6(a1, a3, a5, a7);
#define s_7(n, ii, jj, kk)   mcol3(n)=col3d(P0_INV(n)[3]) 
#define s7(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, __x9, __x10, __x11, __x12, __x13, __x14, __x15, __x16)   s_7(a1, a3, a5, a7);
#define s_8(n, ii, jj, kk, n1, i, j, k)   A3DMCOO(n1,i,j,k) = A3DCOO(n,ii,jj,kk ) 
#define s8(__x0, a1, __x2, a3, __x4, a5, __x6, a7, __x8, a9, __x10, a11, __x12, a13, __x14, a15, __x16)   s_8(a1, a3, a5, a7, a9, a11, a13, a15);

#undef P0_3
#undef col3d_6
#undef mcol3_8
#undef mrow3_7
#undef mz3d_9
#undef row3d_5
#undef z3d_4
#undef col3d_1
#undef row3d_0
#undef z3d_2
#define P0_INV(t0) P0->getInv(t0)
#define P0(t0,t1,t2) P0->get({t0,t1,t2})
#define P0_3(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6, __tv7) P0_INV(__tv1)[3]
#define col3d(t0) EX_COL3D(t0)
#define col3d_6(__tv0, __tv1) col3d(__tv1)
#define mcol3(t0) EX_MCOL3(t0)
#define mcol3_8(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6, __tv7, __tv8, __tv9) mcol3(__tv9)
#define mrow3(t0) EX_MROW3(t0)
#define mrow3_7(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6, __tv7, __tv8, __tv9) mrow3(__tv9)
#define mz3d(t0)  EX_MZ3D(t0)
#define mz3d_9(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6, __tv7, __tv8, __tv9) mz3d(__tv9)
#define row3d(t0) EX_ROW3D(t0)
#define row3d_5(__tv0, __tv1) row3d(__tv1)
#define z3d(t0)   EX_Z3D(t0)
#define z3d_4(__tv0, __tv1) z3d(__tv1)
#define col3d(t0) EX_COL3D(t0)
#define col3d_1(__tv0, __tv1) col3d(__tv1)
#define row3d(t0) EX_ROW3D(t0)
#define row3d_0(__tv0, __tv1) row3d(__tv1)
#define z3d(t0)   EX_Z3D(t0)
#define z3d_2(__tv0, __tv1) z3d(__tv1)

int t1 = 8; 
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

if (NR >= 1 && NZ >= 1 && NC >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NC >= col3d_1(t1,t2)+1 && row3d_0(t1,t2) >= 0 && z3d_2(t1,t2) >= 0 && NR >= row3d_0(t1,t2)+1 && NZ >= z3d_2(t1,t2)+1 && col3d_1(t1,t2) >= 0) {
      t4=row3d_0(t1,t2);
      t6=col3d_1(t1,t2);
      t8=z3d_2(t1,t2);
      s0(0,t2,0,t4,0,t6,0,t8,0,0,0,0,0,0,0,0,0);
    }
  }

}

// time reordering step
auto start = std::chrono::high_resolution_clock::now();
s1(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

auto stop = std::chrono::high_resolution_clock::now();
if (NR >= 1 && NZ >= 1 && NC >= 1) {
  for(t2 = 0; t2 <= NNZ-1; t2++) {
    if (NZ >= z3d_2(t1,t2)+1 && z3d_2(t1,t2) >= 0 && NC >= col3d_1(t1,t2)+1 && NR >= row3d_0(t1,t2)+1 && col3d_1(t1,t2) >= 0 && row3d_0(t1,t2) >= 0) {
      t4=row3d_0(t1,t2);
      t6=col3d_1(t1,t2);
      t8=z3d_2(t1,t2);

      t10=P0_3(t1,t2,t3,t4,t5,t6,t7,t8);
      s3(2,t2,0,t4,0,t6,0,t8,0,0,0,0,0,0,0,0,0);
      s5(2,t2,0,t4,0,t6,0,t8,0,0,0,0,0,0,0,0,0);
      s7(2,t2,0,t4,0,t6,0,t8,0,0,0,0,0,0,0,0,0);
      t14=col3d_6(t1,t2);
      s8(8,t2,0,t4,0,t6,0,t8,0,t10,0,t4,0,t14,0,t8,0);
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
#undef s8
#undef s_8
#undef P0_3
#undef col3d_6
#undef mcol3_8
#undef mrow3_7
#undef mz3d_9
#undef row3d_5
#undef z3d_4
#undef col3d_1
#undef row3d_0
#undef z3d_2
delete P0; 
