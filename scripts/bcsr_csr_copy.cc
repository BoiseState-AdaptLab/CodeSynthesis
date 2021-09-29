#include <iostream>

#define BR 2
#define BC 2
#define NR_BR 3
#define NR 6
#define s_0(ii, kk, jj, hr, hc, p, k)   ACSR[k] = ABCSR[p]; 
#define s0(ii, kk, jj, hr, hc, p, k)   s_0(ii, kk, jj, hr, hc, p, k);

#define bcol(t0) bcol[t0]
#define bcol_0(__tv0, __tv1) bcol(__tv1)
#define browptr(t0) browptr[t0]
#define browptr_1(__tv0) browptr(__tv0 + 1 )
#define col(t0) col[t0]
#define col_2(__tv0, __tv1, __tv2, __tv3, __tv4, __tv5, __tv6) col(__tv6)
#define col_inv(t0,t1) col_inv[t0][t1]
#define col_inv_1(__tv0, __tv1, __tv2, __tv3, __tv4) col_inv(BR * __tv0 + __tv3, BC * __tv2 + __tv4)
#define rowptr(t0) rowptr[t0]
#define rowptr1_0(__tv0, __tv1, __tv2, __tv3) rowptr(BR * __tv0 + __tv3)
#define rowptr1_1(__tv0, __tv1, __tv2, __tv3) rowptr(BR * __tv0 + __tv3 + 1)
#define rowptr_4_11(__tv0, __tv1, __tv2, __tv3) rowptr_4(__tv0, __tv1, __tv2, __tv3)
int main () {

int t1 = 0; 
int t2 = 0; 
int t3 = 0; 
int t4 = 0; 
int t5 = 0; 
int t6 = 0; 
int t7 = 0; 

int browptr[] = {0,1,2,3} ;
int bcol[] = { 0, 2, 0};
double ABCSR[] = {1,0,5,3,2,8,1,0,1,3,2,0};

double ACSR[] = {0,0,0,0,0,0,0,0,0};
int rowptr[] = {0,1,3,5,6,8,9};
int col[] ={0,0,1,4,5,4,0,1,0};
int col_inv[][6]={
	{0,-1,-1,-1,-1,-1 },
	{1,2,-1,-1,-1,-1 },
	{-1,-1,-1,-1,3,4 },
	{-1,-1,-1,-1,5,-1 },
	{6,7,-1,-1,-1,-1 },
	{8,-1,-1,-1,-1,-1 },

};

for(t1 = 0; t1 <= NR_BR-1; t1++) {
  for(t2 = browptr(t1); t2 <= browptr_1(t1)-1; t2++) {
    t3=bcol_0(t1,t2);
    for(t4 = 0; t4 <= BR; t4++) {
      int i_t = BR * t1 + t3;
      if (rowptr1_1(t1,t2,t3,t4) >= rowptr1_0(t1,t2,t3,t4)+1) {
	for(t5 = 0; t5 <= BC; t5++) {
	  if (rowptr1_1(t1,t2,t3,t4) >= col_inv_1(t1,t2,t3,t4,t5)+1 
			  && col_inv_1(t1,t2,t3,t4,t5) >= rowptr1_0(t1,t2,t3,t4)) {
            t7=col_inv_1(t1,t2,t3,t4,t5);
            if (col_2(t1,t2,t3,t4,t5,t6,t7) == BR*bcol_0(t1,t2)+t5) {
              s0(t1,t2,t3,t4,t5,BR*BC*t2+BC*t4+t5,t7);
            }
          }
        }
      }
    }
  }
}
std::cout << "{ ";
for (int i = 0 ; i < 9; i++){
   std::cout<< ACSR[i] << ",";
}
std::cout << "}\n";
}
