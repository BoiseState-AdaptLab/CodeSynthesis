#include "sparse_format.h"
#include "bitops.hpp"
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <iostream>
#include <assert.h>
// This mask is from the paper 4 x 8 x 2
#define MASK1 0b001010
#define MASK2 0b110100
#define MASK3 0b000001


// This mask is for NELL-2
/*
#define MASK1 0b0010010010010010010010010010010010010010010
#define MASK2 0b0001001001001001001001001001001001001001001
#define MASK3 0b1100100100100100100100100100100100100100100
*/
#define R 200
#define SEED 431890943221
#define ITER 10

void morton_sort(coo_d* coo);

unsigned long long linearize(int i, int j , int k){
    return  pdep((long long unsigned int) i ,(long long unsigned int) MASK1) |
		pdep((long long unsigned int) j,(long long unsigned int)MASK2) |
		pdep((long long unsigned int) k,(long long unsigned int)MASK3);

}

void linearization_test(){
    assert (2 == linearize(1,0,0));
    assert (15 == linearize(3,1,1));
    assert (20 == linearize(0,3,0));
    assert (25 == linearize(2,2,1));
    assert (42 == linearize(3,4,0));
    assert (51 == linearize(1,6,1));
     
}

void test_build_alto_paper (Alto* alto){
    assert(6 == alto->nnz);
    assert(2 == alto->pos[0]);
    assert (15 == alto->pos[1]);
    assert (20 == alto->pos[2]);
    assert (25 == alto->pos[3]);
    assert (42 == alto->pos[4]);
    assert (51 == alto->pos[5]);
}

void build_alto(Alto* alto, coo_d* coo){
    alto->pos = new unsigned long long [coo->nnz]();
    alto->vals = new float [coo->nnz] ();
    alto->nnz = coo->nnz;
    alto->nr =  coo->nr;
    alto->nc =  coo->nc;
    alto->nz =  coo->nz;
    for(int n = 0 ; n < alto->nnz ;n++ ){
        unsigned long long pos = 0;
	//encode i coordinate
	pos= linearize(coo->rows[n],coo->cols[n],coo->zs[n]);
	// Perform insertion sort
	int i = n - 1;
	while (i >= 0 && pos < alto->pos[i]){
	        alto->pos[i+1] = alto->pos[i];
	        alto->vals[i+1] = alto->vals[i];
	        i--;
	}
	alto->pos[i+1]  = pos;
	alto->vals[i+1] = coo->vals[n];
    }
}

void alto_mttkrp_trns(Alto &alto, matrix &A, matrix &B, matrix &C){
   for(int r = 0 ; r < R; r++){
      for( int n =0; n < alto.nnz; n++){
          int i = pext((long long unsigned int)alto.pos[n],(long long unsigned int)MASK1);
          int j = pext((long long unsigned int)alto.pos[n],(long long unsigned int)MASK2);
          int k = pext((long long unsigned int)alto.pos[n],(long long unsigned int)MASK3);
          A.vals[r * R + i] += alto.vals[n] * B.vals[r*R + j] * C.vals[r*R+k];
      }
   }
}

void alto_mttkrp(Alto &alto, matrix &A, matrix &B, matrix &C){
   for( int n =0; n < alto.nnz; n++){
      for(int r = 0 ; r < R; r++){
          int i = pext((long long unsigned int)alto.pos[n],(long long unsigned int)MASK1);
          int j = pext((long long unsigned int)alto.pos[n],(long long unsigned int)MASK2);
          int k = pext((long long unsigned int)alto.pos[n],(long long unsigned int)MASK3);
          A.vals[i * A.row + r] += alto.vals[n] * B.vals[j*B.row + r] * C.vals[k * C.row + r];
      }
   }
}

void alto_mttkrp_parallel_1(Alto &alto, matrix &A, matrix &B, matrix &C){
   #pragma omp parallel for collapse(2)
   for(int r = 0 ; r < R; r++){
      for( int n =0; n < alto.nnz; n++){
          int i = pext((long long unsigned int)alto.pos[n],(long long unsigned int)MASK1);
          int j = pext((long long unsigned int)alto.pos[n],(long long unsigned int)MASK2);
          int k = pext((long long unsigned int)alto.pos[n],(long long unsigned int)MASK3);
          A.vals[r * R + i] += alto.vals[n] * B.vals[r*R + j] * C.vals[r*R+k];
      }
   }
}

void alto_mttkrp_parallel_2(Alto &alto, matrix &A, matrix &B, matrix &C, int L){
   int part = alto.nnz / L;
   int nnz = alto.nnz;
   unsigned long long alto_bitmask_size = MASK1 | MASK2 | MASK3;
   float* temp  = new float [alto_bitmask_size * R] ();
   #pragma omp parallel for
   for(int l = 0 ; l < L ; l++){
      int rL= l * part ;
      int uL = l == L -1? nnz : rL + part;
      for(int ll = rL; ll < uL ;ll++){
          int i = pext((long long unsigned int)alto.pos[ll],(long long unsigned int)MASK1);
          int j = pext((long long unsigned int)alto.pos[ll],(long long unsigned int)MASK2);
          int k = pext((long long unsigned int)alto.pos[ll],(long long unsigned int)MASK3);
          for(int r = 0 ; r < R; r++){
             temp[ (i - alto.pos[rL])* alto_bitmask_size +  r] += alto.vals[ll] * B.vals[j*B.row + r] * C.vals[k * C.row + r];

	     #pragma omp critical
             A.vals[i * A.row + r] += alto.vals[ll] * B.vals[j*B.row + r] * C.vals[k * C.row + r];
	  }
      }
   }
   #pragma omp parallel for
   for(int b = 0; b < alto.nr; b++){
      for(int l = 0; l < L; l++){
         
         int rL= l * part ;
         int uL = l == L -1? nnz : rL + part;
	 //Check if b is an element of current l
         for(int ll = rL; ll < uL ;ll++){
            int i = pext((long long unsigned int)alto.pos[ll],(long long unsigned int)MASK1);
            if ( i == b){
		// Accumulation Code
                for(int r = 0 ; r < R; r++){
	           A.vals[i * A.row + r] += temp[ (i - alto.pos[rL])* alto_bitmask_size +  r];
	      	}
	    }
	 } 
      }    
   }
}

void generate_random_matrix(matrix &randMatrix){
   if (randMatrix.vals) 
       delete randMatrix.vals;
   randMatrix.vals = new float[randMatrix.row * randMatrix.col];

   for(int i = 0; i <randMatrix.row * randMatrix.col; i++){
       randMatrix.vals[i] = rand();
   }
}

int main ( int ac, char** argv){
   srand(SEED);
   
   linearization_test();
   std::cerr << "Linearization passed! \n";
   // Tensor from paper
   coo_d tens_paper;
   tens_paper.nnz = 6;
   tens_paper.nr = 4;
   tens_paper.nc = 8;
   tens_paper.nz = 2;
   int  rows []=  {0,1,1,2,3,3};
   int  cols []=  {3,0,6,2,1,4};
   int  zs []  =  {0,0,1,1,1,0};
   float  vals []= {5,6,9,3,20,2}; 
 
   tens_paper.rows = &rows[0];
   tens_paper.cols = &cols[0];
   tens_paper.zs =   &zs[0];
   
   Alto alto_paper;
   std::cerr << "Building alto\n";
   build_alto(&alto_paper,&tens_paper);
   std::cerr << "Finished building alto \n";
   
   test_build_alto_paper(&alto_paper);
   std::cerr << "Build Alto Test passed\n";
      
   /* // This code is for nell-2, uncomment this code 
    * // to when using nell2 tensor and also uncomment 
    * // the mask information.
   std::string file_name(argv[1]);
   // Nell2
   coo_d nell_coo;
   std::cout << "Reading from file: "<< file_name << " \n";
   read_sparse_coo(file_name,nell_coo);  
   std::cout << "Finished reading from file\n";
   Alto alto;

   std::cout << "Building alto \n";
   build_alto(&alto,&nell_coo);
   std::cout << "Finished alto \n";
   matrix B;
   matrix A;
   matrix C;
   B.row = nell_coo.nc;
   B.col = R;

   A.row = nell_coo.nr;
   A.col = R;
   
   A.row = nell_coo.nz;
   A.col = R;
   std::cout << "Generating Factor Matrices \n";
   generate_random_matrix(A);
   generate_random_matrix(B);
   generate_random_matrix(C);
  
   clock_t start, end;
   double cpu_time_used = 0;
   std::cout << "Running Mode-1 MTTKRP: \n"; 
   for (int i  = 0; i < ITER ; i++){
      start = clock();
      // Mode1 MTTKRP
      alto_mttkrp(alto,A,B,C);
      end = clock();
      cpu_time_used+=((double) (end - start)) / CLOCKS_PER_SEC;
   }
   double average_time = cpu_time_used / ITER;
   std::cout << "Average Time for Mode-1 MTTKRP: "<< 
	   average_time <<"\n";
   
   std::cout << "Generating Factor Matrices \n";
   generate_random_matrix(A);
   generate_random_matrix(B);
   generate_random_matrix(C);

   cpu_time_used = 0;
   std::cout << "Running Mode-1 MTTRK-reordered and transpose: \n";
   for (int i  = 0; i < ITER ; i++){
      start = clock();
      // Mode1 MTTKRP
      alto_mttkrp_trns(alto,A,B,C);
      end = clock();
      cpu_time_used+=((double) (end - start)) / CLOCKS_PER_SEC;
   }
   average_time = cpu_time_used / ITER;
   std::cout << "Average Time for Mode-1 MTTKRP, re-ordered & transposed: "<< 
	   average_time <<"\n";*/
}