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

typedef unsigned long long LIType;
typedef unsigned long long LIT;
typedef float FType;

// This mask is for NELL-2
/*
#define MASK1 0b0010010010010010010010010010010010010010010
#define MASK2 0b0001001001001001001001001001001001001001001
#define MASK3 0b1100100100100100100100100100100100100100100
*/
#define R 200
#define SEED 431890943221
#define ITER 10


struct MPair {
    int mode;
    int bits;
};
typedef enum PackOrder_ { LSB_FIRST, MSB_FIRST } PackOrder;
typedef enum ModeOrder_ { SHORT_FIRST, LONG_FIRST, NATURAL } ModeOrder;

// Achieving alto_bits_min requires packing/compression.
static inline void
setup_packed_alto(Alto* at, PackOrder po, ModeOrder mo)
{
    LIT ALTO_MASKS[MAX_NUM_MODES] = {}; //initialized to zeros by default

    int nmode = at->nmode;
    int alto_bits_min = 0, alto_bits_max = 0;
    LIT alto_mask = 0;
    int max_num_bits = 0, min_num_bits = sizeof(IType) * 8;
    
    MPair* mode_bits = (MPair*)malloc(nmode * sizeof(MPair));
    assert(mode_bits);

    //Initial mode values.
    for (int n = 0; n < nmode; ++n) {
        int mbits = (sizeof(IType) * 8) - clz(at->dims[n] - 1);
        mode_bits[n].mode = n;
        mode_bits[n].bits = mbits;
        alto_bits_min += mbits;
        max_num_bits = std::max(max_num_bits, mbits);
        min_num_bits = std::min(min_num_bits, mbits);
        printf("num_bits for mode-%d=%d\n", n + 1, mbits);
    }
    
#ifdef ALT_PEXT
    //Simple prefix sum
    at->mode_pos[0] = 0;
    for (int n = 1; n < nmode; ++n) {
        at->mode_pos[n] = at->mode_pos[n-1] + mode_bits[n-1].bits;
    }
#endif
    
    alto_bits_max = max_num_bits * nmode;
    //printf("range of mode bits=[%d %d]\n", min_num_bits, max_num_bits);
    printf("alto_bits_min=%d, alto_bits_max=%d\n", alto_bits_min, alto_bits_max);

    assert(alto_bits_min <= ((int)sizeof(LIT) * 8));

    //Assuming we use a power-2 data type for ALTO_idx with a minimum size of a byte
    //int alto_bits = pow(2, (sizeof(int) * 8) - __builtin_clz(alto_bits_min));
    int alto_bits = (int)0x1 << std::max<int>(3, (sizeof(int) * 8) - __builtin_clz(alto_bits_min));
    printf("alto_bits=%d\n", alto_bits);

    double alto_storage = 0;
    alto_storage = at->nnz * (sizeof(FType) + sizeof(LIT));
    printf("Alto format storage:    %g Bytes\n", alto_storage);
    
    alto_storage = at->nnz * (sizeof(FType) + (alto_bits >> 3));
    printf("Alto-power-2 format storage:    %g Bytes\n", alto_storage);

    alto_storage = at->nnz * (sizeof(FType) + (alto_bits_min >> 3));
    printf("Alto-opt format storage:    %g Bytes\n", alto_storage);
    
    {//Dilation & shifting.
        int level = 0, shift = 0, inc = 1;

        //Sort modes, if needed.
        if (mo == SHORT_FIRST)
            std::sort(mode_bits, mode_bits + nmode, [](auto& a, auto& b) { return a.bits < b.bits; });
        else if(mo == LONG_FIRST)
            std::sort(mode_bits, mode_bits + nmode, [](auto& a, auto& b) { return a.bits > b.bits; });
          
        if (po == MSB_FIRST) {
            shift = alto_bits_min - 1;
            inc = -1;
        }
        
        bool done;
        do {
            done = true;

            for (int n = 0; n < nmode; ++n) {
                if (level < mode_bits[n].bits) {
                    ALTO_MASKS[mode_bits[n].mode] |= (LIT)0x1 << shift;
                    shift += inc;
                    done = false;
                }
            }
            ++level;
        } while (!done);
        
        assert(level == (max_num_bits+1));
        assert(po == MSB_FIRST ? (shift == -1) : (shift == alto_bits_min));
    }

    for (int n = 0; n < nmode; ++n) {
        at->mode_masks[n] = ALTO_MASKS[n];
        alto_mask |= ALTO_MASKS[n];
#ifdef ALTO_DEBUG        
        printf("ALTO_MASKS[%d] = 0x%llx\n", n, ALTO_MASKS[n]);
#endif
    }
    at->alto_mask = alto_mask;
#ifdef ALTO_DEBUG
    printf("alto_mask = 0x%llx\n", alto_mask);
#endif
    free(mode_bits);
}

void morton_sort(coo_d* coo);
void free_alto (Alto *alto){
   delete [] alto->pos;
   delete [] alto->vals;
}

void free_coo (coo_d *coo){
   delete [] coo->rows;
   delete [] coo->cols;
   delete [] coo->zs;
   delete [] coo->vals;
}

unsigned long long linearize(int i, int j , int k, unsigned long long* mode_masks){
    return  pdep((long long unsigned int) i ,(long long unsigned int)mode_masks[0]) |
		pdep((long long unsigned int) j,(long long unsigned int)mode_masks[1]) |
		pdep((long long unsigned int) k,(long long unsigned int)mode_masks[2]);

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

    assert(6 == alto->vals[0]);
    assert (20 == alto->vals[1]);
    assert (5 == alto->vals[2]);
    assert (3 == alto->vals[3]);
    assert (2 == alto->vals[4]);
    assert (9 == alto->vals[5]);
}

void build_alto(Alto* alto, coo_d* coo){
    alto->pos = new unsigned long long [coo->nnz]();
    alto->vals = new float [coo->nnz] ();
    alto->nnz = coo->nnz;
    alto->nmode = 3;
    alto->mode_masks = new unsigned long long[3] ();
    alto->dims[0] = coo->nr;
    alto->dims[1] = coo->nc;
    alto->dims[2] = coo->nz;
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
          int i = pext((long long unsigned int)alto.pos[n],(long long unsigned int)alto.mode_masks[0]);
          int j = pext((long long unsigned int)alto.pos[n],(long long unsigned int)alto.mode_masks[1]);
          int k = pext((long long unsigned int)alto.pos[n],(long long unsigned int)alto.mode_masks[2]);
          A.vals[r * R + i] += alto.vals[n] * B.vals[r*R + j] * C.vals[r*R+k];
      }
   }
}

void alto_mttkrp(Alto &alto, matrix &A, matrix &B, matrix &C){
   for( int n =0; n < alto.nnz; n++){
      for(int r = 0 ; r < R; r++){
          int i = pext((long long unsigned int)alto.pos[n],(long long unsigned int)alto.mode_masks[0]);
          int j = pext((long long unsigned int)alto.pos[n],(long long unsigned int)alto.mode_masks[1]);
          int k = pext((long long unsigned int)alto.pos[n],(long long unsigned int)alto.mode_masks[2]);
          A.vals[i * A.row + r] += alto.vals[n] * B.vals[j*B.row + r] * C.vals[k * C.row + r];
      }
   }
}

void alto_mttkrp_parallel_1(Alto &alto, matrix &A, matrix &B, matrix &C){
   #pragma omp parallel for collapse(2)
   for(int r = 0 ; r < R; r++){
      for( int n =0; n < alto.nnz; n++){
          int i = pext((long long unsigned int)alto.pos[n],(long long unsigned int)alto.mode_masks[0]);
          int j = pext((long long unsigned int)alto.pos[n],(long long unsigned int)alto.mode_masks[1]);
          int k = pext((long long unsigned int)alto.pos[n],(long long unsigned int)alto.mode_masks[2]);
          A.vals[r * R + i] += alto.vals[n] * B.vals[r*R + j] * C.vals[r*R+k];
      }
   }
}

void alto_mttkrp_parallel_2(Alto &alto, matrix &A, matrix &B, matrix &C, int L){
   int part = alto.nnz / L;
   int nnz = alto.nnz;
   //Todo work on this part.
   unsigned long long alto_bitmask_size = 0;
   float* temp  = new float [alto_bitmask_size * R] ();
   #pragma omp parallel for
   for(int l = 0 ; l < L ; l++){
      int rL= l * part ;
      int uL = l == L -1? nnz : rL + part;
      for(int ll = rL; ll < uL ;ll++){
          int i = pext((long long unsigned int)alto.pos[ll],(long long unsigned int)alto.mode_masks[0]);
          int j = pext((long long unsigned int)alto.pos[ll],(long long unsigned int)alto.mode_masks[1]);
          int k = pext((long long unsigned int)alto.pos[ll],(long long unsigned int)alto.mode_masks[2]);
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
          int i = pext((long long unsigned int)alto.pos[ll],(long long unsigned int)alto.mode_masks[0]);
            if ( i == b){
		// Accumulation Code
                for(int r = 0 ; r < R; r++){
	           A.vals[i * A.row + r] += temp[ (i - alto.pos[rL])* alto_bitmask_size +  r];
	      	}
	    }
	 } 
      }    
   }

   delete temp;
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

   int*  rows=  new int [6] {0,1,1,2,3,3};
   int*  cols= new int [6]  {3,0,6,2,1,4};
   int*  zs  =  new int [6] {0,0,1,1,1,0};
   float*  vals = new float [6] {5,6,9,3,20,2}; 
   
    
   tens_paper.rows = rows;
   tens_paper.cols = cols;
   tens_paper.zs =   zs;
   tens_paper.vals = vals;


   Alto alto_paper;
   std::cerr << "Building alto\n";
   build_alto(&alto_paper,&tens_paper);
   std::cerr << "Finished building alto \n";
   
   std::cerr << "Create Masks\n";
   setup_packed_alto(&alto, LSB_FIRST, SHORT_FIRST);  

   // check if masks is correct
   assert (MASK1,alto.mode_masks[0]);
   assert (MASK2,alto.mode_masks[1]);
   assert (MASK3,alto.mode_masks[2]);


   test_build_alto_paper(&alto_paper);
   std::cerr << "Build Alto Test passed\n";
   
   free_alto(&alto_paper);
   free_coo(&tens_paper);  

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
