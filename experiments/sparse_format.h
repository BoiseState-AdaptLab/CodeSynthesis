#ifndef SPARSE_FORMAT_H
#define SPARSE_FORMAT_H
#include <string>
#include <string.h>
#include <stdlib.h>
#include <fstream>
struct coo_d {
   float * vals;
   int * cols;
   int * rows;
   int * zs;
   int nnz;
   int nr;
   int nc;
   int nz;
};
struct csc_d{
   float *vals;
   int * rows;
   int * cptr;
   int nnz;
};

struct Alto{
   unsigned long long * pos;
   float * vals;
   int nmode;
   unsigned long long* dims;
   unsigned long long alto_mask;
   unsigned long long * mode_masks;
   int nnz;
};

struct csr_d{
   float *vals;
   int * cols;
   int * rptr;
   int nnz;
};
struct matrix{
   int row;
   int col;
   float* vals;
};

struct dia_d {
   int ndiag;
   int * offsets; /*distance*/
   float * vals;
   int * idiag; /*ndiag * idiag = diag matrix values
	      *this is the dia matrix dimensions */
};
void read_sparse_coo (const std::string&, coo_d&);

void split(char *, char*,char* [10] );


#endif
