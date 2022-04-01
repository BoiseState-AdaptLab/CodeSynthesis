// Parts of this code are taken from
// https://github.com/llvm/llvm-project/blob/b682616d1fd1263b303985b9f930c1760033af2c/mlir/lib/ExecutionEngine/SparseTensorUtils.cpp
// Which is part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// These sections have been marked.

#include <cctype>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <synth.h>

/// Taken from LLVM SparseTnesorUitls (see top for details).
static constexpr int kColWidth = 1025;

/// Taken from LLVM SparseTnesorUitls (see top for details).
///
/// Helper to convert string to lower case.
static char *toLower(char *token) {
    for (char *c = token; *c; c++)
        *c = tolower(*c);
    return token;
}

/// Taken from LLVM SparseTnesorUitls (see top for details).
///
/// Read the MME header of a general sparse matrix of type real.
static void readMMEHeader(FILE *file, char *filename, char *line,
                          uint64_t *idata, bool *isSymmetric) {
    char header[64];
    char object[64];
    char format[64];
    char field[64];
    char symmetry[64];
    // Read header line.
    if (fscanf(file, "%63s %63s %63s %63s %63s\n", header, object, format, field,
               symmetry) != 5) {
        fprintf(stderr, "Corrupt header in %s\n", filename);
        exit(1);
    }
    *isSymmetric = (strcmp(toLower(symmetry), "symmetric") == 0);
    // Make sure this is a general sparse matrix.
    if (strcmp(toLower(header), "%%matrixmarket") ||
        strcmp(toLower(object), "matrix") ||
        strcmp(toLower(format), "coordinate") || strcmp(toLower(field), "real") ||
        (strcmp(toLower(symmetry), "general") && !(*isSymmetric))) {
        fprintf(stderr,
                "Cannot find a general sparse matrix with type real in %s\n",
                filename);
        exit(1);
    }
    // Skip comments.
    while (true) {
        if (!fgets(line, kColWidth, file)) {
            fprintf(stderr, "Cannot find data in %s\n", filename);
            exit(1);
        }
        if (line[0] != '%')
            break;
    }
    // Next line contains M N NNZ.
    idata[0] = 2; // rank
    if (sscanf(line, "%" PRIu64 "%" PRIu64 "%" PRIu64 "\n", idata + 2, idata + 3,
            idata + 1) != 3) {
        fprintf(stderr, "Cannot find size in %s\n", filename);
        exit(1);
    }
}

// Print COO data in roughly matrix marget exchanged format. The output is not a correct mtx format just close, this
// function is intended for debugging purposes.
void printCOOAsMTX(uint64_t nnz,
                   uint64_t rank,
                   const std::vector<double> &values,
                   const std::vector<std::vector<uint64_t>> &coord) {
    for (uint64_t k = 0; k < nnz; k++) {
        for (uint64_t r = 0; r < rank; r++) {
            if (r == rank-1) {
                // MTX uses 1 based indexing, hence the + 1.
                printf("%lu: ", coord[r][k] + 1);
            } else {
                printf("%lu,", coord[r][k] + 1);
            }
        }
        printf("%f\n", values[k]);
    }
}

void printCSRAsMTX(const std::vector<uint64_t> dims,
                   const std::vector<int> &col,
                   const std::vector<int> &rowptr,
                   const std::vector<double>  &csrValues) {
    uint64_t nr = dims[0];
    for (int i = 0; i< nr; i++) {
        for (int k=rowptr[i]; k<rowptr[i+1]; k++) {
            int j = col[k];
            printf("%d,%d: %f\n", i+1, j+1, csrValues[k]);
        }
    }

}

struct CSR {
    CSR(int nr, int nnz) {
        col = std::vector<int>(nnz);
        rowptr = std::vector<int>(nr);
        values = std::vector<double>(nnz);
    }

public:
    std::vector<int> col;
    std::vector<int> rowptr;
    std::vector<double> values;
};

CSR* COOToCSR(uint64_t nnz, uint64_t rank,
              const std::vector<uint64_t> &dims,
              std::vector<double> &cooValues,
              std::vector<std::vector<uint64_t>> &coord) {
    int nr = dims[0];
    int nc = dims[1];

    CSR *csr = new CSR(nr, nnz);

#define EX_ROW1(n) coord[0][n]
#define EX_COL1(n) coord[1][n]
#define EX_ACOO(n) cooValues[n]
#define EX_COL2(n) csr->col[n]
#define EX_ROWPTR(n) csr->rowptr[n]
#define EX_ACSR(n) csr->values[n]
#define NR nr
#define NC nc
#define NNZ nnz

#include <coo_csr.h>

#undef EX_ROW1
#undef EX_COL1
#undef EX_ACOO
#undef EX_COL2
#undef EX_ROWPTR
#undef EX_ACSR
#undef NR
#undef NC
#undef NNZ

    return csr;
}

int main(int argc, char * argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <filename of matrix market exchange file>\n", argv[0]);
        exit(1);
    }
    char * filename = argv[1];
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Cannot find %s\n", filename);
        exit(1);
    }
    char line[kColWidth];
    uint64_t idata[512];
    bool isSymmetric = false;
    readMMEHeader(file, filename, line, idata, &isSymmetric);
    uint64_t rank = idata[0];
    uint64_t nnz = idata[1];

    auto dims = std::vector<uint64_t>(rank);
    for(int i = 0; i < rank; i++) {
        dims[i] = idata[i + 2];
    }

    std::vector<std::vector<uint64_t>> coord(rank, std::vector<uint64_t>(nnz));
    std::vector<double> values(nnz);

    // Read file into vectors
    for (uint64_t k = 0; k < nnz; k++) {
        if (!fgets(line, kColWidth, file)) {
            fprintf(stderr, "Cannot find next line of data in %s\n", filename);
            exit(1);
        }
        char *linePtr = line;
        for (uint64_t r = 0; r < rank; r++) {
            uint64_t idx = strtoul(linePtr, &linePtr, 10);
            coord[r][k] = idx - 1;
        }

        double value = strtod(linePtr, &linePtr);
        values[k] = value;
    }

    printf("Before========\n");
    printCOOAsMTX(nnz, rank, values, coord);


    CSR *csr = COOToCSR(nnz, rank, dims, values, coord);

    printf("After=========\n");
    printCSRAsMTX(dims, csr->col, csr->rowptr, csr->values);

    return 0;
}

