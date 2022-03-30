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

// 1 based indexing to match mmx format, this function is intended for testing purposes.
void print_coo_1_based_index(uint64_t nnz, uint64_t rank, double * values, uint64_t **coord) {
    for (uint64_t k = 0; k < nnz; k++) {
        for (uint64_t r = 0; r < rank; r++) {
            if (r == rank-1) {
                printf("%lu: ", coord[r][k] + 1);
            } else {
                printf("%lu,", coord[r][k] + 1);
            }
        }
        printf("%f\n", values[k]);
    }
}

int main(int argc, char * argv[]) {
    // TODO: replace with file read from command line
    char * filename = (char *)"/home/aaron/Dev/c++/CodeSynthesis/test/data/test.mtx";
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Cannot find %s\n", filename);
        exit(1);
    }
    // Perform some file format dependent set up.
    char line[kColWidth];
    uint64_t idata[512];
    bool isSymmetric = false;
    readMMEHeader(file, filename, line, idata, &isSymmetric);
    uint64_t rank = idata[0];
    uint64_t nnz = idata[1];

    // Whether the given col/row has a value
    auto *dims = (uint64_t *) calloc(rank, sizeof(uint64_t));
    for(int i = 0; i < rank; i++) {
        dims[i] = idata[i + 2];
    }

    auto **coord = (uint64_t **) calloc(rank, sizeof(uint64_t *));
    for (int i = 0; i < rank; i++) {
        coord[i] = (uint64_t *) calloc(nnz, sizeof(uint64_t));
    }
    auto *values = (double *)calloc(nnz, sizeof(double));

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

    print_coo_1_based_index(nnz, rank, values, coord);

    return 0;
}

