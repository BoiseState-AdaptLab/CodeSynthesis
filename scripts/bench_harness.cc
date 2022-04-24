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
#include <opt_synth.h>
#include <chrono>
#include <algorithm>  // std::sort
#include <numeric>    // std::iota

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

// TODO: if necessary make this work on arbitrary dimensions
template<typename Fn>
bool verify(const std::vector<std::vector<uint64_t>> &coords,
            const std::vector<double> &values,
            Fn f) {
    uint64_t rank = coords.size();
    uint64_t nnz = coords[0].size();
    for (int k = 0; k < nnz; k++) {
        std::vector<uint64_t> c(rank);
        for (int r = 0; r < rank; r++) {
            c[r] = coords[r][k];
        }
        double out = f(c);
        if (out != values[k]) {
            return false;
        }
    }
    return true;
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
            if (r == rank - 1) {
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
                   const std::vector<double> &csrValues) {
    uint64_t nr = dims[0];
    for (int i = 0; i < nr; i++) {
        for (int k = rowptr[i]; k < rowptr[i + 1]; k++) {
            int j = col[k];
            printf("%d,%d: %f\n", i + 1, j + 1, csrValues[k]);
        }
    }
}

struct CSR {
    CSR(int nr, int nnz) {
        col = std::vector<int>(nnz);
        rowptr = std::vector<int>(nr + 1);
        values = std::vector<double>(nnz);
    }

public:
    std::vector<int> col;
    std::vector<int> rowptr;
    std::vector<double> values;
};

struct COO {
    COO(uint64_t nnz, uint64_t rank) {
        coord = std::vector<std::vector<uint64_t>>(rank, std::vector<uint64_t>(nnz));
        values = std::vector<double>(nnz);
    }

public:
    std::vector<std::vector<uint64_t>> coord;
    std::vector<double> values;
};


std::pair<COO *, uint64_t> COOToSortedCOO(uint64_t nnz, uint64_t rank,
                                          const COO &coo) {
    auto start = std::chrono::high_resolution_clock::now();
    COO *sorted = new COO(nnz, rank);
    sorted->coord = coo.coord;
    sorted->values = coo.values;

    // initialize original index locations
    std::vector<size_t> idx(coo.values.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based row
    std::sort(idx.begin(), idx.end(),
              [&coo](size_t i1, size_t i2) { return coo.coord[0][i1] < coo.coord[0][i2]; });

    // load data into new COO matrix based on row index
    for (int k = 0; k < nnz; k++) {
        sorted->values[k] = coo.values[idx[k]];
        sorted->coord[0][k] = coo.coord[0][idx[k]];
        sorted->coord[1][k] = coo.coord[1][idx[k]];
    }

    auto stop = std::chrono::high_resolution_clock::now();
    uint64_t milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    return {sorted, milliseconds};
}

std::pair<CSR *, double> COOToCSR(uint64_t nnz, uint64_t rank,
                                  const std::vector<uint64_t> &dims,
                                  const COO &coo) {
    int nr = dims[0];
    int nc = dims[1];

    auto start = std::chrono::high_resolution_clock::now();

    CSR *csr = new CSR(nr, nnz);
    std::vector<int> &col = csr->col;
    std::vector<int> &rowptr = csr->rowptr;
    std::vector<double> &values = csr->values;
    const std::vector<std::vector<uint64_t>> &coord = coo.coord;
    const std::vector<double> &cooValues = coo.values;

#define EX_ROW1(n) coord[0][n]
#define EX_COL1(n) coord[1][n]
#define EX_ACOO(n) cooValues[n]
#define EX_COL2(n) col[n]
#define EX_ROWPTR(n) rowptr[n]
#define EX_ACSR(n) values[n]
#define NR nr
#define NC nc
#define NNZ nnz

#include <scoo_csr_opt.h>

#undef EX_ROW1
#undef EX_COL1
#undef EX_ACOO
#undef EX_COL2
#undef EX_ROWPTR
#undef EX_ACSR
#undef NR
#undef NC
#undef NNZ

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = stop - start;

    return {csr, fp_ms.count()};
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <filename (should be in matrix market exchange format)> "
                        "<type of conversion to run options: csr, sort> <validate: t/f>\n", argv[0]);
        exit(1);
    }
    char *filename = argv[1];
    char *conversion = argv[2];
    bool validate = false;
    if (strcmp(argv[3], "t") == 0) {
        validate = true;
    }

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
    for (int i = 0; i < rank; i++) {
        dims[i] = idata[i + 2];
    }

    COO coo = COO(nnz, rank);

    std::vector<std::vector<uint64_t>> &coord = coo.coord;
    std::vector<double> &values = coo.values;

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

    auto output = [validate, &coord, &values, conversion, filename](
            std::function<double(const std::vector<uint64_t> &)> &&check, double milliseconds) {
        if (validate) {
            if (verify(coord, values, check)) {
                printf("[PASS] coo->%s %s, time: %f milliseconds, \n", conversion, filename, milliseconds);
            } else {
                printf("[FAIL] coo->%s %s, time: %f milliseconds, \n", conversion, filename, milliseconds);
            }
        } else {
            printf("%s, %f, \n", filename, milliseconds);
        }
    };

    int n = validate ? 1 : 25;

    if (strcmp(conversion, "csr") == 0) {
        auto p1 = COOToSortedCOO(nnz, rank, coo);
        auto sortedCOO = p1.first;

        CSR *csr;
        double milliseconds = 0;
        for (int i = 0; i < n; i++) {
            auto p = COOToCSR(nnz, rank, dims, *sortedCOO);
            csr = p.first;
            milliseconds += p.second;
            if (i != n - 1) {
                delete (csr);
            }
        }
        milliseconds /= n;
        delete (sortedCOO);

        output([&csr, dims](const std::vector<uint64_t> &cord) -> double {
            uint64_t inI = cord[0];
            uint64_t inJ = cord[1];
            uint64_t nr = dims[0];
            for (uint64_t i = 0; i < nr; i++) {
                for (uint64_t k = csr->rowptr[i]; k < csr->rowptr[i + 1]; k++) {
                    int j = csr->col[k];
                    if (i == inI && j == inJ) {
                        return csr->values[k];
                    }
                }
            }
            return 0;
        }, milliseconds);
        delete (csr);
    } else if (strcmp(conversion, "sort") == 0) {
        auto p = COOToSortedCOO(nnz, rank, coo);
        COO *sortedCOO = p.first;
        double milliseconds = p.second;

        output([sortedCOO, nnz](const std::vector<uint64_t> &cord) -> double {
            uint64_t inI = cord[0];
            uint64_t inJ = cord[1];
            for (int k = 0; k < nnz; k++) {
                if (sortedCOO->coord[0][k] == inI && sortedCOO->coord[1][k] == inJ) {
                    return sortedCOO->values[k];
                }
            }
            return 0;
        }, milliseconds);
        delete (sortedCOO);
    } else {
        printf("unknown conversion: %s\n", conversion);
        exit(1);
    }
}
