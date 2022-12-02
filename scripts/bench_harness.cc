// Parts of this code are taken from
// https://github.com/llvm/llvm-project/blob/b682616d1fd1263b303985b9f930c1760033af2c/mlir/lib/ExecutionEngine/SparseTensorUtils.cpp
// Which is part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// These sections have been marked.
#include <climits>
#include <cctype>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>
#include <synth.h>
#include <Permute.h>
#include <chrono>
#include <algorithm>  // std::sort
#include <numeric>    // std::iota
#include <user_defs.h>

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

struct COO {
    explicit COO(uint64_t nnz, uint64_t rank) : nnz(nnz), rank(rank) {
        coord = std::vector<std::vector<uint64_t>>(rank, std::vector<uint64_t>(nnz));
        values = std::vector<double>(nnz);
    }

public:
    const uint64_t nnz;
    const uint64_t rank;
    std::vector<std::vector<uint64_t>> coord;
    std::vector<double> values;
};

bool equal(const COO &rhs, const COO &lhs) {
    if (rhs.rank != lhs.rank) {
        return false;
    }
    if (rhs.nnz != lhs.nnz) {
        return false;
    }
    for (uint64_t r = 0; r < rhs.rank; r++) {
        for (uint64_t i = 0; i < rhs.nnz; i++) {
            if (rhs.coord[r][i] != lhs.coord[r][i]) {
                return false;
            }
        }
    }
    for (uint64_t i = 0; i < rhs.nnz; i++) {
        if (rhs.values[i] != lhs.values[i]) {
            return false;
        }
    }

    return true;
}

/// Taken from LLVM SparseTnesorUitls (see top for details).
///
/// Read the MME header of a general sparse matrix of type real.
static void readMMEHeader(FILE *file, char *filename, char *line,
                          uint64_t *idata, bool *isPattern, bool *isSymmetric) {
    char header[64];
    char object[64];
    char format[64];
    char field[64];
    char symmetry[64];
    // Read header line.
    if (fscanf(file, "%63s %63s %63s %63s %63s\n", header, object, format, field, symmetry) != 5) {
        fprintf(stderr, "Corrupt header in %s\n", filename);
        exit(1);
    }
    // Set properties
    *isPattern = (strcmp(toLower(field), "pattern") == 0);
    *isSymmetric = (strcmp(toLower(symmetry), "symmetric") == 0);
    // Make sure this is a general sparse matrix.
    if (strcmp(toLower(header), "%%matrixmarket") ||
        strcmp(toLower(object), "matrix") ||
        strcmp(toLower(format), "coordinate") ||
        (strcmp(toLower(field), "real") && !(*isPattern)) ||
        (strcmp(toLower(symmetry), "general") && !(*isSymmetric))) {
        fprintf(stderr, "Cannot find a general sparse matrix in %s\n", filename);
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
    if (sscanf(line, "%" PRIu64 "%" PRIu64 "%" PRIu64 "\n", idata + 2, idata + 3, idata + 1) != 3) {
        fprintf(stderr, "Cannot find size in %s\n", filename);
        exit(1);
    }
}

/// Taken from LLVM SparseTnesorUitls (see top for details).
///
/// Read the "extended" FROSTT header. Although not part of the documented
/// format, we assume that the file starts with optional comments followed
/// by two lines that define the rank, the number of nonzeros, and the
/// dimensions sizes (one per rank) of the sparse tensor.
static void readExtFROSTTHeader(FILE *file, char *filename, char *line,
                                uint64_t *idata) {
    // Skip comments.
    while (true) {
        if (!fgets(line, kColWidth, file)) {
            fprintf(stderr, "Cannot find data in %s\n", filename);
            exit(1);
        }
        if (line[0] != '#')
            break;
    }
    // Next line contains RANK and NNZ.
    if (sscanf(line, "%" PRIu64 "%" PRIu64 "\n", idata, idata + 1) != 2) {
        fprintf(stderr, "Cannot find metadata in %s\n", filename);
        exit(1);
    }
    // Followed by a line with the dimension sizes (one per rank).
    for (uint64_t r = 0; r < idata[0]; r++)
        if (fscanf(file, "%" PRIu64, idata + 2 + r) != 1) {
            fprintf(stderr, "Cannot find dimension size %s\n", filename);
            exit(1);
        }
    fgets(line, kColWidth, file); // end of line
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
    explicit CSR(int nr, int nnz) : nr(nr), nnz(nnz) {
        col = std::vector<int>(nnz);
        rowptr = std::vector<int>(nr + 1);
        values = std::vector<double>(nnz);
    }

public:
    const int nr;
    const int nnz;
    std::vector<int> col;
    std::vector<int> rowptr;
    std::vector<double> values;
};


struct CSC {
    explicit CSC(int nc, int nnz) : nc(nc), nnz(nnz) {
        row = std::vector<int>(nnz);
        colptr = std::vector<int>(nc + 1);
        values = std::vector<double>(nnz);
    }

public:
    const int nc;
    const int nnz;
    std::vector<int> row;
    std::vector<int> colptr;
    std::vector<double> values;
};

struct DIA {
    DIA(int nd, int nr, int nc) {
        off = std::vector<int>(nd);
        values = std::vector<double>((nr < nc ? nr : nc) * nd);
    }

    DIA() {
    }

public:
    std::vector<int> off;
    std::vector<double> values;
};

std::pair<CSC *, double> CSRToCSC(uint64_t nnz, uint64_t rank,
                                  const std::vector<uint64_t> &dims,
                                  const CSR &csr) {

    int nr = dims[0];
    int nc = dims[1];

    CSC *csc = new CSC(nr, nnz);
    std::vector<int> &row = csc->row;
    std::vector<int> &colptr = csc->colptr;
    std::vector<double> &cscValues = csc->values;
    const std::vector<int> &col = csr.col;
    const std::vector<int> &rowptr = csr.rowptr;
    const std::vector<double> &csrValues = csr.values;

#define EX_COL(n) col[n]
#define EX_ROW(n) row[n]
#define EX_ACSR(n) csrValues[n]
#define EX_ROWPTR(n) rowptr[n]
#define EX_COLPTR(n) colptr[n]
#define EX_ACSC(n) cscValues[n]
#define NR nr
#define NC nc
#define NNZ nnz

#include <csr_csc_opt2_auto.h>

#undef EX_ROW
#undef EX_COL
#undef EX_ACSR
#undef EX_COLPTR
#undef EX_ROWPTR
#undef EX_ACSC
#undef NR
#undef NC
#undef NNZ

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = stop - start;

    return {csc, fp_ms.count()};
}

std::pair<COO *, uint64_t> COOToMCOO3D(uint64_t nnz, uint64_t rank,
                                       const std::vector<uint64_t> &dims,
                                       const COO &coo) {

    int nr = dims[0];
    int nc = dims[1];
    int nz = dims[2];
    COO *mcoo = new COO(nnz, rank);
    std::vector<std::vector<uint64_t>> &mCoord = mcoo->coord;
    std::vector<double> &mCooValues = mcoo->values;
    const std::vector<std::vector<uint64_t>> &coord = coo.coord;
    const std::vector<double> &cooValues = coo.values;
#define EX_ROW3D(n) coord[0][n]
#define EX_COL3D(n) coord[1][n]
#define EX_Z3D(n) coord[2][n]
#define EX_A3DCOO(n) cooValues[n]
#define EX_MROW3(n) mCoord[0][n]
#define EX_MCOL3(n) mCoord[1][n]
#define EX_MZ3D(n) mCoord[2][n]
#define EX_A3DMCOO(n) mCooValues[n]
#define NR nr
#define NC nc
#define NZ nz
#define NNZ nnz

#include <coo3d_mcoo3d_opt2.h>

#undef EX_ROW3D
#undef EX_COL3D
#undef EX_Z3D
#undef EX_A3DCOO
#undef EX_MROW3
#undef EX_MCOL3
#undef EX_MZ3D
#undef EX_A3DMCOO
#undef NR
#undef NC
#undef NZ
#undef NNZ

    std::chrono::duration<double, std::milli> fp_ms = stop - start;

    return {mcoo, fp_ms.count()};
}

std::pair<COO *, uint64_t> COOToMCOO(uint64_t nnz, uint64_t rank,
                                     const std::vector<uint64_t> &dims,
                                     const COO &coo) {

    int nr = dims[0];
    int nc = dims[1];
    auto start = std::chrono::high_resolution_clock::now();
    COO *mcoo = new COO(nnz, rank);
    std::vector<std::vector<uint64_t>> &mCoord = mcoo->coord;
    std::vector<double> &mCooValues = mcoo->values;
    const std::vector<std::vector<uint64_t>> &coord = coo.coord;
    const std::vector<double> &cooValues = coo.values;
#define EX_ROW1(n) coord[0][n]
#define EX_COL1(n) coord[1][n]
#define EX_ACOO(n) cooValues[n]
#define EX_ROW3(n) mCoord[0][n]
#define EX_COL3(n) mCoord[1][n]
#define EX_AMCOO(n) mCooValues[n]
#define NR nr
#define NC nc
#define NNZ nnz

#include <coo_mcoo_opt.h>

#undef EX_ROW1
#undef EX_COL1
#undef EX_ACOO
#undef EX_ROW3
#undef EX_COL3
#undef EX_AMCOO
#undef NR
#undef NC
#undef NNZ

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = stop - start;

    return {mcoo, fp_ms.count()};
}

COO *createSorted(const COO &coo) {
    COO *sorted = new COO(coo.nnz, coo.rank);

    // initialize original index locations
    std::vector<size_t> idx(coo.values.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based row
    std::sort(idx.begin(), idx.end(),
              [&coo](size_t i1, size_t i2) {
                  // lexicographic sort based on coordinate values
                  for (uint64_t i = 0; i < coo.rank; i++) {
                      auto one = coo.coord[i][i1];
                      auto two = coo.coord[i][i2];
                      if (one != two) {
                          return one < two;
                      }
                  }
                  // There shouldn't be duplicates, but if there are make sure they end up in a deterministic order.
                  if (coo.values[i1] != coo.values[i2]) {
                      return coo.values[i1] < coo.values[i2];
                  }
                  return false;
              });

    // load data into new COO matrix based on row index
    for (uint64_t k = 0; k < coo.nnz; k++) {
        sorted->values[k] = coo.values[idx[k]];
        for (uint64_t i = 0; i < coo.rank; i++) {
            sorted->coord[i][k] = coo.coord[i][idx[k]];
        }
    }

    return sorted;
}

std::pair<COO *, double> COOToSortedCOO(uint64_t nnz, uint64_t rank,
                                        const COO &coo) {
    auto start = std::chrono::high_resolution_clock::now();
    COO *sorted = createSorted(coo);
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = stop - start;


    return {sorted, fp_ms.count()};
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

std::pair<CSC *, double> COOToCSC(uint64_t nnz, uint64_t rank,
                                  const std::vector<uint64_t> &dims,
                                  const COO &coo) {
    int nr = dims[0];
    int nc = dims[1];


    CSC *csc = new CSC(nr, nnz);
    std::vector<int> &row = csc->row;
    std::vector<int> &colptr = csc->colptr;
    std::vector<double> &values = csc->values;
    const std::vector<std::vector<uint64_t>> &coord = coo.coord;
    const std::vector<double> &cooValues = coo.values;

#define EX_ROW1(n) coord[0][n]
#define EX_COL1(n) coord[1][n]
#define EX_ACOO(n) cooValues[n]
#define EX_ROW(n) row[n]
#define EX_COLPTR(n) colptr[n]
#define EX_ACSC(n) values[n]
#define NR nr
#define NC nc
#define NNZ nnz

#include <coo_csc_opt2.h>

#undef EX_ROW1
#undef EX_COL1
#undef EX_ACOO
#undef EX_COL2
#undef EX_COLPTR
#undef EX_ACSC
#undef NR
#undef NC
#undef NNZ

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = stop - start;

    return {csc, fp_ms.count()};
}


std::pair<DIA *, double> COOToDIA(uint64_t nnz, uint64_t rank,
                                  const std::vector<uint64_t> &dims,
                                  const COO &coo) {
    uint64_t nr = dims[0];
    uint64_t nc = dims[1];

    DIA *dia = new DIA();
    std::vector<int> &offset = dia->off;
    std::vector<double> &values = dia->values;
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t nd = 0;
    const std::vector<std::vector<uint64_t>> &coord = coo.coord;
    const std::vector<double> &cooValues = coo.values;

#define EX_ROW1(n)  coord[0][n]
#define EX_COL1(n)  coord[1][n]
#define EX_ACOO(n)  cooValues[n]
#define EX_ADIA(kd) values[kd]
#define NR nr
#define NC nc
#define NNZ nnz
#define ND nd

#include <coo_dia2.h>

#undef EX_ROW1
#undef EX_COL1
#undef EX_ACOO
#undef EX_ADIA
    // Copy out the values from off array
    auto stop = std::chrono::high_resolution_clock::now();
    for (uint64_t h = 0; h < off->getSize(); h++) {
        offset.push_back(off->get(h));
    }


    std::chrono::duration<double, std::milli> fp_ms = stop - start;

    return {dia, fp_ms.count()};

}

void output(const COO &beforeConversion,
            const COO &afterConversion, double milliseconds, char *conversion, char *filename) {
    auto before = createSorted(beforeConversion);
    auto after = createSorted(afterConversion);
    if (equal(*before, *after)) {
        printf("[PASS] ");
    } else {
        printf("[FAIL] ");
    }
    printf("%s %s, time: %f milliseconds, \n", conversion, filename, milliseconds);
    delete (after);
    delete (before);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <filename (should be in matrix market exchange format)> "
                        "<type of conversion to run options: coo_csr,sort,coo_mcoo,coo_mcoo3d,csr_csc,coo_csc,coo_dia>\n",
                argv[0]);
        exit(1);
    }
    char *filename = argv[1];
    char *conversion = argv[2];

    // Read data out of file =====
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Cannot find %s\n", filename);
        exit(1);
    }
    char line[kColWidth];
    uint64_t idata[512];
    bool isSymmetric = false;
    bool isPattern = false;
    if (strstr(filename, ".mtx")) {
        readMMEHeader(file, filename, line, idata, &isPattern, &isSymmetric);
    } else if (strstr(filename, ".tns")) {
        readExtFROSTTHeader(file, filename, line, idata);
    } else {
        fprintf(stderr, "Unknown format %s\n", filename);
        exit(1);
    }
    uint64_t rank = idata[0];
    uint64_t nnz = idata[1];

    auto dims = std::vector<uint64_t>(rank);
    for (uint64_t i = 0; i < rank; i++) {
        dims[i] = idata[i + 2];
    }

    COO beforeConversion = COO(nnz, rank);


    std::vector<std::vector<uint64_t>> &coord = beforeConversion.coord;
    std::vector<double> &values = beforeConversion.values;

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

    // create validation fuction ===== auto

    // The number of runs to average over
    int n = 25;

    if (strcmp(conversion, "coo_csr") == 0) {
        auto p1 = COOToSortedCOO(nnz, rank, beforeConversion);
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


        COO afterConversion = COO(nnz, rank);
        for (uint64_t i = 0; i < csr->nr; i++) {
            for (uint64_t k = csr->rowptr[i]; k < csr->rowptr[i + 1]; k++) {
                int j = csr->col[k];
                afterConversion.coord[0][k] = i;
                afterConversion.coord[1][k] = j;
                afterConversion.values[k] = csr->values[k];
            }
        }

        output(beforeConversion, afterConversion, milliseconds, conversion, filename);
        delete (csr);
    } else if (strcmp(conversion, "sort") == 0) {
        COO *sortedCoo;
        double milliseconds = 0;
        for (int i = 0; i < n; i++) {
            auto p = COOToSortedCOO(nnz, rank, beforeConversion);
            sortedCoo = p.first;
            milliseconds += p.second;
            if (i != n - 1) {
                delete (sortedCoo);
            }
        }
        milliseconds /= n;

        output(beforeConversion, *sortedCoo, milliseconds, conversion, filename);
        delete (sortedCoo);
    } else if (strcmp(conversion, "coo_mcoo") == 0) {
        COO *mcoo;
        double milliseconds = 0;
        for (int i = 0; i < n; i++) {
            auto p = COOToMCOO(nnz, rank, dims, beforeConversion);
            mcoo = p.first;
            milliseconds += p.second;
            if (i != n - 1) {
                delete (mcoo);
            }
        }
        milliseconds /= n;

        output(beforeConversion, *mcoo, milliseconds, conversion, filename);
        delete (mcoo);
    } else if (strcmp(conversion, "coo_mcoo3d") == 0) {
        assert(dims.size() == 3 && "Dims must be 3D");
	    COO *mcoo;
        double milliseconds = 0;
        for (int i = 0; i < n; i++) {
            auto p = COOToMCOO3D(nnz, rank, dims, beforeConversion);
            mcoo = p.first;
            milliseconds += p.second;
            if (i != n - 1) {
                delete (mcoo);
            }
        }
        milliseconds /= n;

        output(beforeConversion, *mcoo, milliseconds, conversion, filename);
        delete (mcoo);
    } else if (strcmp(conversion, "csr_csc") == 0) {
        auto p1 = COOToSortedCOO(nnz, rank, beforeConversion);
        auto sortedCoo = p1.first;
        // COO to CSR first
        CSR *csr;
        auto p = COOToCSR(nnz, rank, dims, *sortedCoo);
        csr = p.first;

        CSC *csc;
        double milliseconds = 0;
        for (int i = 0; i < n; i++) {
            auto p = CSRToCSC(nnz, rank, dims, *csr);
            csc = p.first;
            milliseconds += p.second;
            if (i != n - 1) {
                delete (csc);
            }
        }
        milliseconds /= n;

        COO afterConversion = COO(nnz, rank);
        for (uint64_t j = 0; j < csc->nc; j++) {
            for (uint64_t k = csc->colptr[j]; k < csc->colptr[j + 1]; k++) {
                int i = csc->row[k];
                afterConversion.coord[0][k] = i;
                afterConversion.coord[1][k] = j;
                afterConversion.values[k] = csc->values[k];
            }
        }
        output(beforeConversion, afterConversion, milliseconds, conversion, filename);
        delete (csc);
        delete (csr);
    } else if (strcmp(conversion, "coo_csc") == 0) {
        auto p1 = COOToSortedCOO(nnz, rank, beforeConversion);
        auto sortedCoo = p1.first;


        CSC *csc;
        double milliseconds = 0;
        for (int i = 0; i < n; i++) {
            auto p = COOToCSC(nnz, rank, dims, *sortedCoo);
            csc = p.first;
            milliseconds += p.second;
            if (i != n - 1) {
                delete (csc);
            }
        }
        milliseconds /= n;

        COO afterConversion = COO(nnz, rank);
        for (uint64_t j = 0; j < csc->nc; j++) {
            for (uint64_t k = csc->colptr[j]; k < csc->colptr[j + 1]; k++) {
                int i = csc->row[k];
                afterConversion.coord[0][k] = i;
                afterConversion.coord[1][k] = j;
                afterConversion.values[k] = csc->values[k];
            }
        }
        output(beforeConversion, afterConversion, milliseconds, conversion, filename);
        delete (csc);
        delete (sortedCoo);
    } else if (strcmp(conversion, "coo_dia") == 0) {
        DIA *dia;
        double milliseconds = 0;
        for (int i = 0; i < n; i++) {
            auto p = COOToDIA(nnz, rank, dims, beforeConversion);
            dia = p.first;
            milliseconds += p.second;
            if (i != n - 1) {
                delete (dia);
            }
        }
        milliseconds /= n;


        COO afterConversion = COO(nnz, rank);
        uint64_t nr = dims[0];
	int n = 0;
        for (uint64_t i = 0; i < nr; i++) {
            for (uint64_t d = 0; d < dia->off.size(); d++) {
                int j = i + dia->off[d];
                long k = dia->off.size() * i + d;
                if (dia->values[k]!=INT_MIN){ 
		afterConversion.coord[0][n] = i;
                afterConversion.coord[1][n] = j;
                afterConversion.values[n] = dia->values[k];
                n++;
	       	}
            }
	}
	// order wont stay the same after conversion so it is necessary to 
	// sort after getting back coo.
        output(beforeConversion, afterConversion, milliseconds, conversion, filename);
        delete (dia);
    } else {
        printf("unknown conversion: %s\n", conversion);
        exit(1);
    }
}
