#include <climits>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

int main()
{
    mkdir("fake_data", 0777);
    for (int n = 100; n <= 10000; n += 100)
    {
        std::stringstream ss;
        ss << "./fake_data/fake_nnz_" << n << ".mtx";
        const std::string tmp = ss.str();
        const char *filename = tmp.c_str();

        FILE *file = fopen(filename, "w");
        if (!file)
        {
            fprintf(stderr, "Cannot find %s\n", filename);
            exit(1);
        }

        fprintf(file, "%%%%MatrixMarket matrix coordinate real general\n");
        fprintf(file, "%d %d %d\n", n / 2, n / 2, n);
        for (int j = 0; j < n; j++)
        {
            int row = (j / 2) + 1;
            int col = row + (j % 2);
            if (col == n / 2 + 1)
            {
                col -= 2;
            }
            fprintf(file, "%d %d %f \n", row, col, row + (double)col / (double)n);
        }

        fclose(file);
    }
}
