#include <iostream>
#include <string>
#include <iegenlib/computation/Computation.h>
#include <CodeSynthesis.h>
#include <map>
#include <iegenlib/set_relation/environment.h>
#include <iegenlib/set_relation/set_relation.h>
#include <iostream>
#include <fstream>

using namespace code_synthesis;
int main(int argc, char**argv) {
    if(argc != 5 && argc != 9) {
        std::cerr << "Usage: synthDriver -src "
                  <<"<formatname>,<dataName> -dest <formatName>,<dataName> "
                  <<" [-fuse <fuse-list> -fuselevel <level>]\n";
        return 0;
    }
    // Load preset optimizations
    std::vector<int> fuseStmts;
    int fuseLevel = 0;
    //Load preset formats
    std::map<std::string,SparseFormat*> supportedFormats;
    SparseFormat * coo = new SparseFormat();
    coo->mapToDense = "{[n,ii ,jj] -> [i,j]:"
                      " row1(n) = ii and 0 <= n and n < NNZ "
                      "and col1(n) = jj and i = ii and"
		      " j = jj and i >= 0 and "
                      " i < NR and j >= 0 and j < NC}";
    coo->knowns = { "NR","NC","NNZ"};
    coo->ufQuants = { UFQuant( "{[x]:0 <= x < NNZ}","{[i]: 0 <= i <= NC}",
                               "col1",false, Monotonic_NONE),
                      UFQuant( "{[x]:0 <= x < NNZ}","{[i]: 0 <= i <= NR}",
                               "row1",false, Monotonic_NONE)
                    };
    coo->dataConstraint = "A[n]!=0";
    coo->dataAccess = "{[n,ii,jj] -> [n]}";
    supportedFormats["COO"] = coo;

    SparseFormat * bcsr = new SparseFormat();
    bcsr->mapToDense = "{[ii,kk,jj,hr,hc, is, js]->[i,j] : 0<= ii < NR_BR &&"
                       " browptr(ii) <= kk < browptr(ii+ 1) && jj= bcol(kk) && 0 <= hr"
                       " < BR && 0 <= hc < BC && ii = int_div(i ,BR)"
		       " && jj = int_div(j , BC)  && i = ii * 99 + hr &&"
                       " j = jj * 999 + hc && i = is && j = js && 0 <= i < NR"
		       " && 0 <= j < NC }";
    bcsr->dataAccess = "{[hr,hc,ii,jj,kk, is , js] -> [p]: p= kk * 9999 + hr * 999 + hc }";

    bcsr->knowns = { "BR","BC","NR", "NC", "NC_BC","NR_BR","BRBC","int_div"};
    bcsr->ufQuants = { UFQuant( "{[ii]: 0 <= ii <= NR_BR}",
                                "{[ii]:0 <= ii < NR_BR}",
                                "browptr",false, Monotonic_Nondecreasing),
                       UFQuant( "{[jj]:0 <= jj < NC_BC}",
                                "{[jj]: 0 <= jj <= NC_BC}",
                                "bcol",false, Monotonic_NONE)
                     };
    bcsr->dataConstraint = "A(i,j) != 0 || ii,jj, x | A(i,j) = 0 && "
                           "BR ∗ x <= i < BR ∗ ( x + 1) && BC ∗ x <= j < BC ∗ (x+ 1) && "
                           "BR∗x <= ii < BR∗ ( x + 1) & BC ∗ x <= jj< BC ∗ (x+ 1) && "
                           "A(ii,jj) != 0";
    supportedFormats["BCSR"] = bcsr;


    SparseFormat * csr = new SparseFormat();
    csr->mapToDense = "{[ii,k, jj]->[i,j]: jj = j and ii= i and i >= 0 and i < NR and"
                      " j >= 0 and j < NC and rowptr(ii) <= k < rowptr(ii+1)"
                      " and jj = col2(k)}";
    csr->dataAccess = "{[ii,k,jj] -> [k]}";
    csr->knowns = { "NR","NC","NNZ"};
    csr->dataConstraint = "A[k]!=0";
    csr->ufQuants = { UFQuant( "{[i]: 0 <= i <= NR}","{[x]:0 <= x < NNZ}",
                               "rowptr",false, Monotonic_Nondecreasing),
                      UFQuant( "{[x]:0 <= x < NNZ}","{[i]: 0 <= i <= NC}",
                               "col2",false, Monotonic_NONE,
                               "{[e1,e2]: e1 < e2}",
                               "{[e1,e2]: i * NR + e1 < j * NC + e2 ")
                    };
    supportedFormats["CSR"] = csr;



    SparseFormat * csc = new SparseFormat();
    csc->mapToDense = "{[jj,k,ii]->[i,j]: j = jj and  i = ii and  i >= 0 and i < NR and"
                      " j >= 0 and j < NC and colptr(jj) <= k < colptr(jj+1)"
                      " and ii = row4(k) and 0 <= jj < NC}";
    csc->dataAccess = "{[jj,k,ii] -> [k]}";
    csc->knowns = { "NR","NC","NNZ"};
    csc->dataConstraint = "A[k]!=0";
    csc->ufQuants = { UFQuant( "{[i]: 0 <= i <= NR}","{[x]:0 <= x < NNZ}",
                               "colptr",false, Monotonic_Nondecreasing),
                      UFQuant( "{[x]:0 <= x < NNZ}","{[i]: 0 <= i <= NR}",
                               "row4",false, Monotonic_NONE,
                               "{[e1,e2]: e1 < e2}",
                               "{[e1,e2]: j * NC + e1 < j * NC + e2 ")
                    };
    supportedFormats["CSC"] = csc;


    SparseFormat * dia = new SparseFormat();
    dia->mapToDense = "{[id,dd,jj] -> [ii,j]: jj = j && id = ii && 0 <= id"
                      " < NR &&  0 <= id < NR && 0<= ii < NR &&"
                      " 0 <= dd < ND "
                      " && j = ii + off(dd) && 0 <= j < NC"
                      "}";
    // 99 is ND, direct replacement is required.
    dia->dataAccess = "{[id,dd, jj] -> [kd]: kd = 99 * id + dd}";
    dia->knowns = { "NR","NC","NNZ"};
    dia->ufQuants = { UFQuant( "{[i]: 0 <= i < ND}","{[x]:0 <= x < NNZ}",
                               "off",true, Monotonic_Increasing)};
    // TODO: represent that there exists some non zero for
    // every zero that exists in the same plane as the diagonal.
    dia->dataConstraint = "∃i′, j′, x|"
                          "A(i′, j′) = 0 ^ i' - j'  = i - j";
    // It is difficult to represent the rage of offset! It is the
    // count of the number of non zeros in each diagonal. So for now skip
    supportedFormats["DIA"] = dia;

    SparseFormat * mortonCoo = new SparseFormat();
    mortonCoo->mapToDense = "{[n1, i , j] -> [i,j]:"
                            " row3(n1) = i and 0 <= n1 and n1 < NNZ "
                            "and col3(n1) = j and  i >= 0 and "
                            " i < NR and j >= 0 and j < NC}";
    mortonCoo->knowns = { "NR","NC","NNZ"};
    mortonCoo->ufQuants = { UFQuant( "{[x]:0 <= x < NNZ}","{[i]: 0 <= i <= NC}",
                                     "col3",false, Monotonic_NONE, "{[e1,e2]: e1 < e2}",
                                     "{[e1,e2]: MORTON(row3(e1),col3(e1)) < "
                                     "MORTON(row3(e2),col3(e2))} "),

                            UFQuant( "{[x]:0 <= x < NNZ}","{[i]: 0 <= i <= NR}",
                                     "row3",false, Monotonic_NONE,"{[e1,e2]: e1 < e2}",
                                     "{[e1,e2]: MORTON(row3(e1),col3(e1)) < "
                                     "MORTON(row3(e2),col3(e2))} ")
                          };

    mortonCoo->dataConstraint = "A[n]!=0";
    mortonCoo->dataAccess = "{[n1, i, j] -> [n1]}";
    supportedFormats["MCOO"] = mortonCoo;
    // Parse command line
    int currIndex = 1;
    SparseFormat* sourceFormat = NULL;
    SparseFormat* destFormat = NULL;


    while(currIndex < argc ) {
        std::string argString (argv[currIndex]);
        if(argString == "-src") {
            std::string srcFormat(argv[++currIndex]);
            std::string sourceName;
            std::string dataName;
            int del = srcFormat.find(",");
            sourceName = srcFormat.substr(0,del);
            dataName = srcFormat.substr(del +1,srcFormat.size());
            sourceFormat = supportedFormats[sourceName];
            assert(sourceFormat && "Source format is not supported");
            sourceFormat->dataName = dataName;
        }
        if(argString == "-dest") {
            std::string srcFormat(argv[++currIndex]);
            std::string destName;
            std::string dataName;
            int del = srcFormat.find(",");
            destName = srcFormat.substr(0,del);
            dataName = srcFormat.substr(del +1,srcFormat.size());

            destFormat = supportedFormats[destName];
            assert(destFormat && "Destination format is not supported");
            destFormat->dataName = dataName;
        }
        if (argString == "-fuse") {
            std::string s(argv[++currIndex]);
            size_t pos = 0;
            std::string token;
            std::string delimiter = ",";
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                fuseStmts.push_back(std::stoi(token));
                s.erase(0, pos + delimiter.length());
            }
        }

        if (argString == "-fuselevel") {
            std::string level(argv[++currIndex]);
            fuseLevel = std::stoi(level);
        }
        currIndex++;
    }
    assert(sourceFormat && destFormat
           && "Unsopported Source or Destination Format");
    CodeSynthesis* synth = new CodeSynthesis(sourceFormat, destFormat);
    std::string code = synth->generateFullCode(fuseStmts,fuseLevel);
    std::ofstream fileOut;
    fileOut.open("synth.h");
    fileOut << synth->GetSupportHeader();
    fileOut.close();
    std::cout << "generated synth.h...\n";
    fileOut.open("synth.c");
    fileOut << "#include <synth.h> \n";
    fileOut << code;
    fileOut.close();
    std::cout << "generated synth.c ..\n";

    for (auto s : supportedFormats) {
        delete s.second;
    }
    return 0;
}
