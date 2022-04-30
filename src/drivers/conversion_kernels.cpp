#include <iostream>
#include <string>
#include <iegenlib/computation/Computation.h>
#include <CodeSynthesis.h>
#include <map>
#include <iegenlib/set_relation/environment.h>
#include <iegenlib/set_relation/set_relation.h>
#include <iostream>
#include <fstream>

using namespace iegenlib;
int main(){
   Computation * comp = new Computation("COO_MCOO");
   comp->addStmt(new Stmt("P0->insert({row1(n), col1(n)})",
			  "{ [n, tv1, tv2] : tv1 - row1(n) = 0 &&"
			  " tv2 - col1(n) = 0 && n >= 0 && col1(n) >= 0"
			  " && row1(n) >= 0 && NC - 1 >= 0 && NNZ - 1 >= 0"
			  " && NR - 1 >= 0 && -n + NNZ - 1 >= 0 && NC - col1(n)"
			  " - 1 >= 0 && NR - row1(n) - 1 >= 0 }",
			  "{ [tv0, tv1, tv2] -> [0, a1, 0, a3, 0, a5, 0] :"
			  " tv0 - a1 = 0 && tv1 - a3 = 0 && tv2 - a5 = 0 }",
{
{"row1","{ [tv0, tv1, tv2] -> [tv3] : tv0 - tv3 = 0 }"},
{    "col1", "{ [tv0, tv1, tv2] -> [tv3] : tv0 - tv3 = 0 }"}
},
{
{"P0", "{ [tv0, tv1, tv2] -> [0] }"}
}));
   comp->addStmt( new Stmt("col3(n)=P0_INV(n)[1]",
 "{ [n] : n >= 0 && P0(row1(n), col1(n)) >= 0 && col1(n) >= 0 && row1(n) >= 0 && NC - 1 >= 0 && NNZ - 1 >= 0 && NR - 1 >= 0 && -n + NNZ - 1 >= 0 && NC - col1(n) - 1 >= 0 && NNZ - P0(row1(n), col1(n)) - 1 >= 0 && NR - row1(n) - 1 >= 0 }",
 "{ [tv0] -> [1, a1, 0, 0, 0, 0, 0] : tv0 - a1 = 0 }",
{
{"col1", "{ [tv0] -> [tv1] : tv0 - tv1 = 0 }"},
{"P0", "{ [tv0] -> [tv1, tv2] : tv1 - row1(tv0) = 0 && tv2 - col1(tv0) = 0 }"},
{"row1", "{ [tv0] -> [tv1] : tv0 - tv1 = 0 }"},
{"col1", "{ [tv0] -> [tv1] : tv0 - tv1 = 0 }"}
},
{
{"col3", "{ [tv0] -> [tv1] : tv1 - P0(row1(tv0), col1(tv0)) = 0 }"}
}));



comp->addStmt( new Stmt("row3(n)=P0_INV(n)[0]",
 "{ [n] : n >= 0 && P0(row1(n), col1(n)) >= 0 && col1(n) >= 0 && row1(n) >= 0 && NC - 1 >= 0 && NNZ - 1 >= 0 && NR - 1 >= 0 && -n + NNZ - 1 >= 0 && NC - col1(n) - 1 >= 0 && NNZ - P0(row1(n), col1(n)) - 1 >= 0 && NR - row1(n) - 1 >= 0 }",
 "{ [tv0] -> [2, a1, 0, 0, 0, 0, 0] : tv0 - a1 = 0 }",
{
{"row1", "{ [tv0] -> [tv1] : tv0 - tv1 = 0 }"},
{"P0", "{ [tv0] -> [tv1, tv2] : tv1 - row1(tv0) = 0 && tv2 - col1(tv0) = 0 }"},
{"row1", "{ [tv0] -> [tv1] : tv0 - tv1 = 0 }"},
{"col1", "{ [tv0] -> [tv1] : tv0 - tv1 = 0 }"}
},
{
{"row3", "{ [tv0] -> [tv1] : tv1 - P0(row1(tv0), col1(tv0)) = 0 }"}
}));


comp->addStmt( new Stmt("AMCOO(n1) = ACOO(n )",
 "{ [n, n1] : n1 - P0(row1(n), col1(n)) = 0 && col1(n) - col3(n1) = 0 && row1(n) - row3(n1) = 0 && n >= 0 && n1 >= 0 && col1(n) >= 0 && row1(n) >= 0 && -n + NNZ - 1 >= 0 && -n1 + NNZ - 1 >= 0 && NC - col1(n) - 1 >= 0 && NR - row1(n) - 1 >= 0 }",
 "{ [tv0, tv1] -> [3, a1, 0, a3, 0, 0, 0] : tv0 - a1 = 0 && tv1 - a3 = 0 }",
{
{"ACOO", "{ [n, n1] -> [n] : n - n = 0 }"},
},
{
{"AMCOO", "{ [n, n1] -> [n1] : n1 - n1 = 0 }"}
}));

std::cout << comp->codeGen();
delete comp;



}
