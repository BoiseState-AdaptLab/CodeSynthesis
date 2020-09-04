#include <iostream>
#include <smt/MinimalSatisfiablity.h>
using namespace code_synthesis::smt;
int main() {
    // row(n)
    Exp* e1 = new Exp();
    e1->addTerm(new VarTerm("n"));
    UFCallTerm * rowUF = new UFCallTerm("row",1);
    rowUF->setParamExp(0,e1);

    // row(n) - i = 0
    Exp* e2 = new Exp();
    e2->addTerm(rowUF);
    e2->addTerm(new VarTerm(-1,"i"));
    e2->setEquality();

    Stmt* stmt = MinimalSatisfiablity::synthStmt(e2,rowUF);

    //EXPECT_EQ("row(n) = i", stmt->toString());
//EXPECT_EQ("row(n) = i", stmt->toString());
    return 0;
}
