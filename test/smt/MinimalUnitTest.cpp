//
// Created by Tobi Popoola on 8/30/20.
//

#include <gtest/gtest.h>
#include <iegenlib.h>
#include <smt/MinimalTrue.h>

using namespace code_synthesis::smt;
class MinimalUnitTest : public::testing::Test{
protected:
    virtual void SetUp(){}
    virtual void TearDown(){}
};


TEST_F(MinimalUnitTest,TEST_EXPRESSION_TREE){
    // a < b + 1 equivalent to b - a >= 0
    // -a + b
    Exp * e1 = new Exp();
    e1->addTerm(new VarTerm("b"));
    e1->addTerm(new VarTerm(-1,"a") );
    Conjunction* conj1 = new Conjunction(4);
    TupleDecl tdecl1(2);
    tdecl1.setTupleElem(0,"a");
    tdecl1.setTupleElem(1,"b");
    conj1->setTupleDecl(tdecl1);
    conj1->addInequality(e1);
    EXPECT_EQ("{ [a, b] : -a + b >= 0 }", conj1->prettyPrintString());
}
TEST_F(MinimalUnitTest,TEST_STATEMENT_SYNTHESIS_EQUALITY){
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

    Stmt* stmt = MinimalSatisfiablity::synthStmt(e2,rowUF->clone());

    EXPECT_EQ("row(n) = i", stmt->toString());
}