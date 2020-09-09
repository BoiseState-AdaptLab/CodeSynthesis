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

TEST_F(MinimalUnitTest,TEST_TERM_LIST){
    Set* dense = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto list = MinimalSatisfiablity::getTermList(dense);
    EXPECT_EQ((*list.begin())->toString(),"i");
}

TEST_F(MinimalUnitTest,TEST_UNKNOWN_TERMS){
    // row(n)
    Set* dense = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    Relation * dense_coo  = new Relation("{[i,j] -> [n]:"
                                         " row(n) = i and col(n) = j and "
                                         " i < NR and j >= 0 and j < NC}");
    auto list =MinimalSatisfiablity::evaluateUnknowns(dense_coo,dense);
    for(auto const & l : list)
        std::cerr << "Hhh->" <<l->toString() << "\n";
}