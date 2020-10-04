//
// Created by Tobi Popoola on 8/30/20.
//

#include <gtest/gtest.h>
#include <iegenlib.h>
#include <smt/MinimalTrue.h>

using namespace code_synthesis::smt;
class CodeSynthesisUnitTest : public::testing::Test{
protected:
    virtual void SetUp(){}
    virtual void TearDown(){}
};


TEST_F(CodeSynthesisUnitTest, TEST_EXPRESSION_TREE){
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

TEST_F(CodeSynthesisUnitTest, TEST_TERM_LIST){
    Set* dense = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto list = MinimalSatisfiablity::getTermList(dense);
    EXPECT_EQ((*list.begin())->prettyPrintString(dense->getTupleDecl()),"i");
}

TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS){
    // row(n)
    auto dense = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto dense_coo  = new Relation("{[i,j] -> [n]:"
                                         " row(n) = i and col(n) = j and  i >= 0 and "
                                         " i < NR and j >= 0 and j < NC}");
    auto list =MinimalSatisfiablity::evaluateUnknowns(dense_coo,dense);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(dense_coo->getTupleDecl()),"row(n)");
    EXPECT_EQ((*(++list.begin()))->
        prettyPrintString(dense_coo->getTupleDecl()),"col(n)");
    EXPECT_EQ(list.size(),2);
}

TEST_F(CodeSynthesisUnitTest, TEST_MINIMAL_TRUE){
    // -a + b >= 0
    Exp * e1 = new Exp();
    e1->addTerm(new VarTerm("b"));
    e1->addTerm(new VarTerm(-1,"a") );
    e1->setInequality();

    //-a + b = 0
    Exp * minTrueExp = MinimalSatisfiablity::getMinTrueExpr(e1);

    EXPECT_TRUE(minTrueExp->isEquality());

    // Testing with conjunction to have a better view.
    Conjunction* conj1 = new Conjunction(4);
    TupleDecl tdecl1(2);
    tdecl1.setTupleElem(0,"a");
    tdecl1.setTupleElem(1,"b");
    conj1->setTupleDecl(tdecl1);
    conj1->addEquality(minTrueExp);

    // Note, IEGENlib normalizes equations. For example,
    // -a + b = 0 becomes a - b = 0. How much does this
    // affect our minimally true statement ? Can we go
    // ahead with this ?
    EXPECT_EQ("{ [a, b] : a - b = 0 }", conj1->prettyPrintString());


}

TEST_F(CodeSynthesisUnitTest, TEST_CONTAINS_TERM){
    Exp * param1 = new Exp();
    param1->addTerm(new TupleVarTerm(-1,0));
    Exp * param2 = new Exp();
    param2->addTerm(new TupleVarTerm(1,1));
    Exp * e1 = new Exp();
    e1->addTerm(new VarTerm("b"));
    e1->addTerm(new VarTerm(-1,"a") );
    UFCallTerm * uf = new UFCallTerm(-1,"row",2);
    uf->setParamExp(0,param1);
    uf ->setParamExp(1,param2);
    e1->addTerm(uf);
    Conjunction* conj1 = new Conjunction(4);
    TupleDecl tdecl1(2);
    tdecl1.setTupleElem(0,"i");
    tdecl1.setTupleElem(1,"j");
    conj1->setTupleDecl(tdecl1);
    conj1->addInequality(e1);

    TupleVarTerm * tupleVarTerm = new TupleVarTerm(1,0);
    TupleVarTerm * tupleVarTerm2 = new TupleVarTerm(1,5);

    EXPECT_EQ(false,(*tupleVarTerm)==(*tupleVarTerm2));
    EXPECT_EQ(false,(*(*uf->getParamExp(0)->getTermList().begin())) ==(*tupleVarTerm2));

    EXPECT_EQ(true,MinimalSatisfiablity::
            containsTerm(e1->getTermList(),tupleVarTerm));

    EXPECT_EQ(false,MinimalSatisfiablity::
                containsTerm(e1->getTermList(),tupleVarTerm2));

    EXPECT_EQ("{ [i, j] : -a + b - row(-i, j) >= 0 }",conj1->prettyPrintString());
}


TEST_F (CodeSynthesisUnitTest, TEST_DOMAIN_EXTRACT){
    auto dense = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto dense_coo  = new Relation("{[i,j] -> [n]:"
                                   " row(n) = i and col(n) = j and  i >= 0 and "
                                   " i < NR and j >= 0 and j < NC}");
    auto list =MinimalSatisfiablity::evaluateUnknowns(dense_coo,dense);
    Relation * restrictDomain = dense_coo->Restrict(dense);

    Set* rowDomain = MinimalSatisfiablity::
            getDomain(restrictDomain,(*list.begin()),list);
    
    EXPECT_EQ("{ [i, j, tv2] : i >= 0 && j >= 0"
              " && Ad(i, j) - 1 >= 0 && -i + NR - 1"
              " >= 0 && -j + NC - 1 >= 0 }",rowDomain->prettyPrintString());

    // Column Domain
    Set* colDomain = MinimalSatisfiablity::
        getDomain(restrictDomain,(*(++list.begin())),list);
    EXPECT_EQ("{ [i, j, tv2] : i >= 0 && "
              "j >= 0 && Ad(i, j) - 1 >= 0 && "
              "-j + NC - 1 >= 0 }",colDomain->prettyPrintString());
}