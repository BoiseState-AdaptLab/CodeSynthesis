//
// Created by Tobi Popoola on 8/30/20.
//

#include <gtest/gtest.h>
#include <iegenlib.h>
#include <smt/MinimalTrue.h>
#include <CodeSynthesis.h>


using namespace code_synthesis::smt;
using namespace code_synthesis;
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

TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_COO){
    // row(n)
    auto denseIterationSpace = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto mapFromDenseToCoo  = new Relation("{[i,j] -> [n]:"
                                         " row(n) = i and col(n) = j and  i >= 0 and "
                                         " i < NR and j >= 0 and j < NC}");
    auto list =MinimalSatisfiablity::evaluateUnknowns(mapFromDenseToCoo, denseIterationSpace);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromDenseToCoo->getTupleDecl()), "row(n)");
    EXPECT_EQ((*(++list.begin()))->
        prettyPrintString(mapFromDenseToCoo->getTupleDecl()), "col(n)");
    EXPECT_EQ(list.size(),2);
}

TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_CSR){
    auto denseIterationSpace = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto mapFromDenseToCsr  = new Relation("{[i,j] -> [k]:"
                                         " rowptr(i) <= k and rowptr(i+1) > k and col(k) = j and i >= 0 and "
                                         " i < NR and j >= 0 and j < NC}");
    auto list =MinimalSatisfiablity::evaluateUnknowns(mapFromDenseToCsr, denseIterationSpace);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromDenseToCsr->getTupleDecl()), "col(k)");
    EXPECT_EQ((*(++list.begin()))->
        prettyPrintString(mapFromDenseToCsr->getTupleDecl()), "k");
    EXPECT_EQ((*(++(++list.begin())))->
        prettyPrintString(mapFromDenseToCsr->getTupleDecl()), "rowptr(i)");
    EXPECT_EQ((*(++(++(++list.begin()))))->
        prettyPrintString(mapFromDenseToCsr->getTupleDecl()), "rowptr(i + 1)");
    EXPECT_EQ(list.size(),4);
}

TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_BCSR){
    //what about b? 
    //naming for i' and j'

    auto denseIterationSpace = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto mapFromDenseToBcsr  = new Relation("{[i,j] -> [b, ii, jj]:"
                                         " (ip, jp) = blocks(b) and ii = i - ip and jj = j - jp}");
    auto list =MinimalSatisfiablity::evaluateUnknowns(mapFromDenseToBcsr, denseIterationSpace);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromDenseToBcsr->getTupleDecl()), "blocks(b)");
    EXPECT_EQ((*(++list.begin()))->
        prettyPrintString(mapFromDenseToBcsr->getTupleDecl()), "( ip, jp )");
    EXPECT_EQ((*(++(++list.begin())))->
        prettyPrintString(mapFromDenseToBcsr->getTupleDecl()), "ii");
    EXPECT_EQ((*(++(++(++list.begin()))))->
        prettyPrintString(mapFromDenseToBcsr->getTupleDecl()), "ip");
    EXPECT_EQ((*(++(++(++(++list.begin())))))->
        prettyPrintString(mapFromDenseToBcsr->getTupleDecl()), "jj");
    EXPECT_EQ((*(++(++(++(++(++list.begin()))))))->
        prettyPrintString(mapFromDenseToBcsr->getTupleDecl()), "jp");
    EXPECT_EQ(list.size(),6);

}


TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_DIA){
    //what about d? 
    //why does it think that j is unknown? 

    auto denseIterationSpace = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto mapFromDenseToDia  = new Relation("{[i,j] -> [j, d]:"
                                         " (i - j) = diags(d)}");
    auto list =MinimalSatisfiablity::evaluateUnknowns(mapFromDenseToDia, denseIterationSpace);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromDenseToDia->getTupleDecl()), "j");
    EXPECT_EQ((*(++list.begin()))->
        prettyPrintString(mapFromDenseToDia->getTupleDecl()), "diags(d)");
    EXPECT_EQ(list.size(),2);
}

TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_ELL){
    //what about k? 
    //why does it think that i is unknown? 

    auto denseIterationSpace = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto mapFromDenseToEll  = new Relation("{[i,j] -> [i, k]:"
                                         " j = cols(i, k)}");
    auto list =MinimalSatisfiablity::evaluateUnknowns(mapFromDenseToEll, denseIterationSpace);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromDenseToEll->getTupleDecl()), "i");
    EXPECT_EQ((*(++list.begin()))->
        prettyPrintString(mapFromDenseToEll->getTupleDecl()), "cols(i, k)");
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
    auto denseIterationSpace = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto mapFromDenseToCoo  = new Relation("{[i,j] -> [n]:"
                                   " row(n) = i and col(n) = j and  i >= 0 and "
                                   " i < NR and j >= 0 and j < NC}");
    auto list =MinimalSatisfiablity::evaluateUnknowns(mapFromDenseToCoo, denseIterationSpace);
    Relation * restrictDomain = mapFromDenseToCoo->Restrict(denseIterationSpace);



    Set* rowDomain = MinimalSatisfiablity::
            getDomain(restrictDomain,(*list.begin()),list);
    
    EXPECT_EQ("{ [i, j] : i >= 0 && j >= 0 && Ad(i, j) - 1"
              " >= 0 && -i + NR - 1 >= 0 && -j + NC - 1 >= 0 }",
              rowDomain->prettyPrintString());

    // Column Domain
    Set* colDomain = MinimalSatisfiablity::
                    getDomain(restrictDomain,(*(++list.begin())),list);

    EXPECT_EQ("{ [i, j] : i >= 0 && j >= 0 && Ad(i, j) - 1"
              " >= 0 && -i + NR - 1 >= 0 && -j + NC - 1 >= 0 }",
              colDomain->prettyPrintString());




}

TEST_F (CodeSynthesisUnitTest, TEST_INSPECTOR_GENERATION_DENSE_TO_COO) {
    std::string denseSpace  = "{[i,j]: i >= 0 and i < NR and"
                              " j >= 0 and j < NC and Ad(i,j) > 0}";
    std::string mapFromDenseToCoo = "{[i,j] -> [n]:"
                          " row(n) = i and col(n) = j and  i >= 0 and "
                          " i < NR and j >= 0 and j < NC}";

    CodeSynthesis* synth = new CodeSynthesis(mapFromDenseToCoo, denseSpace);
    Computation *comp = synth->generateInspectorComputation();

    EXPECT_EQ(comp->getStmtSource(0),"row = newUF(1);");
    EXPECT_EQ(comp->getIterSpace(0), "{  }");
    EXPECT_EQ(comp->getExecSched(0), "{ [0, 0, 0, 0, 0] }");


    EXPECT_EQ(comp->getStmtSource(1),"row.insert(i);");
    EXPECT_EQ(comp->getIterSpace(1),
              "{ [i, j] : i >= 0 && j >= 0 && Ad(i, j) - 1 >= 0"
              " && -i + NR - 1 >= 0 && -j + NC - 1 >= 0 }");
    EXPECT_EQ(comp->getExecSched(1),
              "{ [i, j] -> [1, i, 0, j, 0] : i - i = 0 && j - j = 0 }");

    EXPECT_EQ(comp->getStmtSource(2),"col = newUF(1);");
    EXPECT_EQ(comp->getIterSpace(2), "{  }");
    EXPECT_EQ(comp->getExecSched(2),
              "{ [2, 0, 0, 0, 0] }");

    EXPECT_EQ(comp->getStmtSource(3),"col.insert(j);");
    EXPECT_EQ(comp->getIterSpace(3),
              "{ [i, j] : i >= 0 && j >= 0 && Ad(i, j) - 1"
              " >= 0 && -i + NR - 1 >= 0 && -j + NC - 1 >= 0 }");
    EXPECT_EQ(comp->getExecSched(3),
              "{ [i, j] -> [3, i, 0, j, 0] : i - i = 0 && j - j = 0 }");

    comp->printInfo();

}
