//
// Created by Tobi Popoola on 8/30/20.
//

#include <gtest/gtest.h>
#include <iegenlib.h>
#include <CodeSynthesis.h>
#include <Permute.h>
#include <stdlib.h>
using namespace code_synthesis;
class CodeSynthesisUnitTest : public::testing::Test {
protected:
    virtual void SetUp() {}
    virtual void TearDown() {}
    void solveForOutput(std::string rel, std::string expected) {
        SCOPED_TRACE(rel);
        iegenlib::Relation* relO = new Relation(rel);
        iegenlib::Relation* res  = CodeSynthesis::solveForOutputTuple(relO);
        std::string resString = res->prettyPrintString();
        delete res;
        delete relO;
        EXPECT_EQ(resString,expected);
    }
};


TEST_F(CodeSynthesisUnitTest, TEST_EXPRESSION_TREE) {
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
/*
 * TODO: Switch test to for private functions in CodeSynthesis
*/
TEST_F(CodeSynthesisUnitTest, TEST_TERM_LIST) {
    Set* dense = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto list = CodeSynthesis::getTermList(dense);
    EXPECT_EQ((*list.begin())->prettyPrintString(dense->getTupleDecl()),"i");
}

TEST_F(CodeSynthesisUnitTest, TEST_CONTAINS_TERM) {
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

    EXPECT_EQ(true,CodeSynthesis::
              containsTerm(e1->getTermList(),tupleVarTerm));

    EXPECT_EQ(false,CodeSynthesis::
              containsTerm(e1->getTermList(),tupleVarTerm2));

    EXPECT_EQ("{ [i, j] : -a + b - row(-i, j) >= 0 }",conj1->prettyPrintString());
}



TEST_F(CodeSynthesisUnitTest, TEST_CONSTRAINT_TO_STATEMENT) {
    // UF(x,y) = n
    // {[x,y] -> [n]}
    TupleDecl td(3);
    td.setTupleElem(0,"x");
    td.setTupleElem(1,"y");
    td.setTupleElem(2,"n");
    Exp * e1 = new Exp();
    UFCallTerm* uf = new UFCallTerm("UF",2);
    Exp* eArg1 = new Exp();
    eArg1->addTerm(new TupleVarTerm(1,0));
    Exp* eArg2 = new Exp();
    eArg2->addTerm(new TupleVarTerm(1,1));
    uf->setParamExp(0,eArg1);
    uf->setParamExp(1,eArg2);
    e1->addTerm(uf);
    e1->addTerm(new TupleVarTerm(-1,2));
    e1->setEquality();

    SynthExpressionCase caseR =
        CodeSynthesis::GetUFExpressionSynthCase(e1,"UF",2,3);
    std::string statement = CodeSynthesis::
                            constraintToStatement(e1,"UF",td,caseR);
    EXPECT_EQ(statement,"UF->insert({x, y})");

    //{[ii,kk,jj,hr,hc] -> [k]}
    //rowptr(5 * ii + hr) = k
    TupleDecl td2(6);
    td2.setTupleElem(0,"ii");
    td2.setTupleElem(1,"kk");
    td2.setTupleElem(2,"jj");
    td2.setTupleElem(3,"hr");
    td2.setTupleElem(4,"hc");
    td2.setTupleElem(5,"k");
    e1 = new Exp();
    UFCallTerm* rowptr = new UFCallTerm("rowptr",1);
    Exp *paramExp = new Exp();
    paramExp->addTerm(new TupleVarTerm(5,0));
    paramExp->addTerm(new TupleVarTerm(1,3));
    rowptr->setParamExp(0,paramExp);
    e1->addTerm(rowptr);
    e1->addTerm(new TupleVarTerm(-1,5));
    e1->setEquality();
    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"rowptr",5,6);
    statement = CodeSynthesis::constraintToStatement(e1,"rowptr",td2,caseR);
    EXPECT_EQ(statement,"rowptr->insert({5 ii + hr})");



    // col(colinv(5ii+hr,5jj+hc)) = 5jj + hc // case 2
    e1 = new Exp();
    UFCallTerm* colInv = new UFCallTerm("colinv",2);
    Exp* colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    Exp* colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);


    UFCallTerm* colUF = new UFCallTerm("col",1);
    Exp* colUFArg1 = new Exp();
    colUFArg1->addTerm(colInv);
    colUF->setParamExp(0,colUFArg1);
    e1->addTerm(colUF);
    e1->addTerm( new TupleVarTerm(-5,2));
    e1->addTerm( new TupleVarTerm(-1,4));
    e1->setEquality();

    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"col",5,6);
    statement = CodeSynthesis::constraintToStatement(e1,"col",td2,caseR);

    EXPECT_EQ(statement,"col(colinv(5 ii + hr, 5 jj + hc))=5 jj + hc");
    delete e1;

    // rowptr(5ii + hr) <= colinv(5ii+hr,5jj+hc)
    // case 3
    colInv = new UFCallTerm("colinv",2);
    colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);

    rowptr = new UFCallTerm(-1,"rowptr",1);
    paramExp = new Exp();
    paramExp->addTerm(new TupleVarTerm(5,0));
    paramExp->addTerm(new TupleVarTerm(1,3));
    rowptr->setParamExp(0,paramExp);

    e1 = new Exp();
    e1->setInequality();
    e1->addTerm(rowptr);
    e1->addTerm(colInv);

    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"rowptr",5,6);
    statement = CodeSynthesis::constraintToStatement(e1,"rowptr",td2,caseR);

    EXPECT_EQ(statement,
              "rowptr(5 ii + hr) = min("
              "rowptr(5 ii + hr),colinv(5 ii + hr, 5 jj + hc))");
    delete e1;

    // rowptr(5ii + hr) >= colinv(5ii+hr,5jj+hc) + 1
    // case 4
    colInv = new UFCallTerm(-1,"colinv",2);
    colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);

    rowptr = new UFCallTerm("rowptr",1);
    paramExp = new Exp();
    paramExp->addTerm(new TupleVarTerm(5,0));
    paramExp->addTerm(new TupleVarTerm(1,3));
    rowptr->setParamExp(0,paramExp);

    e1 = new Exp();
    e1->setInequality();
    e1->addTerm(rowptr);
    e1->addTerm(colInv);
    e1->addTerm(new Term(-1));


    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"rowptr",5,6);
    statement = CodeSynthesis::constraintToStatement(e1,"rowptr",td2,caseR);
    EXPECT_EQ(statement,
              "rowptr(5 ii + hr) = max(rowptr(5 ii + hr),"
              "colinv(5 ii + hr, 5 jj + hc) + 1)");
    delete e1;

}



TEST_F(CodeSynthesisUnitTest, EXECUTION_SCHEDULE_SYNTHESIS) {
    Set* set = new Set("{[n,k]: n <= 0 and k >= M}");
    Relation * rel = CodeSynthesis::getExecutionSchedule(set,0);
    EXPECT_EQ("{ [tv0, tv1] -> [0, tv3, 0, tv5, 0] : tv0 - tv3 = 0 && tv1 - tv5 = 0 }",
              rel->prettyPrintString());
    delete rel;
    rel = CodeSynthesis::getExecutionSchedule(set,1);
    EXPECT_EQ("{ [tv0, tv1] -> [1, tv3, 0, tv5, 0] : tv0 - tv3 = 0 && tv1 - tv5 = 0 }",
              rel->prettyPrintString());

    delete set;
    set = new Set("{[n,k,q]: n <= 0 and k >= M and q <= 20}");
    rel = CodeSynthesis::getExecutionSchedule(set,0);
    EXPECT_EQ("{ [tv0, tv1, tv2] -> [0, tv4, 0, tv6, 0, tv8, 0] :"
              " tv0 - tv4 = 0 && tv1 - tv6 = 0 && tv2 - tv8 = 0 }",
              rel->prettyPrintString());
    delete set;
    delete rel;

    set = new Set("{ [tv0, tv1, tv2]: tv0 >= N}");
    rel = CodeSynthesis::getExecutionSchedule(set,0);

    EXPECT_EQ("{ [tv0, tv1, tv2] -> [0, tv4, 0, tv6, 0, tv8, 0] :"
              " tv0 - tv4 = 0 && tv1 - tv6 = 0 && tv2 - tv8 = 0 }",
              rel->prettyPrintString());
    delete rel;
    delete set;
}

TEST_F(CodeSynthesisUnitTest, CASE_TEST) {
    // UF(x,y) = n
    // {[x,y] -> [n]}
    Exp * e1 = new Exp();
    UFCallTerm* uf = new UFCallTerm("UF",2);
    Exp* eArg1 = new Exp();
    eArg1->addTerm(new TupleVarTerm(1,0));
    Exp* eArg2 = new Exp();
    eArg2->addTerm(new TupleVarTerm(1,1));
    uf->setParamExp(0,eArg1);
    uf->setParamExp(1,eArg2);
    e1->addTerm(uf);
    e1->addTerm(new TupleVarTerm(-1,2));
    e1->setEquality();

    SynthExpressionCase caseR =
        CodeSynthesis::GetUFExpressionSynthCase(e1,"UF",2,3);
    EXPECT_EQ(caseR,CASE1);
    delete e1;

    //{[ii,kk,jj,hr,hc] -> [k]}
    //rowptr(5 * ii + hr) = k
    e1 = new Exp();
    UFCallTerm* rowptr = new UFCallTerm("rowptr",1);
    Exp *paramExp = new Exp();
    paramExp->addTerm(new TupleVarTerm(5,0));
    paramExp->addTerm(new TupleVarTerm(1,3));
    rowptr->setParamExp(0,paramExp);
    e1->addTerm(rowptr);
    e1->addTerm(new TupleVarTerm(-1,5));
    e1->setEquality();
    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"rowptr",5,6);
    EXPECT_EQ(caseR,CASE1);
    delete e1;


    // col(colinv(5ii+hr,5jj+hc)) = 5jj + hc // case 2
    e1 = new Exp();
    UFCallTerm* colInv = new UFCallTerm("colinv",2);
    Exp* colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    Exp* colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);


    UFCallTerm* colUF = new UFCallTerm("col",1);
    Exp* colUFArg1 = new Exp();
    colUFArg1->addTerm(colInv);
    colUF->setParamExp(0,colUFArg1);
    e1->addTerm(colUF);
    e1->addTerm( new TupleVarTerm(-5,2));
    e1->addTerm( new TupleVarTerm(-1,4));
    e1->setEquality();
    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"col",5,6);

    EXPECT_EQ(caseR,CASE2);
    delete e1;

    // rowptr(5ii + hr) <= colinv(5ii+hr,5jj+hc)
    // case 3
    colInv = new UFCallTerm("colinv",2);
    colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);

    rowptr = new UFCallTerm(-1,"rowptr",1);
    paramExp = new Exp();
    paramExp->addTerm(new TupleVarTerm(5,0));
    paramExp->addTerm(new TupleVarTerm(1,3));
    rowptr->setParamExp(0,paramExp);

    e1 = new Exp();
    e1->setInequality();
    e1->addTerm(rowptr);
    e1->addTerm(colInv);

    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"rowptr",5,6);

    EXPECT_EQ(caseR, CASE3);
    delete e1;

    // rowptr(5ii + hr + 1) >= colinv(5ii+hr,5jj+hc) + 1
    // case 4
    colInv = new UFCallTerm(-1,"colinv",2);
    colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);

    rowptr = new UFCallTerm("rowptr",1);
    paramExp = new Exp();
    paramExp->addTerm(new TupleVarTerm(5,0));
    paramExp->addTerm(new TupleVarTerm(1,3));
    rowptr->setParamExp(0,paramExp);

    e1 = new Exp();
    e1->setInequality();
    e1->addTerm(rowptr);
    e1->addTerm(colInv);
    e1->addTerm(new Term(-1));


    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"rowptr",5,6);
    EXPECT_EQ(caseR, CASE4);
    delete e1;

    // UF(n) = y
    // {[x,y] -> [n]}
    e1 = new Exp();
    uf = new UFCallTerm("UF",1);
    eArg1 = new Exp();
    eArg1->addTerm(new TupleVarTerm(1,2));
    uf->setParamExp(0,eArg1);
    e1->addTerm(uf);
    e1->addTerm(new TupleVarTerm(-1,1));
    e1->setEquality();
    std::cerr << e1->toString() << "\n";
    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"UF",2,3);
    EXPECT_EQ(caseR,CASE5);
    delete e1;
}


TEST_F(CodeSynthesisUnitTest, EXPRESSION_DATA_ACCESS_SYNTHESIS) {
    // UF(x,y) = n
    // {[x,y] -> [n]}
    Exp * e1 = new Exp();
    UFCallTerm* uf = new UFCallTerm("UF",2);
    Exp* eArg1 = new Exp();
    eArg1->addTerm(new TupleVarTerm(1,0));
    Exp* eArg2 = new Exp();
    eArg2->addTerm(new TupleVarTerm(1,1));
    uf->setParamExp(0,eArg1);
    uf->setParamExp(1,eArg2);
    e1->addTerm(uf);
    e1->addTerm(new TupleVarTerm(-1,2));
    e1->setEquality();

    SynthExpressionCase caseR =
        CodeSynthesis::GetUFExpressionSynthCase(e1,"UF",2,3);
    auto writeAccess = CodeSynthesis::GetWrites("UF",e1,caseR,2);
    auto readAccess = CodeSynthesis::GetReads("UF",e1,caseR,2);
    EXPECT_EQ(1,writeAccess.size());
    EXPECT_EQ(0,readAccess.size());

    auto write = writeAccess[0];
    EXPECT_EQ("UF",write.first);
    EXPECT_EQ("{[tv0, tv1]->[0]}",write.second);
    delete e1;


    //{[ii,kk,jj,hr,hc] -> [k]}
    //rowptr(5 * ii + hr) = k
    e1 = new Exp();
    UFCallTerm* rowptr = new UFCallTerm("rowptr",1);
    Exp *paramExp = new Exp();
    paramExp->addTerm(new TupleVarTerm(5,0));
    paramExp->addTerm(new TupleVarTerm(1,3));
    rowptr->setParamExp(0,paramExp);
    e1->addTerm(rowptr);
    e1->addTerm(new TupleVarTerm(-1,5));
    e1->setEquality();
    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"rowptr",5,6);
    writeAccess = CodeSynthesis::GetWrites("rowptr",e1,caseR,5);
    readAccess = CodeSynthesis::GetReads("rowptr",e1,caseR,5);
    EXPECT_EQ(1,writeAccess.size());
    EXPECT_EQ(0,readAccess.size());

    write = writeAccess[0];
    EXPECT_EQ("rowptr",write.first);
    EXPECT_EQ("{[tv0, tv1, tv2, tv3, tv4]->[0]}",write.second);
    delete e1;


    // col(colinv(5ii+hr,5jj+hc)) = 5jj + hc + NR// case 2
    e1 = new Exp();
    UFCallTerm* colInv = new UFCallTerm("colinv",2);
    Exp* colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    Exp* colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);


    UFCallTerm* colUF = new UFCallTerm("col",1);
    Exp* colUFArg1 = new Exp();
    colUFArg1->addTerm(colInv);
    colUF->setParamExp(0,colUFArg1);
    e1->addTerm(colUF);
    //RHS 5jj + hc + NR
    e1->addTerm( new TupleVarTerm(-5,2));
    e1->addTerm( new TupleVarTerm(-1,4));
    e1->addTerm(new VarTerm(-1,"NR"));

    e1->setEquality();
    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"col",5,6);

    writeAccess = CodeSynthesis::GetWrites("col",e1,caseR,5);
    readAccess = CodeSynthesis::GetReads("col",e1,caseR,5);
    EXPECT_EQ(1,writeAccess.size());
    EXPECT_EQ(2,readAccess.size());

    write = writeAccess[0];
    EXPECT_EQ("col",write.first);
    EXPECT_EQ("{ [tv0, tv1, tv2, tv3, tv4] -> [tv5] :"
              " tv5 - colinv(5 tv0 + tv3, 5 tv2 + tv4) = 0 }",
              write.second);


    auto read1 = readAccess[0];
    EXPECT_EQ("NR",read1.first);
    EXPECT_EQ("{ [tv0, tv1, tv2, tv3, tv4] -> [tv5] : tv5 = 0 }",read1.second);

    auto read2 = readAccess[1];
    EXPECT_EQ("colinv",read2.first);
    EXPECT_EQ("{ [tv0, tv1, tv2, tv3, tv4] -> [tv5, tv6] :"
              " 5 tv0 + tv3 - tv5 = 0 && 5 tv2 + tv4 - tv6 = 0 }",read2.second);
    delete e1;


    //TODO: complete tests for Cases 3 & 4
    // rowptr(5ii + hr) <= colinv(5ii+hr,5jj+hc)
    // case 3
    colInv = new UFCallTerm("colinv",2);
    colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);

    rowptr = new UFCallTerm(-1,"rowptr",1);
    paramExp = new Exp();
    paramExp->addTerm(new TupleVarTerm(5,0));
    paramExp->addTerm(new TupleVarTerm(1,3));
    rowptr->setParamExp(0,paramExp);

    e1 = new Exp();
    e1->setInequality();
    e1->addTerm(rowptr);
    e1->addTerm(colInv);

    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"rowptr",5,6);

    delete e1;

    // rowptr(5ii + hr + 1) >= colinv(5ii+hr,5jj+hc) + 1
    // case 4
    colInv = new UFCallTerm(-1,"colinv",2);
    colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);

    rowptr = new UFCallTerm("rowptr",1);
    paramExp = new Exp();
    paramExp->addTerm(new TupleVarTerm(5,0));
    paramExp->addTerm(new TupleVarTerm(1,3));
    rowptr->setParamExp(0,paramExp);

    e1 = new Exp();
    e1->setInequality();
    e1->addTerm(rowptr);
    e1->addTerm(colInv);
    e1->addTerm(new Term(-1));


    caseR = CodeSynthesis::GetUFExpressionSynthCase(e1,"rowptr",5,6);
    delete e1;
}

TEST_F(CodeSynthesisUnitTest, TEST_REMOVE_SYMBOLIC_CONSTRAINTS) {

    iegenlib::Relation* map1 =
        new iegenlib::Relation("{[i,j]->[k]: i >= 0 and i < NR and"
                               " j >= 0 and j < NC and rowptr(i) <= k < rowptr(i+1)"
                               " and j = col1(k) and P(i,j) = k }");
    std::vector<std::string> ufs = { "rowptr","col1"};
    CodeSynthesis::RemoveSymbolicConstraints(ufs,map1);
    EXPECT_EQ("{ [i, j] -> [k] : k - P(i, j) = 0 && "
              "i >= 0 && j >= 0 && -i + NR - 1 >= 0"
              " && -j + NC - 1 >= 0 }",map1->prettyPrintString());
}

TEST_F(CodeSynthesisUnitTest,TEST_GETCASEDOMAIN) {
    iegenlib::Set* s =  new Set("{ [n,k] : k - P(row(n), col2(n)) = 0"
                                " && col1(k) - col2(n) = 0 && n >= 0"
                                " && col1(k) >= 0 && col2(n) >= 0 && row(n) >= 0"
                                " && k - rowptr(row(n)) >= 0 && NC - 1 >= 0 "
                                "&& NNZ - 1 >= 0 && NR - 1 >= 0 && P(row(n),"
                                " col2(n)) - rowptr(row(n)) >= 0 && -n + NNZ"
                                " - 1 >= 0 && -k + rowptr(row(n) + 1) - 1 >= 0"
                                " && NC - col1(k) - 1 >= 0 && NC - col2(n) - 1 >= 0"
                                " && NR - row(n) - 1 >= 0 && -P(row(n) + 1, col2(n))"
                                " + P(row(n), col2(n)) - 1 >= 0 && -P(row(n), col2(n))"
                                " + rowptr(row(n) + 1) - 1 >= 0 && -rowptr(row(n)) +"
                                " rowptr(row(n) + 1) - 1 >= 0 }");
    // creating constraint col1(k) - col2(n) = 0
    Exp * col2Const = new Exp();
    col2Const->setEquality();
    UFCallTerm *col1 = new UFCallTerm(1,"col1",1);
    Exp* arg1Exp = new Exp();
    arg1Exp->addTerm(new TupleVarTerm(1,1));
    col1->setParamExp(0,arg1Exp);

    UFCallTerm* col2 = new UFCallTerm(-1,"col2",1);
    Exp* argCol2Exp = new Exp();
    argCol2Exp->addTerm(new TupleVarTerm(1,0));
    col2->setParamExp(0,argCol2Exp);
    col2Const->addTerm(col1);
    col2Const->addTerm(col2);
    auto res = CodeSynthesis::GetCaseDomain("col2",s,col2Const,CASE2);

    EXPECT_EQ("{ [n, k] : k - P(row(n), col2(n)) = 0 && col1(k)"
              " - col2(n) = 0 && n >= 0 && col1(k) >= 0 && col2(n)"
              " >= 0 && row(n) >= 0 && k - rowptr(row(n)) >= 0 &&"
              " NC - 1 >= 0 && NNZ - 1 >= 0 && NR - 1 >= 0 &&"
              " P(row(n), col2(n)) - rowptr(row(n)) >= 0 && -n"
              " + NNZ - 1 >= 0 && -k + rowptr(row(n) + 1) - 1 >="
              " 0 && NC - col1(k) - 1 >= 0 && NC - col2(n) - 1 >="
              " 0 && NR - row(n) - 1 >= 0 && -P(row(n) + 1, col2(n))"
              " + P(row(n), col2(n)) - 1 >= 0 && -P(row(n), col2(n))"
              " + rowptr(row(n) + 1) - 1 >= 0 && -rowptr(row(n)) +"
              " rowptr(row(n) + 1) - 1 >= 0 }",res->prettyPrintString());
}



TEST_F(CodeSynthesisUnitTest,TEST_ADD_TO_DATASPACE) {
    std::vector<std::pair<std::string,std::string>> access =
    { { "row", "{[i,j]->[j]}"},{ "col", "{[i,j]->[i]}"}};
    Computation comp;
    CodeSynthesis::addToDataSpace(comp,access,"double");
    auto dataSpaces = comp.getUndelimitedDataSpaces();

    ASSERT_EQ(2,access.size());
    auto itR = dataSpaces.find("row");
    auto itC = dataSpaces.find("col");

    ASSERT_TRUE(itR != dataSpaces.end());
    ASSERT_TRUE(itC != dataSpaces.end());

    EXPECT_EQ("double*",dataSpaces["row"]);
    EXPECT_EQ("double*",dataSpaces["col"]);


}


TEST_F(CodeSynthesisUnitTest,TEST_GET_SUPPPORT_MACRO) {
    EXPECT_EQ("#define min(a,b) a < b ? a : b\n"
              "#define max(a,b) a > b ? a: b\n",
              code_synthesis::CodeSynthesis::getSupportingMacros());
}


TEST_F(CodeSynthesisUnitTest,TEST_MONOTONIC_DOMAIN) {
    iegenlib::Set * rowptrDomain = new iegenlib::Set("{[i]:0 <= i <= NR}");
    iegenlib::Set * rowptrRange = new iegenlib::Set("{[x]:0 <= x < NNZ}");

    Set* synthDomain =  code_synthesis::CodeSynthesis::GetMonotonicDomain(
                            "rowptr",iegenlib::Monotonic_Nondecreasing,
                            rowptrDomain);

    EXPECT_EQ("{ [e1, e2] : e1 >= 0 && e2 >= 0 &&"
              " -e1 + NR >= 0 && -e2 + NR >= 0 &&"
              " -e1 + e2 - 1 >= 0 }",
              synthDomain->prettyPrintString());

    delete synthDomain;
    delete rowptrDomain;
    delete rowptrRange;
}

TEST_F(CodeSynthesisUnitTest,TEST_MONOTONIC_STMT) {

    std::string synthStmt =  code_synthesis::CodeSynthesis::getMonotonicStmt(
                                 "rowptr",iegenlib::Monotonic_Nondecreasing);

    EXPECT_EQ("if ( not (rowptr(e1) <= rowptr(e2)))"
              "{rowptr(e2) = rowptr(e1);}",
              synthStmt);
}


TEST_F(CodeSynthesisUnitTest,TEST_MONOTONIC_DATA_ACCESSES) {

    auto reads =  code_synthesis::CodeSynthesis::getMonotonicReadAccess(
                      "rowptr",iegenlib::Monotonic_Nondecreasing);
    ASSERT_EQ(2, reads.size());
    EXPECT_EQ("rowptr", reads[0].first);
    EXPECT_EQ("{[e1,e2] -> [e2]}", reads[0].second);


    EXPECT_EQ("rowptr", reads[1].first);
    EXPECT_EQ("{[e1,e2] -> [e1]}", reads[1].second);

    auto writes =  code_synthesis::CodeSynthesis::getMonotonicWriteAccess(
                       "rowptr",iegenlib::Monotonic_Nondecreasing);
    ASSERT_EQ(1, writes.size());
    EXPECT_EQ("rowptr", writes[0].first);
    EXPECT_EQ("{[e1,e2] -> [e2]}", writes[0].second);
}

TEST_F(CodeSynthesisUnitTest, TEST_ADD_PERMUTATION_CONSTRAINT) {
    iegenlib::Relation* map1 =
        new iegenlib::Relation("{[i,j]->[k]: i >= 0 and i < NR and"
                               " j >= 0 and j < NC and rowptr(i) <= k < rowptr(i+1)"
                               " and j = col1(k)}");
    auto permutes =
        code_synthesis::CodeSynthesis::AddPermutationConstraint(map1);
    ASSERT_EQ(1, permutes.size());
    EXPECT_EQ("{ [i, j] -> [k] : j - col1(k) = 0 && k - P0(i, j) = 0 &&"
              " i >= 0 && j >= 0 && k - rowptr(i) >= 0 &&"
              " -i + NR - 1 >= 0 && -j + NC - 1 >= 0 &&"
              " -k + rowptr(i + 1) - 1 >= 0 }",
              map1->prettyPrintString());

    delete map1;

    map1 = new iegenlib::Relation("{[i,j]->[i,k]: i >= 0 and i < NR and"
                                  " j >= 0 and j < NC and rowptr(i) <= k < rowptr(i+1)"
                                  " and j = col1(k)}");
    permutes =
        code_synthesis::CodeSynthesis::AddPermutationConstraint(map1);
    ASSERT_EQ(2, permutes.size());
    EXPECT_EQ("{ [i, j] -> [i, k] : i - i = 0 && j - col1(k) = 0 &&"
              " i - P0(i, j) = 0 && k - P1(i, j) = 0 && i >= 0 && j >= 0 &&"
              " k - rowptr(i) >= 0 && -i + NR - 1 >= 0 && -j + NC - 1 >= 0 &&"
              " -k + rowptr(i + 1) - 1 >= 0 }",map1->prettyPrintString());


}


TEST_F(CodeSynthesisUnitTest, TEST_SUBSTITUTE_DIRECT_EQUALITIES) {
    iegenlib::Relation* map1 =
        new iegenlib::Relation("{[n]->[i,k]: i = row1(n) and i >= 0"
                               " and i < NR and j >= 0 and j < NC and"
                               " rowptr(i) <= k < rowptr(i+1) and j = col1(k)}");
    auto res1 = code_synthesis::CodeSynthesis::substituteDirectEqualities(map1);

    EXPECT_EQ("{ [n] -> [i, k] : i - row1(n) = 0 && j - col1(k) = 0"
              " && i >= 0 && j >= 0 && row1(n) >= 0 && k - rowptr(i) >= 0"
              " && k - rowptr(row1(n)) >= 0 && -i + NR - 1 >= 0 && -k +"
              " rowptr(i + 1) - 1 >= 0 && -k + rowptr(row1(n) + 1) - 1"
              " >= 0 && NC - j - 1 >= 0 && NR - row1(n) - 1 >= 0 }",
              res1->prettyPrintString());

    delete map1;
    delete res1;
}


TEST_F(CodeSynthesisUnitTest, TEST_PERMUTE) {
    std::vector<std::vector<int>> test = {{0,0},{0,1},{4,2},{8,0},{4,4}};
    int NNZ = test.size();
    int NR = 10;
    auto comparator = [&NR](const std::vector<int>& a, const std::vector<int>& b)
    {
        return a[0]* NR + a[1] < b[0] * NR + b[1];
    };
    auto linearization =  [&NR](const std::vector<int>& a)
    {
        return a[0]* NR + a[1];
    };
    Permute<int,decltype(linearization),decltype(comparator)>* P0 =
        new Permute<int,decltype(linearization),decltype(comparator)>
    (linearization,comparator,8);
    for(int i = 0 ; i  < NNZ ; i++) {
        P0->insert(test[i]);
    }
    int prevI = -1;
    for (int i = 0 ; i < NNZ; i++) {
        EXPECT_TRUE( prevI <= P0->getInv(i)[0]);
        prevI = P0->getInv(i)[0];
    }
    for (int i = 0 ; i < NNZ; i++) {
        EXPECT_EQ( i, P0->get( {P0->getInv(i)[0],P0->getInv(i)[1]}));
    }
    delete P0;
}


TEST_F(CodeSynthesisUnitTest, TEST_PERMUTESL) {
    std::vector<std::vector<int>> test = {{0,0},{0,1},{4,2},{8,0},{4,4}};
    int NNZ = test.size();
    int NR = 10;
    auto comparator = [&NR](const std::vector<int>& a, const std::vector<int>& b)
    {
        return a[0]* NR + a[1] < b[0] * NR + b[1];
    };
    auto linearization =  [&NR](const std::vector<int>& a)
    {
        return a[0]* NR + a[1];
    };
    PermuteSL<int,decltype(linearization),decltype(comparator)>* P0 =
        new PermuteSL<int,decltype(linearization),decltype(comparator)>
    (linearization,comparator);
    for(int i = 0 ; i  < NNZ ; i++) {
        P0->insert(test[i]);
    }
    std::cerr << "Permute Dump: \n"
              << P0->toString();
    int prevI = -1;
    for (int i = 0 ; i < NNZ; i++) {
        EXPECT_TRUE( prevI <= P0->getInv(i)[0]);
        prevI = P0->getInv(i)[0];
    }
    for (int i = 0 ; i < NNZ; i++) {
        EXPECT_EQ( i, P0->get( {P0->getInv(i)[0],P0->getInv(i)[1]}));
    }
    delete P0;
}

TEST_F(CodeSynthesisUnitTest, TEST_INVERSE_ITERATION_SPACE){
    Set* s = new Set("{[i,k,j,k1]: 0 <= i < NR && rowptr(i) <= k < rowptr(i+1)"
		    " && j = col(k) && k1 = P1(i,j)}");
    // P(i,j)

    UFCallTerm* p1Dim = new UFCallTerm(1,"P1",2);
    Exp* arg1 = new Exp();
    Exp* arg2 = new Exp();
    arg1->addTerm(new TupleVarTerm(0));
    arg2->addTerm(new TupleVarTerm(2));
    p1Dim->setParamExp(0,arg1);
    p1Dim->setParamExp(1,arg2);
    Set* inverseIterSpace = code_synthesis::
	    CodeSynthesis::GetInverseIterationSpace(s,p1Dim);
    EXPECT_EQ("{ [_n, _no, i, k, j, k1] : _n - k1 = 0 && _no - k = 0 &&"
	      " _no - P1MAP(_n) = 0 && i - P1DIM0(_no) = 0 &&"
	      " j - P1DIM1(_no) = 0 && _n >= 0 && -_n + P1SIZE"
	      " - 1 >= 0 }",inverseIterSpace->prettyPrintString());  
    delete s;
    delete inverseIterSpace;
    delete p1Dim;
    s = new Set("{[n,i,j,n1]: 0 <= i < NR && 0 <= j < NC && 0 <= n < NNZ "
		    " && i = row(n) && j=col(n) && n1 = P1(i,j)}");

    p1Dim = new UFCallTerm(1,"P1",2);
    arg1 = new Exp();
    arg2 = new Exp();
    arg1->addTerm(new TupleVarTerm(1));
    arg2->addTerm(new TupleVarTerm(2));
    p1Dim->setParamExp(0,arg1);
    p1Dim->setParamExp(1,arg2);
    inverseIterSpace = code_synthesis::
	    CodeSynthesis::GetInverseIterationSpace(s,p1Dim);


    EXPECT_EQ("{ [_n, _no, n, i, j, n1] : _n - n1 = 0 &&"
	      " _no - n = 0 && _no - P1MAP(_n) = 0 &&"
	      " i - P1DIM0(_no) = 0 && j - P1DIM1(_no) ="
	      " 0 && _n >= 0 && -_n + P1SIZE - 1 >= 0 }",
	      inverseIterSpace->prettyPrintString());  

    delete s;
    delete inverseIterSpace;
    delete p1Dim;
    // Test with aliases
    s = new Set("{[i,k,j,jj,k1]: j = jj && 0 <= i < NR && rowptr(i) <= k < rowptr(i+1)"
		    " && j = col(k) && k1 = P1(i,j)}");

    p1Dim = new UFCallTerm(1,"P1",2);
    arg1 = new Exp();
    arg2 = new Exp();
    arg1->addTerm(new TupleVarTerm(0));
    arg2->addTerm(new TupleVarTerm(2));
    p1Dim->setParamExp(0,arg1);
    p1Dim->setParamExp(1,arg2);
    inverseIterSpace = code_synthesis::
	    CodeSynthesis::GetInverseIterationSpace(s,p1Dim);
    EXPECT_EQ("{ [_no, _n, i, k, j, jj, k1] : _no - k = 0 && _n - k1 = 0 &&"
	      " _n - P1MAP(_no) = 0 && i - P1DIM0(_no) = 0 && j - jj = 0 &&"
              "j - P1DIM1(_no) = 0 && _no >= 0 && -_no + P1SIZE - 1 >= 0 }",
	      inverseIterSpace->prettyPrintString());  
    delete s;
    delete inverseIterSpace;
    delete p1Dim;
}


TEST_F(CodeSynthesisUnitTest, TEST_GET_RESOLVABLE_OUTPUT_TUPLE){
    iegenlib::Relation* map1 =
        new iegenlib::Relation("{[n]->[i,j,k]: i = row1(n) and i >= 0"
                               " and i < NR and j >= 0 and j < NC and"
                               " rowptr(i) <= k < rowptr(i+1) and P(i,j) = k"
			       " and j = col1(k)}");
    std::vector<int> res = code_synthesis::CodeSynthesis
	    ::GetResolvedOutputTuples(map1,{"rowptr","P"});
    ASSERT_EQ(res.size(),2);    
    EXPECT_EQ(1,res[0]);
    EXPECT_EQ(2,res[1]);
}
