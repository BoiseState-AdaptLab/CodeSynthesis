//
// Created by Tobi Popoola on 8/30/20.
//

#include <gtest/gtest.h>
#include <iegenlib.h>
#include <CodeSynthesis.h>


using namespace code_synthesis;
class CodeSynthesisUnitTest : public::testing::Test{
protected:
    virtual void SetUp(){}
    virtual void TearDown(){}
    void solveForOutput(std::string rel, std::string expected){
        SCOPED_TRACE(rel);
	iegenlib::Relation* relO = new Relation(rel);
	iegenlib::Relation* res  = CodeSynthesis::solveForOutputTuple(relO);
	std::string resString = res->prettyPrintString();
	delete res;
	delete relO;
	EXPECT_EQ(resString,expected);
    }
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
/*
 * TODO: Switch test to for private functions in CodeSynthesis
*/
TEST_F(CodeSynthesisUnitTest, TEST_TERM_LIST){
    Set* dense = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto list = CodeSynthesis::getTermList(dense);
    EXPECT_EQ((*list.begin())->prettyPrintString(dense->getTupleDecl()),"i");
}

////////////////////////////////////////////
//                                        //
//          TEST_UNKNOWN_TERMS            //
//                                        //
////////////////////////////////////////////
/*
TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_COO){
    // row(n)
    std::string denseIterationSpace = 
        "{[i,j]: i >= 0 and i < NR and j >= 0 and j < NC and Ad(i,j) > 0}";

    std::string mapFromDenseToCoo  = 
	"{[i,j] -> [n]: row(n) = i and col(n) = j and  i >= 0 and "
        " i < NR and j >= 0 and j < NC}";
    CodeSynthesis codeSynthesis = new CodeSynthesis(
        mapFromDenseToCoo,denseIterationSpace);
    auto list =codeSynthesis->evaluateUnknowns(mapFromDenseToCoo, denseIterationSpace);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromDenseToCoo->getTupleDecl()), "row(n)");
    EXPECT_EQ((*(++list.begin()))->
        prettyPrintString(mapFromDenseToCoo->getTupleDecl()), "col(n)");
    EXPECT_EQ(list.size(),2);
    delete codeSynthesis;
}

TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_CSR){
    const char* denseIterationSpace =
        "{[i,j]: i >= 0 and i < NR and"
        " j >= 0 and j < NC and Ad(i,j) > 0}";
    const char* mapFromDenseToCsr  = 
	"{[i,j] -> [k]: rowptr(i) <= k and rowptr(i+1)"
        " > k and col(k) = j and i >= 0 and i < NR and"
	" j >= 0 and j < NC}";
    CodeSynthesis codeSynthesis = new CodeSynthesis(
        mapFromDenseToCsr,denseIterationSpace);
    auto list =codeSynthesis->evaluateUnknowns();
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromDenseToCsr->getTupleDecl()), "col(k)");
    EXPECT_EQ((*(++list.begin()))->
        prettyPrintString(mapFromDenseToCsr->getTupleDecl()), "k");
    EXPECT_EQ((*(++(++list.begin())))->
        prettyPrintString(mapFromDenseToCsr->getTupleDecl()), "rowptr(i)");
    EXPECT_EQ((*(++(++(++list.begin()))))->
        prettyPrintString(mapFromDenseToCsr->getTupleDecl()), "rowptr(i + 1)");
    EXPECT_EQ(list.size(),4);
    delete codeSynthesis;
}
/*
TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_BCSR){



    auto denseIterationSpace = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto mapFromDenseToBcsr  = new Relation("{[i,j] -> [b, ii, jj]:"
                                         " (ip, jp) = blocks(b) and ii = i - ip and jj = j - jp}");
    auto list =CodeSynthesis::evaluateUnknowns(mapFromDenseToBcsr, denseIterationSpace);
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

    //why does it think that j is unknown? 

    auto denseIterationSpace = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto mapFromDenseToDia  = new Relation("{[i,j] -> [j, d]:"
                                         " (i - j) = diags(d)}");
    auto list =CodeSynthesis::evaluateUnknowns(mapFromDenseToDia, denseIterationSpace);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromDenseToDia->getTupleDecl()), "diags(d)");
    EXPECT_EQ(list.size(),1);
}

TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_ELL){
    //why does it think that i is unknown? 

    auto denseIterationSpace = new Set("{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}");
    auto mapFromDenseToEll  = new Relation("{[i,j] -> [i, k]:"
                                         " j = cols(i, k)}");
    auto list =CodeSynthesis::evaluateUnknowns(mapFromDenseToEll, denseIterationSpace);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromDenseToEll->getTupleDecl()), "cols(i, k)");
    EXPECT_EQ(list.size(),1);
}

TEST_F(CodeSynthesisUnitTest, TEST_UNKNOWN_TERMS_COO_CSR){
    //why does it think that i and j are unknown? 

    auto cooIterationSpace = new Set("{[n]: n >= 0 and n < NNZ and"
                         " row(n) = i and col(n) = j and Ad(i,j) != 0}");
    auto mapFromCooToCsr  = new Relation("{[n] -> [k, i, j]:"
                                         " rowptr(i) <= k and k < rowptr(i + 1) and "
                                         "col_csr(k) = j and row(n) = i and col(n) = j}");

    auto list =CodeSynthesis::evaluateUnknowns(mapFromCooToCsr, cooIterationSpace);
    EXPECT_EQ((*list.begin())->
        prettyPrintString(mapFromCooToCsr->getTupleDecl()), "col_csr(k)");
    EXPECT_EQ((*(++list.begin()))->
        prettyPrintString(mapFromCooToCsr->getTupleDecl()), "k");
    EXPECT_EQ((*(++(++list.begin())))->
        prettyPrintString(mapFromCooToCsr->getTupleDecl()), "rowptr(i)");
    EXPECT_EQ((*(++(++(++list.begin()))))->
        prettyPrintString(mapFromCooToCsr->getTupleDecl()), "rowptr(i + 1)");

    EXPECT_EQ(list.size(),4);
}
*/

TEST_F(CodeSynthesisUnitTest, TEST_MINIMAL_TRUE){
    // -a + b >= 0
    Exp * e1 = new Exp();
    e1->addTerm(new VarTerm("b"));
    e1->addTerm(new VarTerm(-1,"a") );
    e1->setInequality();

    //-a + b = 0
    Exp * minTrueExp = CodeSynthesis::getMinTrueExpr(e1);

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

    EXPECT_EQ(true,CodeSynthesis::
            containsTerm(e1->getTermList(),tupleVarTerm));

    EXPECT_EQ(false,CodeSynthesis::
                containsTerm(e1->getTermList(),tupleVarTerm2));

    EXPECT_EQ("{ [i, j] : -a + b - row(-i, j) >= 0 }",conj1->prettyPrintString());
}


TEST_F (CodeSynthesisUnitTest, TEST_DOMAIN_EXTRACT){
    auto denseIterationSpace = "{[i,j]: i >= 0 and i < NR and"
                         " j >= 0 and j < NC and Ad(i,j) > 0}";
    auto mapFromDenseToCoo  = "{[i,j] -> [n]:"
                                   " row(n) = i and col(n) = j and  i >= 0 and "
                                   " i < NR and j >= 0 and j < NC}";
    CodeSynthesis * codeSynthesis = new CodeSynthesis(mapFromDenseToCoo,
        denseIterationSpace);
    auto list =codeSynthesis->evaluateUnknowns();

    Set* rowDomain = codeSynthesis->
            getDomain((*list.begin()),list);
    
    EXPECT_EQ("{ [i, j] : i >= 0 && j >= 0 && Ad(i, j) - 1"
              " >= 0 && -i + NR - 1 >= 0 && -j + NC - 1 >= 0 }",
              rowDomain->prettyPrintString());

    // Column Domain
    Set* colDomain = codeSynthesis->
                    getDomain((*(++list.begin())),list);

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

    EXPECT_EQ(comp->getStmt(0)->getStmtSourceCode(),"row = newUF(1);");
    EXPECT_EQ(comp->getStmt(0)->getIterationSpace()->prettyPrintString()
		    , "{  }");
    EXPECT_EQ(comp->getStmt(0)->getExecutionSchedule()->prettyPrintString()
		    , "{ [0, 0, 0, 0, 0] }");


    EXPECT_EQ(comp->getStmt(1)->getStmtSourceCode(),"row.insert(i);");
    EXPECT_EQ(comp->getStmt(1)->getIterationSpace()->prettyPrintString(), 
		    "{ [i, j] : i >= 0 && j >= 0 && Ad(i, j) - 1 >= 0"
		    " && -i + NR - 1 >= 0 && -j + NC - 1 >= 0 }");
    EXPECT_EQ(comp->getStmt(1)->getExecutionSchedule()->prettyPrintString(), 
		    "{ [i, j] -> [1, i, 0, j, 0] : i - i = 0 && j - j = 0 }");

    EXPECT_EQ(comp->getStmt(2)->getStmtSourceCode(),"col = newUF(1);");
    EXPECT_EQ(comp->getStmt(2)->getIterationSpace()->prettyPrintString(),
		    "{ [i, j] : i >= 0 && j >= 0 && Ad(i, j) - 1"
		    " >= 0 && -i + NR - 1 >= 0 && -j + NC - 1 >= 0 }");
    EXPECT_EQ(comp->getStmt(3)->getStmtSourceCode(),"col.insert(j);");
    EXPECT_EQ(comp->getStmt(3)->getIterationSpace()->prettyPrintString(),
		    "{ [i, j] : i >= 0 && j >= 0 && Ad(i, j) - 1"
		    " >= 0 && -i + NR - 1 >= 0 && -j + NC - 1 >= 0 }");
    EXPECT_EQ(comp->getStmt(3)->getExecutionSchedule()->prettyPrintString(), 
		    "{ [i, j] -> [3, i, 0, j, 0] : i - i = 0 && j - j = 0 }");
   
    comp->printInfo();
    EXPECT_EQ("",comp->codeGen());
}

TEST_F (CodeSynthesisUnitTest, TEST_INSPECTOR_GENERATION_DENSE_TO_CSR) {
    auto denseIterationSpace = "{[i,j]: i >= 0 and i < NR and"
                                       " j >= 0 and j < NC and Ad(i,j) > 0}";
    auto mapFromDenseToCsr  = 
        "{[i,j] -> [k]: rowptr(i) <= k and rowptr(i+1)"
	" > k and col(k) = j and i >= 0 and "
        " i < NR and j >= 0 and j < NC}";

    CodeSynthesis* synth = new CodeSynthesis(mapFromDenseToCsr , denseIterationSpace);
}

TEST_F (CodeSynthesisUnitTest, TEST_SOLVE_FOR_OUTPUT){
    std::string rel = 
        "{[i,j] -> [k]: A(i,j) > 0 and rowptr(i) <= k"
	" and k < rowptr(i+ 1) and col(k) =j and 0 <= i"
	" and i < NR and 0 <= j and j < NC}";
    std::string res =  
        "{[i,j] -> [k]: A(i,j) > 0 and rowptr(i) <= k"
	" and k < rowptr(i+ 1) and k = col_inv(i,j) and 0 <= i"
	" and i < NR and 0 <= j and j < NC}";
    solveForOutput(rel,res);
    

    solveForOutput("{[i,j] -> [n] : 0 <= n and n < NNZ and row(n) = i"
		   " and col(n) = j and 0 <= i and i < NR and 0 <= j "
		   " and j < NC} ",
		   "{[i,j] -> [n] : 0 <= n and n < NNZ and rowcol_inv(i,j) = n"
		   " and 0 <= i and i < NR and 0 <= j "
		   " and j < NC }");

}

TEST_F(CodeSynthesisUnitTest, TEST_CONSTRAINT_TO_STATEMENT){
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

    std::string statement = CodeSynthesis::constraintToStatement(e1,"UF",2);
    EXPECT_EQ(statement,"UF.insert(__tv0, __tv1)");

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
    statement = CodeSynthesis::constraintToStatement(e1,"rowptr",5);
    EXPECT_EQ(statement,"rowptr.insert(5 __tv0 + __tv3)");
    


    // col(colinv(5ii+hr,5jj+hc)) = 5jj + hc // case 2
    e1 = new Exp();
    UFCallTerm* colUF = new UFCallTerm("col",1);
    UFCallTerm* colInv = new UFCallTerm("colinv",2);
    Exp* colUFArg1 = new Exp();
    colUFArg1->addTerm(colInv);
    
    Exp* colInvArg1 = new Exp();
    colInvArg1->addTerm(new TupleVarTerm(5,0));
    colInvArg1->addTerm(new TupleVarTerm(1,3));
    Exp* colInvArg2 = new Exp();
    colInvArg2->addTerm(new TupleVarTerm(5,2));
    colInvArg2->addTerm(new TupleVarTerm(1,4));
    colInv->setParamExp(0,colInvArg1);
    colInv->setParamExp(1,colInvArg2);

    colUF->setParamExp(0,colUFArg1); 
    e1->addTerm(colUF);
    e1->addTerm( new TupleVarTerm(-5,2));
    e1->addTerm( new TupleVarTerm(-1,4));
    e1->setEquality(); 
    statement = CodeSynthesis::constraintToStatement(e1,"col",5);

    EXPECT_EQ(statement,"col(colinv(5 __tv0 + __tv3, 5 __tv2 + __tv4))=5 __tv2 + __tv4");
}


