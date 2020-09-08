//
// Created by Tobi Popoola on 8/30/20.
//

#include "MinimalTrue.h"
#include <algorithm>
using namespace code_synthesis::smt;
using iegenlib::Exp;
Stmt * MinimalSatisfiablity::synthStmt
        (Exp *constraint, Term *unknownTerm)
{
    Stmt* stmt = new Stmt();
    stmt->lhs = new Exp();
    Exp * rhsExpr  = constraint->solveForFactor(unknownTerm->clone());
    stmt->rhs = rhsExpr;
    return stmt;
}

std::list<Term*>  getTermList(
        SparseConstraints * sparseConstraint){
    std::list<Term*> terms;
    for(auto it =
            sparseConstraint->conjunctionBegin();it !=
            sparseConstraint->conjunctionEnd();it++){
        for(auto jt = (*it)->equalities().begin();
            jt != (*it)->equalities().end(); jt++ ){

            terms.merge((*jt)->getTermList());
        }
    }
    return terms;
}

std::list<Term*> MinimalSatisfiablity::evaluateUnknowns
    (Relation * transformRelation,Set* set){
    std::list<Term*> unknownTerms;

}


std::string Stmt::toString() const {
    std::stringstream ss;
    ss<<lhs->toString() << " = " << rhs->toString();
    return ss.str();
}