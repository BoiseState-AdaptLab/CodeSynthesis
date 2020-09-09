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

std::list<Term*> MinimalSatisfiablity::getTermList(
        SparseConstraints * sparseConstraint){
    std::list<Term*> terms;

    for (auto conj : sparseConstraint->mConjunctions) {
        for (auto eq: conj->equalities()) {
            for (auto term: eq->getTermList()) {
                if (!term->isConst() &&
                    std::find(terms.begin(), terms.end(), term)
                    == terms.end()) {
                    terms.push_back(term);
                }
            }
        }

        for (auto ineq: conj->inequalities()) {
            for (auto term: ineq->getTermList()) {
                if (!term->isConst() &&
                    std::find(terms.begin(), terms.end(), term)
                    == terms.end()) {
                    terms.push_back(term);
                }
            }
        }
    }
    return terms;
}

std::list<Term*> MinimalSatisfiablity::evaluateUnknowns
    (Relation * transformRelation,Set* set){
    std::list<Term*> setTerms = getTermList(set);
    std::list<Term*> relationTerms= getTermList(transformRelation);
    // Terms in relations and not in set are unknowns.
    for(auto term :relationTerms){
        if (std::find(setTerms.begin(),
                      setTerms.end(),(term))!= setTerms.end()){
            relationTerms.remove(term);
        }
    }
    return relationTerms;

}


std::string Stmt::toString() const {
    std::stringstream ss;
    ss<<lhs->toString() << " = " << rhs->toString();
    return ss.str();
}