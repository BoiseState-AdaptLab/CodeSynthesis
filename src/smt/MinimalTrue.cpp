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
                auto t = term->clone();
                t->setCoefficient(abs(t->coefficient()));
                if (t->isConst() || std::find_if(terms.begin(), terms.end(), [&t](const Term *a) {
                    return a->toString() == t->toString();
                }) != terms.end()){
                    delete t;
                    continue;
                }

                terms.push_back(t);
            }
        }

        for (auto ineq: conj->inequalities()) {
            for (auto term: ineq->getTermList()) {
                auto t = term->clone();
                t->setCoefficient(abs(t->coefficient()));
                if (t->isConst() || std::find_if(terms.begin(), terms.end(), [&t](const Term *a) {
                    return a->toString() == t->toString();
                }) != terms.end()){
                    delete t;
                    continue;
                }

                terms.push_back(t);

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
    auto i = relationTerms.begin();
    while (i!=relationTerms.end()){
        auto term = *i;
        auto termVar = term->prettyPrintString(transformRelation->getTupleDecl());
        if (std::find_if(setTerms.begin(),
                         setTerms.end(),[&termVar,&set](const Term *a){
                    auto aVar = a->prettyPrintString(set->getTupleDecl());
                    return termVar == aVar;
                })!= setTerms.end()){
            i = relationTerms.erase(i);
        }else{
            i++;
        }
    }
    return relationTerms;

}


std::string Stmt::toString() const {
    std::stringstream ss;
    ss<<lhs->toString() << " = " << rhs->toString();
    return ss.str();
}