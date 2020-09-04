//
// Created by Tobi Popoola on 8/30/20.
//

#include "MinimalSatisfiablity.h"
#include <algorithm>
using namespace code_synthesis::smt;
using iegenlib::Exp;
Stmt * MinimalSatisfiablity::synthStmt
        (Exp *constraint, Term *unknownTerm)
{
    Stmt* stmt = new Stmt();
    stmt->lhs = std::unique_ptr<Exp>(new Exp());
    stmt->rhs = std::unique_ptr<Exp>(new Exp());
    auto mTerms  = constraint->getTermList();

    for (std::list<Term*>::iterator i=mTerms.begin(); i != mTerms.end(); ++i) {
        Term *t = *i;
        if (unknownTerm->toString()==t->toString()){
            stmt->lhs->addTerm(t->clone());
        }else{
            Term *nt = t->clone();
            nt->multiplyBy(-1);
            stmt->rhs->addTerm(nt);
        }
    }

    return stmt;
}
std::vector<Exp*> MinimalSatisfiablity::getRHS(
        Exp* currentExpr, Term *unknownTerm){
    // base case



}
std::string Stmt::toString() const {
    std::stringstream ss;
    ss<<lhs->toString() << " = " << rhs->toString();
    return ss.str();
}