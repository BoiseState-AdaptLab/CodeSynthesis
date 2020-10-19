//
// Created by Tobi Popoola on 10/7/20.
//

#include "CodeSynthesis.h"


#include <algorithm>
#include <tuple>
#include <unordered_set>
#include <map>
#include <utils/Utils.h>
using namespace code_synthesis;
using iegenlib::Exp;
Stmt * CodeSynthesis::synthStmt
        (Exp *constraint, Term *unknownTerm)
{
    Stmt* stmt = new Stmt();
    stmt->lhs = new Exp();
    Exp * rhsExpr  = constraint->solveForFactor(unknownTerm->clone());
    stmt->rhs = rhsExpr;
    return stmt;
}

/// Function flattens a sparse constraint : set, relation
/// to individual terms
/// \param sparseConstraint
/// \return
std::list<Term*> CodeSynthesis::getTermList(
        SparseConstraints * sparseConstraint){
    std::list<Term*> terms;
    for (auto conj : sparseConstraint->mConjunctions) {
        for (auto eq: conj->equalities()) {
            for (auto term: eq->getTermList()) {
                auto t = term->clone();
                t->setCoefficient(abs(t->coefficient()));
                if (t->isConst() || std::find_if(terms.begin(), terms.end(),
                                                 [&t](const Term *a) {
                                                     return (*a)==(*t);
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
                if (t->isConst() || std::find_if(terms.begin(), terms.end(),
                                                 [&t](const Term *a) {
                                                     return (*a)== (*t);
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

std::list<Term *> CodeSynthesis::evaluateUnknowns
        () {
    std::list<Term*> setTerms = getTermList(originalSpace);
    std::list<Term*> relationTerms= getTermList(mapToNewSpace);
    // Terms in relations and not in set are unknowns.
    auto i = relationTerms.begin();
    while (i!=relationTerms.end()){
        auto term = *i;

        if (std::find_if(setTerms.begin(),
                         setTerms.end(),[&term](const Term *a){
                    return (*term) == (*a);
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

Exp *CodeSynthesis::getMinTrueExpr(Exp *expr) {
    Exp * copy = new Exp(*expr);
    if (copy->isEquality())
        return copy;
    copy->setEquality();
    return copy;

}

std::list<Term *> CodeSynthesis::getDependents(Conjunction *conjunction, Term *term) {

    return std::list<Term *>();
}

int CodeSynthesis::getTupleVarCount(std::list<Term *> &terms) {
    return std::count_if(terms.begin(),terms.end(),[](Term * a){
        return a->isTupleExp();
    });

}

/// This gets the list of all expressions in a conjunction.
/// \param conjunction
/// \return
std::list<Exp *> CodeSynthesis::getExprs(Conjunction *conjunction) {
    std::list<Exp *> exps;
    for (auto eq: conjunction->equalities()) {
        exps.push_back(eq);
    }

    for (auto ineq: conjunction->inequalities()) {
        exps.push_back(ineq);
    }

    return  exps;
}

/// This function gets the domain of an unknown
/// Term in a relation
/// \param unknownTerm unkown term currently being investigated
/// \param unkownTerms unknown terms in the relation
/// \throw Exception if relation has no constraint
/// \return
Set *CodeSynthesis::getDomain(Term *unknownTerm,
                              std::list<Term *> &unkownTerms) {
    // Restrict the original space using the map to new space.
    Relation * relation = mapToNewSpace->Restrict(originalSpace);

    if(relation->getNumConjuncts()==0){
        throw assert_exception("getDomain: Relation should have constraints");
    }
    auto  conj = relation->
            conjunctionBegin();
    auto  end = relation->conjunctionEnd();
    std::list<Conjunction*> setConjunctions;

    while(conj!=end){
        Conjunction * c = (*conj);
        std::list<Exp*> constraints = getExprs(c);
        std::stack<Term*> dependentStack;
        std::unordered_set<Exp*> domainConstraints;
        std::set<TupleVarTerm*,TupleLexOrder> tupleConstraintVar;
        auto it = constraints.begin();
        while(it !=constraints.end() ){
            auto exp = *it;
            // Check if current expression contains the
            // unknown term. If it does select the depeneden terms
            if(containsTerm(exp->getTermList(),unknownTerm)){
                for(auto t : exp->getTermList()){

                    // If term is not part of the unkown and
                    // term is not some integer literal,
                    if (!containsTerm(unkownTerms,t) && !t->isConst()){
                        dependentStack.push(t);

                    }
                }
                // remove this constraint from the list
                // of constraint.
                it = constraints.erase(it);
                continue;
            }
            it++;
        }

        // This section uses the dependenceStack to
        // to get constraints for the domain been
        // extracted.
        while (!dependentStack.empty()){
            auto top = dependentStack.top();
            dependentStack.pop();
            auto topString = top->toString(false,false);
            if(auto tupleTerm = dynamic_cast<TupleVarTerm*>
            (top)){

                tupleConstraintVar.insert(tupleTerm);
            }
            auto it2 = constraints.begin();
            while(it2 != constraints.end() ){
                auto exp = *it2;
                auto st =exp->prettyPrintString(relation->getTupleDecl());
                // Check if the current expression does not
                // intersect with any of the unknown term.
                if (intersectLists(exp->getTermList(),unkownTerms).empty() &&
                    containsTerm(exp->getTermList(),top)){
                    for(auto q : exp->getTermList()){
                        // If term element q is current
                        // dependence on the top of the stack,
                        // take don't consider it.
                        if (compareAbsTerms(q, top)){
                            continue;
                        }
                        if (q->isUFCall()){
                            UFCallTerm* ufTerm = dynamic_cast<UFCallTerm*>(q);

                            // At this point we can add the terms in
                            // the expressions in the stack if it UF
                            // contains dependents
                            for(int i=0; i < ufTerm->numArgs();i++){
                                // Check if the argument contains
                                // the current dependent. If it does
                                // we pick up other arguments in that
                                // expression as well as the rest of the
                                // arguments
                                for (auto t : ufTerm->getParamExp(i)->getTermList()){
                                    if (!t->isConst() && !compareAbsTerms
                                            (t, top)){
                                        dependentStack.push(t);
                                    }
                                }
                            }

                        }else if (!q->isConst()){
                            dependentStack.push(q);
                        }
                    }
                    // Add expression to list of constraints
                    domainConstraints.insert(exp);

                    // Take the constraint out of the
                    it2 = constraints.erase(it2);
                    continue;
                }
                it2++;
            }


        }

        // Pack up the result of the algorithm as a
        // new conjunction
        TupleDecl tupleDecl (tupleConstraintVar.size());

        //TODO: Work on this below.
        // This is necessary to be able to remap previous
        // location of tuple variables to new location of
        // tuple variables in the new set / domain.
        std::vector<int> remapLocation(relation->getTupleDecl().size());


        int tupleID = 0;
        for (auto var : tupleConstraintVar){
            remapLocation[var->tvloc()] = tupleID;
            var->remapLocation(remapLocation);
            tupleDecl.setTupleElem(tupleID,var->prettyPrintString
                    (relation->getTupleDecl(),true));
            tupleID++;
        }
        Conjunction * resConj = new Conjunction(tupleDecl);
        for(auto it : domainConstraints){
            if ((*it).isEquality()){
                resConj->addEquality(it);
            }else{
                resConj->addInequality(it);
            }

        }
        setConjunctions.push_back(resConj);

        conj++;
    }

    Set* set = new Set(setConjunctions.front()->arity());
    for(auto c : setConjunctions){
        set->addConjunction(c);
    }
    return set;


}


/// Check if a term is contained in a list of terms
/// \param terms list of terms
/// \param term  been searched for
/// \return true if term is in there and returns false otherwise.
bool CodeSynthesis::containsTerm(const std::list<Term *> &terms,
                                        const Term *term) {

    return std::find_if(terms.begin(),
                        terms.end(),[&term](Term *a){
                bool res = compareAbsTerms(a,term);
                if (!res && a->isUFCall()){
                    UFCallTerm* ufCallTerm = dynamic_cast<UFCallTerm*>(a);
                    std::list<Term *> ufTerms = getParamTermList(ufCallTerm);
                    return containsTerm(ufTerms,term);

                }
                return res;

            })!= terms.end();
}

std::list <Term*> CodeSynthesis::
getParamTermList(const UFCallTerm *ufCallTerm) {
    std::list<Term*> ufTerms;
    for(int i = 0; i < ufCallTerm->numArgs();i++ ){
        for(auto t : ufCallTerm->getParamExp(i)
                ->getTermList()){
            ufTerms.push_back(t);
        }
    }
    return ufTerms;
}


/// Intersects two lists. Ignores sign of coefficients of
/// the lists
/// \tparam T type of the list
/// \param a first list
/// \param b second list
/// \return returns a list intersection
template<typename T>
std::list<T*> CodeSynthesis::intersectLists(const std::list<T*> &a,
                                                   const std::list<T*> &b) {
    // TODO: Work on this to dereference a pointer
    std::list<T*> ret;
    std::unordered_set<T*>set;
    std::for_each(a.begin(),a.end(),[&set](T* h){
        set.insert(h);
    });
    std::for_each(b.begin(),b.end(),[&set,&ret,this](T* l){

        auto iter = std::find_if(set.begin(),set.end(),[&l,this](T* a){

            return compareAbsTerms(a,l);
        });
        if(iter!=set.end()){
            ret.push_back(*iter);
            set.erase(iter);
        }
    });
    return ret;
}
/// Compares terms absolutely without regards for
//  their coefficients.
/// \param a term a
/// \param b term b
/// \return true if abs(a) is equal to abs(b)
bool CodeSynthesis::compareAbsTerms(const Term *a, const Term *b) {
    Term* a_copy = a->clone();
    Term* b_copy = b->clone();
    a_copy->setCoefficient(1);
    b_copy->setCoefficient(1);
    bool res= (*a_copy) == (*b_copy);
    delete  a_copy;
    delete b_copy;
    return res;
}

CodeSynthesis::~CodeSynthesis() {
    delete mapToNewSpace;
    delete originalSpace;
}

Computation* CodeSynthesis::generateInspectorComputation() {
   // First extract unknowns the information provided.
   Computation * comp = new Computation();
   std::list<Term*> unknownTerms =
           evaluateUnknowns();

   // Extract statement, domain and data spaces for each unknown.
   int stmtID = 0;
   int maxSchedule = 0;
   // TODO: work on generating execution schedule
   std::map<int,std::list<std::string>> executionSchedules;
   for(auto t : unknownTerms){
       // Have statement for creating unknownTerm
       std::string allocStmt = getAllocationStmt(t);
       comp->stmtsInfoMap[stmtID] =
               StmtInfo(allocStmt,"{[]}","{[]->[]}",
                            std::vector<std::pair<std::string,std::string>>(),
                            std::vector<std::pair<std::string,std::string>>());




       // Have statement for inserting into unknown Term.
       stmtID++;
       Set * insertStmtDomain = getDomain(t,unknownTerms);

       comp->stmtsInfoMap[stmtID] =
               StmtInfo("---",insertStmtDomain->prettyPrintString(),"{"+
                        insertStmtDomain->getTupleDecl().toString(true)+"->[]}",
                        std::vector<std::pair<std::string,std::string>>(),
                        std::vector<std::pair<std::string,std::string>>());





       stmtID++;
   }

   return comp;

}
/// This returns a string to allocate memory for an unknown
/// term.
/// \param unknownTerm
/// \return
std::string CodeSynthesis::getAllocationStmt(Term *unknownTerm) {
    if (!unknownTerm->isUFCall())
        throw assert_exception("Unknown term must be a UF");
    UFCallTerm * ufCallTerm = (UFCallTerm*) unknownTerm;
    std::stringstream  allocationString;
    allocationString <<ufCallTerm->name()<<" = newUF(" <<
        ufCallTerm->numArgs()<<");";

    return allocationString.str() ;
}
