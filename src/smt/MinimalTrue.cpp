//
// Created by Tobi Popoola on 8/30/20.
//

#include "MinimalTrue.h"
#include <algorithm>
#include <unordered_set>
#include <utils/Utils.h>
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

/// Function flattens a sparse constraint : set, relation
/// to individual terms
/// \param sparseConstraint
/// \return
std::list<Term*> MinimalSatisfiablity::getTermList(
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

std::list<Term*> MinimalSatisfiablity::evaluateUnknowns
    (Relation * transformRelation,Set* set){
    std::list<Term*> setTerms = getTermList(set);
    std::list<Term*> relationTerms= getTermList(transformRelation);
    // Terms in relations and not in set are unknowns.
    auto i = relationTerms.begin();
    while (i!=relationTerms.end()){
        auto term = *i;

        if (std::find_if(setTerms.begin(),
                         setTerms.end(),[&term,&set](const Term *a){
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

Exp *MinimalSatisfiablity::getMinTrueExpr(Exp *expr) {
    Exp * copy = new Exp(*expr);
    if (copy->isEquality())
        return copy;
    copy->setEquality();
    return copy;

}

std::list<Term *> MinimalSatisfiablity::getDependents(Conjunction *conjunction, Term *term) {

    return std::list<Term *>();
}

int MinimalSatisfiablity::getTupleVarCount(std::list<Term *> &terms) {
    return std::count_if(terms.begin(),terms.end(),[](Term * a){
        return a->isTupleExp();
    });

}

/// This gets the list of all expressions in a conjunction.
/// Each expression in this list is newly allocated and
/// the caller is responsible for deallocating
/// \param conjunction
/// \return
std::list<Exp *> MinimalSatisfiablity::getExprs(Conjunction *conjunction) {
    std::list<Exp *> exps;
    for (auto eq: conjunction->equalities()) {
        exps.push_back(eq->clone());
    }

    for (auto ineq: conjunction->inequalities()) {
        exps.push_back(ineq->clone());
    }

    return  exps;
}

/// This function gets the domain of an unknown
/// Term in a relation
/// \param relation containing domain information
/// \param unknownTerm unkown term currently being investigated
/// \param unkownTerms unknown terms in the relation
/// \throw Exception if relation has no constraint
/// \return
Set *MinimalSatisfiablity::getDomain(Relation *relation, Term *unknownTerm,
                                     std::list<Term *>& unkownTerms) {
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
            auto it = constraints.begin();
            while( it!=constraints.end() ){
                auto exp = *it;
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
                    it = constraints.erase(it);
                    continue;
                }
                it++;
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
bool MinimalSatisfiablity::containsTerm(const std::list<Term *> &terms,
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

std::list <Term*> MinimalSatisfiablity::
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
std::list<T*> MinimalSatisfiablity::intersectLists(const std::list<T*> &a,
                                                  const std::list<T*> &b) {
    // TODO: Work on this to dereference a pointer
    std::list<T*> ret;
    std::unordered_set<T*>set;
    std::for_each(a.begin(),a.end(),[&set](T* h){
        set.insert(h);
    });
    std::for_each(b.begin(),b.end(),[&set,&ret](T* l){

        auto iter = std::find_if(set.begin(),set.end(),[&l](T* a){

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
bool MinimalSatisfiablity::compareAbsTerms(const Term *a, const Term *b) {
    Term* a_copy = a->clone();
    Term* b_copy = b->clone();
    a_copy->setCoefficient(1);
    b_copy->setCoefficient(1);
    bool res= (*a_copy) == (*b_copy);
    delete  a_copy;
    delete b_copy;
    return res;
}

/// This function creates a mapping from a sourceMap
/// to a destination map. The output arity of sourceMap
/// and output arity of the destinationMap must be equal.
/// \param sourceMaptoDense
/// \param destinationMaptToDense
/// \return
Relation *MinimalSatisfiablity::mapToNewSpace(Relation *sourceMap,
                                              Relation *destinationMap) {
    if (sourceMap->outArity()!=destinationMap->outArity()){
        throw assert_exception("source and destination map must have "
                               "the same output arity.");
    }
    Relation * destinationInverse = destinationMap->Inverse();
    Relation* result =  sourceMap->Compose(destinationInverse);
    delete  destinationInverse;
    return result;
}
