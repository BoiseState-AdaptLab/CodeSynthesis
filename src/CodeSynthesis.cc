//
// Created by Tobi Popoola on 10/7/20.
//

#include "CodeSynthesis.h"

#include "Visitors.h"
#include <algorithm>
#include <tuple>
#include <unordered_set>
#include <map>
#include <utils/Utils.h>
#include <assert.h>
#include <sstream>
#include <list>
#include <numeric>
#include <bitset>
#include <stack>
// TODO: Work on memory management.

using namespace code_synthesis;
using iegenlib::Exp;
/// Function flattens a sparse constraint : set, relation
/// to individual terms
/// \param sparseConstraint
/// \return
std::list<Term*> CodeSynthesis::getTermList(
    SparseConstraints * sparseConstraint) {
    TermVisitor termVisitor;
    sparseConstraint->acceptVisitor(&termVisitor);
    return termVisitor.getTerms();
}






std::list<Term *> CodeSynthesis::getDependents(
    Conjunction *conjunction, Term *term) {

    return std::list<Term *>();
}

int CodeSynthesis::getTupleVarCount(std::list<Term *> &terms) {
    return std::count_if(terms.begin(),terms.end(),[](Term * a) {
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



/// Check if a term is contained in a list of terms
/// \param terms list of terms
/// \param term  been searched for
/// \return true if term is in there and returns false otherwise.
bool CodeSynthesis::containsTerm(const std::list<Term *> &terms,
                                 const Term *term) {

    return std::find_if(terms.begin(),
    terms.end(),[&term](Term *a) {
        bool res = compareAbsTerms(a,term);
        if (!res && a->isUFCall()) {
            UFCallTerm* ufCallTerm = dynamic_cast<UFCallTerm*>(a);
            std::list<Term *> ufTerms = getParamTermList(ufCallTerm);
            return containsTerm(ufTerms,term);

        }
        return res;

    })!= terms.end();
}

std::list <Term*> CodeSynthesis:: getParamTermList(
    const UFCallTerm *ufCallTerm) {
    std::list<Term*> ufTerms;
    for(int i = 0; i < ufCallTerm->numArgs(); i++ ) {
        for(auto t : ufCallTerm->getParamExp(i)->getTermList()) {
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
    std::for_each(a.begin(),a.end(),[&set](T* h) {
        set.insert(h);
    });
    std::for_each(b.begin(),b.end(),[&set,&ret](T* l) {

        auto iter = std::find_if(set.begin(),set.end(),[&l](T* a) {

            return compareAbsTerms(a,l);
        });
        if(iter!=set.end()) {
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
    delete sourceMapR;
    delete destMapR;
    delete composeRel;
    delete transRel;
    delete transRelExpanded;
}

Computation* CodeSynthesis::generateInspectorComputation() {
    Computation* inspector = new Computation();
    // transRelExpanded contains more constraints
    // that are certainly true and expands on all equalities
    // good candidate for generating statements.
    Conjunction * conj = *transRelExpanded->conjunctionBegin();
    std::list<iegenlib::Exp*> expList = getExprs(conj);
    // Convert compose relation to set.
    auto composeSet = composeRel->ToSet();

    // Convert trans relation to a set
    Set* transSet = transRel->ToSet();

    Set* noPermuteSet = new Set(*composeSet);

    RemoveSymbolicConstraints(permutes,noPermuteSet);
    
    // Reorder UF.
    UFCallTerm* reorderUF = NULL;


    int executionScheduleIndex  = 0;
    std::vector<std::string> removedPermutes;
    auto permIt = permutes.begin();
    while (permIt != permutes.end()) {
        auto permute  = *permIt;
        // TODO: Mege permutes with CASE1 and CASE2
        // CASE1, p0->insert({t1,t2})
        // CASE2, p0(t1,t2) = t3
        // Merged Permute
        // P0 = Permute(2);
        // p0->insert({t1,t2,t3})
        //
        bool permuteSelfRef = false;
        for(auto e : expList) {
            if (findCallTerm(e,permute)!=NULL) {
                auto caseP =
                    GetUFExpressionSynthCase(e,permute,
                                             transRel->inArity(),transRel->arity());
                if (caseP == SELF_REF) {
                    selfRefs.push_back({permute,e});
                    permuteSelfRef = true;
                }
            }
        }



        // Remove all self referential references to permute
        // from the transRel. This information is only important for
        // memory allocation.
        for(auto selfRef : selfRefs) {
            if (selfRef.first == permute) {
                RemoveConstraint(transRel,selfRef.second);
            }
        }

        // Gather all Permute PExpCandidate
        std::vector<std::pair<Exp*,SynthExpressionCase>>
                permuteExpCandidate;

        for(auto pExp : expList) {
            if (findCallTerm(pExp,permute)!=NULL) {
                auto caseP = GetUFExpressionSynthCase(pExp,permute,
                             transRel->inArity(),transRel->arity());
                if (caseP == CASE1 || caseP == CASE2)
                    permuteExpCandidate.push_back({pExp,caseP});

            }
        }

        // If there is just a single candidate. use
        // the generic case insert.
        if (permuteExpCandidate.size() == 1) {
            auto caseP = permuteExpCandidate[0].second;
            auto pExp = permuteExpCandidate[0].first;


            // Check if this permutation will affect reordering
            if (!permuteSelfRef &&
                    GeneratePermuteConditions(permute,composeRel,ufQuants).empty()) {
                // Remove every instance of this permute as
                // it does not affect ordering and should not
                // be used in synthesis processs
                RemoveSymbolicConstraints({permute},transSet);
                RemoveSymbolicConstraints({permute},transRel);
                RemoveSymbolicConstraints({permute},transRelExpanded);
                RemoveSymbolicConstraints({permute},composeRel);
                //Referesh expression list in case some
                //redundant permutes have been removed.
                conj = *transRelExpanded->conjunctionBegin();
                expList = getExprs(conj);
                removedPermutes.push_back(permute);
                continue;
            }

            std::string pStmt =
                constraintToStatement(pExp,
                                      permute,composeRel->getTupleDecl(),
                                      caseP);
            // Get Domain for P
            iegenlib::Set* pDomain =
                GetCaseDomain(permute,noPermuteSet,pExp,caseP);

            // remove constraints involving unknown UFs
            RemoveSymbolicConstraints(unknowns,pDomain);

            // Remove constraints involving permutes
            RemoveSymbolicConstraints(permutes,pDomain);

            // Get execution schedule
            iegenlib::Relation* pExecutionSchedule =
                getExecutionSchedule(
                    pDomain,executionScheduleIndex++);


            // Get reads and writes.
            auto writes =
                GetWrites(permute,pExp,caseP,pDomain->arity());

            auto reads =
                GetReads(permute,pExp,caseP,pDomain->arity());


            addToDataSpace(*inspector,reads, "double");
            // Writes to P is considered a single data space
            // but read from P is a 2d data space
            //addToDataSpace((*inspector),
            //				writes, "double");

            inspector->addStmt(new Stmt(pStmt,pDomain->
                                        prettyPrintString(),pExecutionSchedule->
                                        prettyPrintString(),reads,writes));

            UFCallTerm* pTerm = dynamic_cast<UFCallTerm*>(findCallTerm(pExp,permute));
            if(pTerm!=NULL) {
                CreateSortIRComponent(pTerm,inspector,executionScheduleIndex++);
                reorderUF = pTerm;
            }
        } else {
            //Remove all instances of permutes that also
            //fall into CASE 2 and references input tuple
            //they do not need to be generated. In the
            //future permutes that do not directly affect
            //the order of the output tensor will also be removed.
            RemoveSymbolicConstraints({permute},transSet);
            RemoveSymbolicConstraints({permute},transRel);
            RemoveSymbolicConstraints({permute},transRelExpanded);
            RemoveSymbolicConstraints({permute},composeRel);
            //Referesh expression list in case some
            //redundant permutes have been removed.
            conj = *transRelExpanded->conjunctionBegin();
            expList = getExprs(conj);
            removedPermutes.push_back(permute);
        }

        permIt++;
    }
    // Remove deleted permutes from
    // permute list.
    for(auto permute: removedPermutes) {
        auto it = std::find(permutes.begin(),
                            permutes.end(),permute);
        if (it != permutes.end())
            permutes.erase(it);

    }
    auto unknownsCopy = unknowns;
    // TryCount for halting the infinite loop
    // when an unknown UF cant be solved.
    int tryCount = 0 ;
    while(unknownsCopy.size() != 0) {
        if (tryCount >= MAX_TRIES) {
            break;
            //throw assert_exception("Synthesis Failed!");
        }
        // These are tuple terms that are assigned
        // to some known UFs and input tuple variable.
	// Get all tuple variables that are resolvable as
	// a function of input tuple variables
        std::vector<int> resolvedOutputTuples=
		GetResolvedOutputTuples(transRel,unknownsCopy);
        std::string currentUF = unknownsCopy.back();
        // Get all the list of viable candidates
        // for the current UF
	typedef std::pair<iegenlib::Exp*,SynthExpressionCase> ExpCasePair;
        std::list<ExpCasePair> expUfs;
        for(auto e : expList) {
            if(findCallTerm(e,currentUF)!=NULL) {
                // Get case a uf in a constraint falls into.
                auto ufCase =
                    GetUFExpressionSynthCase(e,
                     currentUF,transRel->inArity(),transRel->arity(),
		     resolvedOutputTuples);
		std::cout << "uf: " << currentUF << ", " << e->toString() << ", Case: "<< ufCase << "\n";
		// IF UF satisifies synthesis case
                if(ufCase !=UNDEFINED && ufCase != SELF_REF) {
                    expUfs.push_back({e,ufCase});
                }
            }
        }
        if (expUfs.size() == 0 ) {
            tryCount++;
            // Pop and put to the front of the list
            unknownsCopy.pop_back();
            unknownsCopy.insert(unknownsCopy.begin(),currentUF);
            continue;
        }
        // Check for prefered conditions among candidate
        // If equality & rhs is a function of known, we
        // care about such conditions and ignore other candidates
        auto it  = std::find_if(expUfs.begin(),expUfs.end(),
        [&](std::pair<iegenlib::Exp*,SynthExpressionCase> a) {
            if (a.second == CASE3 || a.second == CASE4)
                return false;
            Exp* constr = a.first;
            Term* term = findCallTerm(constr,currentUF);
            Term* tClone = term->clone();
            tClone->setCoefficient(1);
            Exp* solvedFor = constr->solveForFactor(tClone);
            if (solvedFor== NULL)
                return false;
            // Contains unknownsCopy
            bool containsUnknown = false;
            for(auto unknown: unknownsCopy) {
                if (findCallTerm(solvedFor,unknown) != NULL)
                    containsUnknown = true;
            }
            delete solvedFor;
	    bool domainBounded = IsDomainBoundedByUnknown(term,
			    unknownsCopy,composeRel);
            // This candidate is only viable if the current UF Term
            // domain is bounded by unknownsCopy.
            return !containsUnknown && domainBounded;
        });

        if ( it != expUfs.end()) {
            // Ignore and move forward if this fails and keep going
	    try{
	        CreateIRComponent(currentUF,inspector,executionScheduleIndex++, it->second,
                                  it->first,unknownsCopy,transSet,reorderUF);
            
	        // Remove from unknown list since this has been solved.
                unknownsCopy.pop_back();
                if (it->second == CASE5) {
                    permutes.push_back(currentUF);
                    UFCallTerm* ut = dynamic_cast<UFCallTerm*>(
                                     findCallTerm(it->first,currentUF));

                    if(ut!=NULL && queryMonoTypeEnv(currentUF)!= Monotonic_NONE) {
                        CreateSortIRComponent(ut,inspector,executionScheduleIndex++);
                    }
                }
	    } 
	    catch(...){
	    
	    }
            continue;
        }
        for(auto ufExpPair: expUfs) {
            // Avoid Generating code for Case5 that is bounded by
            if (ufExpPair.second == CASE5 ) {
                continue;
            }
            // Ignore and move forward if this fails and keep going
	    try{
               CreateIRComponent(currentUF,inspector,executionScheduleIndex++, 
			    ufExpPair.second,
                              ufExpPair.first,unknownsCopy,transSet,reorderUF);

	       // If what was solved for is case 5, case 2 or Case 1
	       // perform early optimization and dont generate code for any other
	       // constraints on this uf.
	       // early optimizations to remove multiple 
	       // statements doing the same thing
               if (ufExpPair.second == CASE5 || ufExpPair.second == CASE2 || ufExpPair.second == CASE5 ) {
                  break;
	       }
	    } 
	    catch(...){
	    
	    }
        }

        unknownsCopy.pop_back();
    }
    // Generate code to ensure universal constraint
    for(auto uf : unknowns) {
        iegenlib::MonotonicType type = iegenlib::queryMonoTypeEnv(uf);

        if(type==Monotonic_NONE)  continue;
        iegenlib::Set* domain = iegenlib::queryDomainCurrEnv(uf);


        UniQuantRule* uQ =  iegenlib::
                            getUQRForFuncDomainRange(uf);

        // Domain sorrunding statement for montonic
        // synthesis
        iegenlib::Set* stmtDomain =
            GetMonotonicDomain(uf,type,domain);


        // Get reads and writes accesses

        auto ufWrites =
            getMonotonicWriteAccess(uf,type);


        auto ufReads =
            getMonotonicReadAccess(uf,type);


        std::string monStmt =
            getMonotonicStmt(uf,type);
        // Get execution schedule
        iegenlib::Relation* execSched =
            getExecutionSchedule(
                stmtDomain,executionScheduleIndex++);

        inspector->addStmt(new Stmt(monStmt,stmtDomain->prettyPrintString()
                                    ,execSched->prettyPrintString()
                                    ,ufReads,ufWrites));
        delete stmtDomain;
        delete execSched;
    }
    // Skip copy code for now
    // CodeGen (RS2->S1(I)) - Data copy Code
    iegenlib::Set* copyDomain = composeRel->ToSet();

    // IF this space uses a reorder function
    // change the iteration spaace to loop through
    // the reorder function
    if (reorderUF){
       Set* inv = GetInverseIterationSpace(copyDomain,reorderUF);
       delete copyDomain;
       copyDomain = inv; 
    }
    std::string copyStmt = GetCopyStmt(sourceDataName,destDataName,destMapR,
                                       sourceMapR);
    iegenlib::Relation* execSchedule =
        getExecutionSchedule(
            copyDomain,executionScheduleIndex++);
    auto copyReads = getCopyReadAccess();
    auto copyWrites = getCopyWriteAccess();
    inspector->addStmt(new Stmt(copyStmt,copyDomain->prettyPrintString(),
                                execSchedule->prettyPrintString(),copyReads,
                                copyWrites));
    inspector->padExecutionSchedules();
    return inspector;
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

Term* CodeSynthesis::findCallTerm(Exp* exp, std::string ufName) {
    Term* res = nullptr;
    for(auto t : exp->getTermList()) {
        if ((t->isUFCall() && ((UFCallTerm*)t)->name() == ufName)
                ||
                (t->type() == "VarTerm" &&((VarTerm*)t)->symbol()
                 == ufName)) {
            res = t;
            break;
        }
    }
    return res;
}

std::string CodeSynthesis::
constraintToStatement(Exp* constraint,
                      std::string unknownUF, const TupleDecl&
                      tupDecl, SynthExpressionCase expCase) {
    Term* term = findCallTerm(constraint,unknownUF);
    if (term == NULL) {
        throw assert_exception("UFCallTerm must exist in expression");
    }

    int tupleSize = tupDecl.size();
    // Solve for UF term.
    Term* ufClone = term->clone();
    ufClone->setCoefficient(1);
    Exp* solvedUFConst = constraint->solveForFactor(ufClone);

    std::stringstream ss;
    if(expCase == CASE1) {
        UFCallTerm* ufTerm = (UFCallTerm*)term;
        ss << unknownUF << "->insert({";
        bool firstArg = true;
        for (int i = 0; i <ufTerm->numArgs(); ++i) {
            if (not firstArg) {
                ss << ", ";
            }
            if (ufTerm->getParamExp(i)) {
                ss << ufTerm->getParamExp(i)->prettyPrintString(tupDecl);
            }
            firstArg = false;
        }
        ss << "})";
    } else if (expCase == CASE2) {
        ss << term->prettyPrintString(tupDecl,true) << "="
           << solvedUFConst->prettyPrintString(tupDecl);

    } else if (expCase == CASE3) {
        ss << term->prettyPrintString(tupDecl,true) << " = "
           << "min(" << term->prettyPrintString(tupDecl,true)<< ","
           << solvedUFConst->prettyPrintString(tupDecl)
           << ")";

    } else if (expCase == CASE4) {

        ss << term->prettyPrintString(tupDecl,true) << " = "
           << "max(" << term->prettyPrintString(tupDecl,true)<< ","
           <<solvedUFConst->prettyPrintString(tupDecl)
           << ")";
    } else if (expCase == CASE5) {
        ss << unknownUF << "->insert({";
        ss << solvedUFConst->prettyPrintString(tupDecl);
        ss << "})";
    }
    else if(expCase == MERGECASE) {
        UFCallTerm* ufTerm = (UFCallTerm*)term;
        ss << unknownUF << "->insert({";
        bool firstArg = true;
        for (int i = 0; i <ufTerm->numArgs(); ++i) {
            if (not firstArg) {
                ss << ", ";
            }
            if (ufTerm->getParamExp(i)) {
                ss << ufTerm->getParamExp(i)->prettyPrintString(tupDecl);
            }
            firstArg = false;
        }
        ss << ",";
        ss << solvedUFConst->prettyPrintString(tupDecl);
        ss << "})";
    }

    delete solvedUFConst;
    return ss.str();
}



std::string CodeSynthesis::getFormattedTupleString(const std::list<std::string>& list) {
    std::stringstream ss;
    ss << "[";
    bool first = true;
    for(auto l : list) {
        if(first) {
            ss << l;
            first = false;
        } else {
            ss << "," <<l;
        }

    }
    ss << "]";
    return ss.str();
}

iegenlib::Set* CodeSynthesis::GetCaseDomain(std::string ufName,Set* s,
        Exp * constraint, SynthExpressionCase expCase) {

    Set * res;
    // Get the maximum tuple variable
    // present in this expression.
    Exp * constr = constraint->clone();
    // if this is CASE 5 evaluate the RHS instead.
    if (expCase == CASE5) {
        Term* t = findCallTerm(constr,ufName);
        if (t!= NULL) {
            Term* tClone = t->clone();
            tClone->setCoefficient(1);
            // No longer needed.
            delete constr;
            constr  = constraint->solveForFactor(tClone);
        }
    } else if (expCase == CASE1) {
        // For this case the domain
        // is the domain of u in UF(u)
        delete constr;
        // Create a constraint that contains
        // the current UF.
        constr = new Exp();
        Term* t = findCallTerm(constraint,ufName);
        constr->addTerm(t->clone());
        constr->setEquality();
    }
    res = new Set(*s);
    int maxTup = -1;
    for(int i = 0; i < s->arity(); i++) {
        TupleVarTerm t(1,i);
        if(constr->dependsOn(t)) {
            maxTup =  std::max(i,maxTup);
        }
    }

    if (maxTup == -1) {
        // In this case the domain is a a single instance.
        delete constr;
	delete res;
	return new Set("{[0]}");
    }
    // Project out tuple variables after
    // maxTuple.
    while(res->arity() - 1> maxTup ) {
        Set * temp = res->projectOut(res->arity() - 1);
        // We can project out anymore and we sshould just stop here
	// if temp  is null. This avoids complicated projection situation
	// and failure. This is a HACK
	if (temp == nullptr) break;
	delete res;
        res = temp;
    }
    delete constr;
    return res;

}

CodeSynthesis::CodeSynthesis(SparseFormat* source,
                             SparseFormat* dest) {
    destMapR = new Relation(dest->mapToDense);
    sourceMapR = new Relation(source->mapToDense);
    if(sourceMapR->outArity()!= destMapR->outArity()) {
        throw assert_exception("CodeSynthesis:: Format"
			" Descriptor map must"
                               " have the same output arity");
    }
    auto invDestMap = destMapR->Inverse();
    // Add the Permutation constraint.
    permutes = AddPermutationConstraint(invDestMap);

    composeRel = invDestMap->Compose(sourceMapR);
    transRel = composeRel->TransitiveClosure();
    // Expanded candidates for statement selections.
    transRelExpanded = substituteDirectEqualities(transRel);
    sourceDataName = source->dataName;
    destDataName = dest->dataName;
    
    sourceDataConstraint = source->dataConstraint;
    destDataConstraint = source->dataConstraint;

    sourceDataAccessMap = new Relation(source->dataAccess);
    destDataAccessMap = new Relation(dest->dataAccess);



    // Setup UFquantifiers
    iegenlib::setCurrEnv();
    for(auto uf: source->ufQuants) {
        Set* ufDomain = new Set(uf.domain);
        Set* ufRange = new Set(uf.range);
        iegenlib::appendCurrEnv(uf.name,ufDomain,ufRange,false,
                                uf.type);
    }

    for(auto uf: dest->ufQuants) {
        Set* ufDomain = new Set(uf.domain);
        Set* ufRange = new Set(uf.range);
        iegenlib::appendCurrEnv(uf.name,ufDomain,ufRange,false,
                                uf.type);
        ufQuants.push_back(uf);
    }


    for(auto known : dest->knowns) {
        knowns.push_back(known);
    }

    for(auto known : source->knowns) {
        knowns.push_back(known);
    }

    // Setup unknowns.

    iegenlib::StringIterator* iter = destMapR->getSymbolIterator();
    while (iter->hasNext()) {
        std::string symb = iter->next();
        // Exclude symbols that have already been specified to be
        // known.
        if (std::find(knowns.begin(),knowns.end(),symb) == knowns.end()) {
            unknowns.push_back( symb );
        }
    }
}


void CodeSynthesis::RemoveConstraint(SparseConstraints* sc, Exp *e ) {
    Conjunction* conj  = *sc->conjunctionBegin();
    if( e->isEquality()) {

        auto it = std::find_if(conj->equalities().begin(),
                               conj->equalities().end(),
        [e](Exp* e1) {
            return e1->toString() == (e)->toString();
        }
                              );
        if (it != conj->equalities().end()) {
            conj->equalities().erase(it);
        }
    } else if (e->isInequality()) {

        auto it = std::find_if(conj->inequalities().begin(),
                               conj->inequalities().end(),
        [e](Exp* e1) {
            return e1->toString() == (e)->toString();
        }
                              );
        if (it != conj->inequalities().end()) {
            conj->inequalities().erase(it);
        }
    }
}

std::vector<std::string> CodeSynthesis::
AddPermutationConstraint(Relation* rel) {
    // TODO: we can't add permutes whose
    // right hand side equals to some function
    // in the input tuple.
    std::vector<std::string> permutes;
    for(int i = 0; i < rel->outArity(); i++) {
        Exp* e  = new Exp();
        std::string name = PERMUTE_NAME + std::to_string(i);
        UFCallTerm* pUF = new UFCallTerm(1,name,rel->inArity());
        for(int i =0 ; i < rel->inArity(); i++) {
            TupleVarTerm* t = new TupleVarTerm(i);
            Exp* argi = new Exp();
            argi->addTerm(t);
            pUF->setParamExp(i,argi);
        }
        e->addTerm(pUF);
        TupleVarTerm* t = new TupleVarTerm(-1,rel->inArity()+ i);
        e->addTerm(t);
        e->setEquality();
        for(auto it = rel->conjunctionBegin(); it!=rel->conjunctionEnd();
                ++it) {
            (*it)->addEquality(e);
        }
        permutes.push_back(name);
    }
    return permutes;
}

iegenlib::Relation* CodeSynthesis::solveForOutputTuple(iegenlib::Relation* r) {
    assert(r->outArity()==1 && "Output arity must be 1");
    // Get Expressions involving output tuple variables
    TupleVarTerm * term= new TupleVarTerm(r->inArity());
    ExpTermVisitor expVisit(term);
    r->acceptVisitor(&expVisit);
    auto exps = expVisit.getExpressions();

    for(auto &e : exps) {
        auto solution = e->solveForFactor(term);
        if( not solution) {
            continue;
        }

    }
    return r;

}

/// Function creates an execution schedule from a set
/// and sets the outermost loop to position pos
/// Example:
///   {[n,k]: C1 ^ C2}, pos=1
///   {[n,k] -> [1,n,0,k,0]}
/// This functionality helps with code synthesis.
iegenlib::Relation* CodeSynthesis::getExecutionSchedule(iegenlib::Set* s,
        int pos) {
    int size = s->getTupleDecl().size();
    int execTupSize= size * 2;
    Relation* result = new Relation(size,execTupSize+1);
    Conjunction* conj = new Conjunction(size + execTupSize + 1,size);
    result->addConjunction(conj);
    Exp* e1 = new Exp();
    TupleVarTerm* posTup = new TupleVarTerm(1,size);
    Term * posTerm = new Term(-pos);
    e1->addTerm(posTerm);
    e1->addTerm(posTup);
    e1->setEquality();
    conj->addEquality(e1);

    for(int i = 0; i < execTupSize; i++) {
        Exp* e = new Exp();
        e->setEquality();
        TupleVarTerm* t2 = new TupleVarTerm(1,size+i+1);
        e->addTerm(t2);
        if ( i %2 == 0) {
            int originalPos = i / 2;
            //Construct an expression with tuple var terms.
            TupleVarTerm* t1 = new TupleVarTerm(-1,originalPos);
            e->addTerm(t1);
        }
        conj->addEquality(e);
    }
    conj->pushConstConstraintsToTupleDecl();
    return result;
}




// class DataAccessVisitor recursively visits
// an expression and extracts data access
// relations for UFCalls and Symbolic Constants.
class DataAccessVisitor: public Visitor {
private:
    std::vector<std::pair<std::string,std::string>> dAccess;
    int arity;
public:
    DataAccessVisitor(unsigned int arity): arity(arity) {}
    void preVisitUFCallTerm(UFCallTerm * t);
    void preVisitVarTerm(VarTerm* t);
    std::vector<std::pair<std::string,std::string>> getDataAccess() {
        return dAccess;
    }
    void clearDataAccesses() {
        dAccess.clear();
    }
};

void DataAccessVisitor::preVisitVarTerm(VarTerm* t) {
    // Constant variable term becomes
    // N
    // {"N", "{[Z]->[0]}"}
    std::string dataName = t->symbol();
    Relation* rel = new Relation(arity,1);
    Conjunction* conj = new Conjunction(arity+1,arity);
    rel->addConjunction(conj);
    Exp* e = new Exp();
    TupleVarTerm *tupTerm = new TupleVarTerm(1,arity);
    e->setEquality();
    e->addTerm(tupTerm);
    conj->addEquality(e);
    // Fix this issue
    // TODO
    std::string relString = rel->prettyPrintString();
    delete rel;
    dAccess.push_back({dataName,relString});
}
void DataAccessVisitor::preVisitUFCallTerm(UFCallTerm* t) {
    // The goal of this code section is to add
    // a constraint to the create relation such that
    // UF(x)
    //
    // {UF, {[Z]->[y]: y = x }
    // and x can contain tuple variables in Z
    std::string ufName = t->name();
    Relation* rel = new Relation(arity,t->numArgs());
    Conjunction* conj = new Conjunction(arity + t->numArgs(),arity);
    rel->addConjunction(conj);
    for(int i = 0; i < t->numArgs(); i++) {
        Exp* e = new Exp();
        e->setEquality();
        TupleVarTerm* tupTerm = new TupleVarTerm(1,arity+i);
        Exp* paramClone = t->getParamExp(i)->clone();
        paramClone->multiplyBy(-1);
        e->addTerm(tupTerm);
        e->addExp(paramClone);
        conj->addEquality(e);
    }
    std::string relString = rel->prettyPrintString();
    dAccess.push_back({ufName,relString});
    delete rel;
}

/// Function returns a list of write accesses.
std::vector<std::pair<std::string,std::string>> CodeSynthesis::GetWrites(
            std::string uf, iegenlib::Exp* constraint,
SynthExpressionCase expCase,int arity) {
    std::vector<std::pair<std::string,std::string>> result;
    Term* term = findCallTerm(constraint,uf);
    if (term == NULL) {
        throw assert_exception("UFCallTerm must exist in expression");
    }
    // Solve for UF term.
    Term* ufClone = term->clone();
    ufClone->setCoefficient(1);
    if (expCase == CASE1 || expCase == CASE5 || expCase == MERGECASE) {
        // This is for case 1 where
        // we have an insert abstraction
        TupleDecl tdl(arity);
        result.push_back({uf,"{"+tdl.toString(true)+
                          "->[0]}"});
    } else if (expCase ==  CASE2 ||
               expCase == CASE3 ||
               expCase == CASE4|| expCase == SELF_REF ) {
        DataAccessVisitor dV(arity);
        term->acceptVisitor(&dV);
        // We know there is only one UF write,
        // so we only pick up data accesses
        // that involves UF on LHS
        for(auto dAccess: dV.getDataAccess()) {
            if (dAccess.first == uf) {
                result.push_back(dAccess);
            }
        }
    }
    return result;
}

/// This function returns a list of read access made by an expression
/// based on its case.
std::vector<std::pair<std::string,std::string>> CodeSynthesis::GetReads(
            std::string uf,iegenlib::Exp* constraint,SynthExpressionCase expCase,
int arity) {
    Term* term = findCallTerm(constraint,uf);
    if (term == NULL) {
        throw assert_exception("UFCallTerm must exist in expression");
    }

    // Solve for UF term.
    Term* ufClone = term->clone();
    ufClone->setCoefficient(1);
    Exp* solvedUFConst = constraint->solveForFactor(ufClone);
    DataAccessVisitor dV(arity);
    solvedUFConst->acceptVisitor(&dV);
    auto result = dV.getDataAccess();
    dV.clearDataAccesses();
    // There might be some read accesses inside of the
    // ufTerm.
    // Do not check reads on the term if this
    // is CASE1 or CASE5. This is because
    // these cases use insert abstraction
    // and are not termed as reads
    if (expCase != CASE5 && expCase != CASE1) {
        term->acceptVisitor(&dV);
    }
    auto ufTermDataAccesses = dV.getDataAccess();
    for(auto dataAccess : ufTermDataAccesses) {
        if(dataAccess.first != uf) {
            result.push_back(dataAccess);
        }
    }
    return result;
}



/// Function returns case of expression as regards a UF.
SynthExpressionCase CodeSynthesis::GetUFExpressionSynthCase(Exp* constraint,
        std::string unknownUF, int inputArity, int tupleSize,
	const std::vector<int> resolvedOutputTuples) {
    SynthExpressionCase caseResult = UNDEFINED;
    Term* term = findCallTerm(constraint,unknownUF);
    if (term == NULL) {
        throw assert_exception("UFCallTerm must exist in expression");
    }
    // Solve for UF term.
    Term* ufClone = term->clone();
    ufClone->setCoefficient(1);
    Exp* solvedUFConst = constraint->solveForFactor(ufClone);

    if (solvedUFConst == NULL)
        return UNDEFINED;

    // if UF is self referential
    if (findCallTerm(solvedUFConst,unknownUF)!= NULL) {
        return SELF_REF;
    }
    
    // Compute UF and F dependency on output 
    // tuple
    // UF(x) (opp) F(x)
    bool FdependsOnOutput  = false;
    for(int i = inputArity; i < tupleSize; i++) {
       TupleVarTerm t(1,i);
       // check if this tuple is asigned to something
       // that solely depends on the input or some constant
       // then we know thist tuple is resolvable and can 
       // be assumed not to depend on output. We do this
       // by checking the resolvable output tuple list.
       bool outputTupleIsResolvable = std::find(resolvedOutputTuples.begin(),
		      resolvedOutputTuples.end(),i) !=
	       resolvedOutputTuples.end(); 
       if (solvedUFConst->dependsOn(t)
            && !outputTupleIsResolvable){
        	FdependsOnOutput = true;
                break;
           }
        }
        bool UFDependsOnOutput = false;

        if (term->isUFCall()) {
            UFCallTerm* ufTerm = (UFCallTerm*)term;
            for(int i = inputArity; i < tupleSize; i++) {
                TupleVarTerm t(1,i);
                bool outputTupleIsResolvable = std::find(resolvedOutputTuples.begin(),
		      resolvedOutputTuples.end(),i) != resolvedOutputTuples.end();
                for(int arg = 0; arg < ufTerm->numArgs(); arg++) {
                    if (ufTerm->getParamExp(arg)->dependsOn(t) && !outputTupleIsResolvable ) {
                        UFDependsOnOutput = true;
                        break;
                    }
                }
            }
        }
    
    std::stringstream ss;
    if (constraint->isEquality()) {
        if (!FdependsOnOutput && !UFDependsOnOutput) {
            caseResult = CASE2;
        } else if (UFDependsOnOutput && !FdependsOnOutput) {
            caseResult = CASE5;
        } else if (!UFDependsOnOutput && FdependsOnOutput) {
            caseResult = CASE1;
        }
    } else {
        // All cases in this section
        // arity(y) > arity(x)
        // Relax y_arity >= x_arity constraints
        if (!FdependsOnOutput && !UFDependsOnOutput) {
            // Case 3
            // UF(x) <= F(y)
            if (term->coefficient() < 0) {
                caseResult = CASE3;
            } else if (term->coefficient() > 0) {
                caseResult = CASE4;
            }
        }

    }
    delete solvedUFConst;
    return caseResult;
}

std::string CodeSynthesis::getSupportingMacros () {
    std::stringstream ss;
    ss << "#define min(a,b) a < b ? a : b\n"
       << "#define max(a,b) a > b ? a: b\n";
    return ss.str();
}

void CodeSynthesis::addToDataSpace(Computation& comp,
                                   std::vector<std::pair<std::string,std::string>> access,std::string baseType) {
    for(auto a : access) {
        Relation* accessRel = new Relation(a.second);
        int outArity = accessRel->outArity();
        delete accessRel;
        stringstream ss;
        ss << baseType;
        for(int i = 0 ; i < outArity; i++) ss << "*";
        std::string accessType = ss.str();
        comp.addDataSpace(a.first,accessType);
    }
}

/* Helper Class to return all constraints involving a certain
 * symbol name.
 * */

class GetAllConstraintOnSymbolName: public Visitor{
    std::vector<Exp*> constraints;
    std::stack<std::pair<int,bool>> s;
    std::vector<std::string> symbNames;
    public:
        GetAllConstraintOnSymbolName(std::vector<std::string> symbNames):
		symbNames(symbNames){}
        void preVisitExp(Exp* e) override;
	void postVisitExp(Exp* e) override;
	void preVisitUFCallTerm(UFCallTerm* ut) override;
	void preVisitVarTerm(VarTerm* vt) override;
	std::vector<Exp*> getConstraints();
};

void GetAllConstraintOnSymbolName::preVisitExp(Exp* e){
    if (!s.empty()){
       s.push({s.top().first+1,false});
    }else
       s.push({0,false});
}

void GetAllConstraintOnSymbolName::preVisitUFCallTerm(UFCallTerm* t){
    auto itU = std::find(symbNames.begin(),symbNames.end(),t->name());
    if(itU != symbNames.end()) { 
        // Pop whatever was stored for this expression 
	// and push true on the stack
	// Update every stack information backwards to true
	int currlevel = s.top().first;
	int count = 0;
	while(count <= currlevel){
	   s.pop();
	   count++;
	}
        count = 0;
	while(count <= currlevel){
	   s.push({count,true});
	   count++;
	}
    }
}

void GetAllConstraintOnSymbolName::preVisitVarTerm(VarTerm* vt){
    auto itU = std::find(symbNames.begin(),symbNames.end(),vt->symbol());
    if(itU != symbNames.end()) {
        // Pop whatever was stored for this expression 
	// and push true on the stack
	int currlevel = s.top().first;
	int count = 0;
	while(count <= currlevel){
	   s.pop();
	   count++;
	}
        count = 0;
	while(count <= currlevel){
	   s.push({count,true});
	   count++;
	}

    }
}

void GetAllConstraintOnSymbolName::postVisitExp(Exp* e){
    if (s.top().second){
        constraints.push_back(e);
    }
    s.pop();
}


std::vector<Exp*> GetAllConstraintOnSymbolName::getConstraints(){
    return constraints;
}


void CodeSynthesis::RemoveSymbolicConstraints(const std::vector<std::string>& symbNames,
        SparseConstraints* sc) {
    //Instantiate get all coonstraint on symbol name visitor	
    GetAllConstraintOnSymbolName GCV(symbNames);
    sc->acceptVisitor(&GCV);
    std::vector<Exp*> allConstraints = GCV.getConstraints();
    for (auto constraint : allConstraints){
        RemoveConstraint(sc,constraint); 
    }    
}

// Function returns read accesses for code generated
// in a montonic statement.
std::vector<std::pair<std::string,std::string>> CodeSynthesis::
getMonotonicReadAccess(std::string uf,MonotonicType type) {
    if (type == Monotonic_NONE) {
        throw assert_exception("getMonotonicStmt: none monotonic type"
                               " is not supported.");
    }
    // Reads from uf(e1) and uf(e2)
    std::vector<std::pair<std::string,std::string>>res;
    res.push_back({uf,"{[e1,e2] -> [e2]}"});
    res.push_back({uf,"{[e1,e2] -> [e1]}"});
    return res;
}


// Function returns write accesses for code generated
// in a montonic statement.
std::vector<std::pair<std::string,std::string>>
        CodeSynthesis::getMonotonicWriteAccess(std::string uf,
MonotonicType type) {
    if (type == Monotonic_NONE) {
        throw assert_exception("getMonotonicStmt: none monotonic type"
                               " is not supported.");
    }
    // Only writes to uf(e2) in all cases
    std::vector<std::pair<std::string,std::string>>res;
    res.push_back({uf,"{[e1,e2] -> [e2]}"});
    return res;
}



std::string CodeSynthesis::getMonotonicStmt(std::string uf,MonotonicType type) {
    if (type == Monotonic_NONE) {
        throw assert_exception("getMonotonicStmt: none monotonic type"
                               " is not supported.");
    }
    std::stringstream ss;
    ss <<  "if (";
    if (type == Monotonic_Increasing) {
        // if ( not (uf (e1) < uf(e2)) ){
        ss << " not (" << uf << "(e1) < "<< uf << "(e2))){";
        //   UF(e2) = UF(e1) + 1
        ss << uf << "(e2) = " << uf << "(e1) + 1;";
    } else if (type == Monotonic_Nondecreasing) {
        // if ( not (uf (e1) <= uf(e2)) ){
        ss << " not (" << uf << "(e1) <= "<< uf << "(e2))){";
        //   UF(e2) = UF(e1)
        ss << uf << "(e2) = " << uf << "(e1);";
    } else if (type == Monotonic_Decreasing) {
        // if ( not (uf (e1) > uf(e2)) ){
        ss << " not (" << uf << "(e1) > "<< uf << "(e2))){";
        //   UF(e2) = UF(e1) - 1
        ss << uf << "(e2) = " << uf << "(e1) - 1;";
    } else if (type ==Monotonic_Nonincreasing) {
        // if ( not (uf (e1) >= uf(e2)) ){
        ss << " not (" << uf << "(e1) >= "<< uf << "(e2))){";
        //   UF(e2) = UF(e1)
        ss << uf << "(e2) = " << uf << "(e1) - 1;";
    }
    ss << "}";
    return ss.str();
}




Set* CodeSynthesis::GetMonotonicDomain(std::string uf, MonotonicType type,
                                       Set* ufDomain) {
    if (ufDomain->arity()!=1) {
        throw assert_exception("GetMonotonicDomain: Uf domain"
                               " must have an arity of 1");
    }
    if (type == Monotonic_NONE) {
        throw assert_exception("GetMonotonicDomain: none monotonic type"
                               " is not supported.");
    }
    std::stringstream ss;
    TupleVarTerm tVar(0);

    std::list<Exp*> upperBounds = ufDomain->GetUpperBounds(tVar);
    std::list<Exp*> lowerBounds = ufDomain->GetLowerBounds(tVar);
    // Early optimization rather than have e1 < e2,
    // we can add e2 = e1 + 1 which holds true for the 
    // relationship between e1 and e2
    ss << "{[e1,e2]: e1 + 1 = e2 ";
    for(Exp* e : upperBounds) {
        ss << " && e1 <= " << e->toString();
        ss << " && e2 <= " << e->toString();
        // DELETE e Since we own e. I do not think this is a good idea
        // I suggest we refactor such that we don't own the
        // expression. TODO
        delete e;
    }


    for(Exp* e : lowerBounds) {
        ss << " && e1 >= " << e->toString();
        ss << " && e2 >= " << e->toString();

        // DELETE e Since we own e. I do not think this is a good idea
        // I suggest we refactor such that we don't own the
        // expression. TODO
        delete e;
    }
    ss << " }";
    return new Set(ss.str());
}



Exp* CodeSynthesis::getMonotonicDiff(std::string uf,Exp* ex) {
    // TODO: revisit this.
    UFCallTerm* ufTerm = (UFCallTerm*)findCallTerm(ex,uf);
    if (ufTerm == NULL) {
        throw assert_exception("getMonotonicDiff:"
                               " UFCallTerm must exist in expression");
    }
    // Solve for UF term.
    Term* ufClone = ufTerm->clone();
    ufClone->setCoefficient(1);
    Exp* solvedUFConst = ex->solveForFactor(ufClone);
    UFCallTerm* ufTerm2 = (UFCallTerm*)findCallTerm(solvedUFConst,uf);

    if (ufTerm2 == NULL) {
        throw assert_exception("getMonotonicDiff: Expression"
                               " is not self referential");
    }

    if(ufTerm2->numArgs() !=  1 || ufTerm->numArgs() != 1) {
        throw assert_exception("getMonotonicDiff: UFCallTerm must"
                               " have a single argument");
    }

    Exp* e = new Exp();
    // TODO: Figure out which what to subtract from what.
    Exp* e1Clone = ufTerm->getParamExp(0)->clone();
    Exp* e2Clone = ufTerm2->getParamExp(0)->clone();
    // eDiff = e2 - e1

    e1Clone->multiplyBy(-1);
    e->addExp(e2Clone);
    e->addExp(e1Clone);

    return e;

}


// Function returns statement for data copy
std::string CodeSynthesis::GetCopyStmt(std::string sourceDataName, std::string destDataName,
                                       Relation* destMap, Relation* sourceMap) {
    std::stringstream ss;
    ss << destDataName << "(";
    bool isFirst =true;
    for(int i =0; i < destMap->inArity(); i++) {
        if(isFirst) {
            ss << destMap->getTupleDecl().elemVarString(i);
            isFirst = false;
        } else {
            ss << "," <<  destMap->getTupleDecl().elemVarString(i);
        }
    }
    ss << ") = ";
    isFirst = true;
    ss << sourceDataName << "(";
    for(int i =0; i < sourceMap->inArity(); i++) {
        if(isFirst) {
            ss << sourceMap->getTupleDecl().elemVarString(i);
            isFirst = false;
        } else {
            ss << "," <<  sourceMap->getTupleDecl().elemVarString(i);
        }
    }
    ss << " )";
    return ss.str();
}

std::vector<std::pair<std::string,std::string>>
CodeSynthesis::getCopyWriteAccess() {
    std::vector<std::pair<std::string,std::string>>res;
    res.push_back({destDataName,destDataAccessMap->prettyPrintString()});
    return res;
}


std::vector<std::pair<std::string,std::string>>
CodeSynthesis::getCopyReadAccess() {
    // This is now the source data access
    std::vector<std::pair<std::string,std::string>>res;
    res.push_back({sourceDataName,sourceDataAccessMap->prettyPrintString()});
    return res;
}


std::string CodeSynthesis::generateFullCode(std::vector<int>& fuseStmts,
        int level) {
    Computation* comp = generateInspectorComputation();
    ReadReductionFusionOptimization(comp,fuseStmts,level);
     
    
    std::stringstream ss;
    ss << getSupportingMacros();
    for(auto permute : permutes ) {
        auto mergeIT = std::find_if(mergedPermutes.begin(),
                                    mergedPermutes.end(),
        [&permute](std::pair<std::string,int>& val) {
            return val.first == permute;
        });
        if (mergeIT != mergedPermutes.end()) {
            ss << "PermuteSimp<int> * " << permute <<
               " = new PermuteSimp<int>("<<mergeIT->second <<");\n";
            continue;
        }
        std::string permInit = "PermuteSimp<int> * "+permute +
                               " = new PermuteSimp<int>();\n";
        auto it = std::find_if(selfRefs.begin(),
                               selfRefs.end(),
        [&permute](std::pair<std::string,iegenlib::Exp*>& val) {
            return val.first == permute;
        });
	
	if(std::find(unknowns.begin(),unknowns.end(),permute)
			!= unknowns.end()){
	    // this is a growth function. A growth function is 
	    // one that fits into case 5 and part of unknown.
            ss <<"GrowthFunc<decltype(" << permute << "Comp)>* "
		   << permute << " = new GrowthFunc"<<
		   "<decltype(" << permute << "Comp)> ("<< permute
		   << "Comp);\n";
	    continue;
 
	}
        
	//Instantiate reorder stream constructor
	ss << "ReorderStream* "<< permute <<
		"= new ReorderStream("<< sourceMapR->outArity()
		<< ");\n";
        // Generate Comparator code.
        ss << "ComparatorInt "<< permute <<"Comp = ";
        ss << "[&](const int a,const int b){\n";
        if (it!=selfRefs.end()) {
            ss << GenerateSelfRefPermuteConditions(it->second, it->first);
        }
        ss << GeneratePermuteConditions(permute,composeRel,ufQuants);
        ss << "return false;\n";
        ss << "}; \n";
        ss << permute << "->setComparator("<<permute<<"Comp);\n";
        	
    }
    // Add Datamacros for Source and Destination
    // #define <sourceDataName>(i) <sourceDataName>[i]
    // #define <destDataName>(i) <destDataName>[
    bool isFirst = true;
    ss << "#define ";
    ss << sourceDataName << "(";
    for(int i =0; i < sourceDataAccessMap->inArity(); i++) {
        if(isFirst) {
            ss <<  sourceDataAccessMap->getTupleDecl().elemVarString(i);
            isFirst = false;
        } else {
            ss << ","  <<  sourceDataAccessMap->getTupleDecl().elemVarString(i);
        }
    }
    ss << ")";

    ss << " " << sourceDataName;
    for(int i=sourceDataAccessMap->inArity(); i < sourceDataAccessMap->arity(); i++) {
        TupleVarTerm tv (i);
        auto resolvedExp = getEqualExpr(sourceDataAccessMap,&tv);

        ss << "[" << (resolvedExp!=NULL? resolvedExp->prettyPrintString
                      (sourceDataAccessMap->getTupleDecl()):
                      sourceDataAccessMap->getTupleDecl().elemVarString(i))   << "]";
    }
    ss << "\n";
    isFirst = true;
    ss << "#define ";
    ss << destDataName << "(";
    for(int i =0; i < destDataAccessMap->inArity(); i++) {
        if(isFirst) {
            ss <<  destDataAccessMap->getTupleDecl().elemVarString(i);
            isFirst = false;
        } else {
            ss << ","  <<  destDataAccessMap->getTupleDecl().elemVarString(i);
        }
    }
    ss << ")";

    ss << " " << destDataName;
    for(int i=destDataAccessMap->inArity(); i < destDataAccessMap->arity(); i++) {
        TupleVarTerm tv (i);
        auto resolvedExp = getEqualExpr(destDataAccessMap,&tv);

        ss << "[" << (resolvedExp!=NULL? resolvedExp->prettyPrintString
                      (destDataAccessMap->getTupleDecl()):
                      destDataAccessMap->getTupleDecl().elemVarString(i))   << "]";
    }
    ss << "\n";
    std::string code = comp->codeGen();

    for(auto permute : permutes) {
        // This has to be modified to use
	// P_INV has to be done today.
	// TODO:
	// Replace P[][] access to P->get({});
        // #define P(i,j) P[i][j] becomes:
        // #define P(i,j) p->get({i,j})
        std::string p1 = permute+ "(";
        std::string p2 = permute;
        std::string p3 = permute+ "->get({";
        isFirst = true;
        for (int i =0 ; i < destMapR->outArity(); i++) {
            if(isFirst) {
                p1+="t" + std::to_string(i);
                p3+="t" + std::to_string(i);
                isFirst = false;
            } else {
                p1+=",t"+ std::to_string(i);
                p3+=",t"+ std::to_string(i);
            }
            p2+="[t" + std::to_string(i)+ "]";
        }
        p1 += ")";
        p3 += "})";
        std::string toReplace  = "#define "+ p1 + " "+ p2;
        std::string replacement = "#define "+ p1 + " "+ p3;
        Utils::replaceAllString(code, toReplace,replacement);


	// Replace
	// #define P1MAP(n) P1MAP[n]
	// with
	// #define P1MAP(n) P1->getMap(n)
	//
        toReplace   = "#define "+permute+"MAP(t0) "+permute+"MAP[t0]";
        replacement = "#define "+permute+"MAP(idx) "+permute+"->getMap(idx)";
        Utils::replaceAllString(code,toReplace,replacement);
        
	// Replace 
	// #define P1DIM0..N(idx) P1DIM0..N(idx) 
	// with 
	// #define P1DIM0..N(idx) P1->getDim(0..N,idx)
	for(int i = 0 ; i < sourceMapR->outArity(); i++){
 	    
            toReplace   = "#define "+permute+"DIM"+std::to_string(i)
		    +"(t0) "+permute+"DIM"+std::to_string(i)+"[t0]";
            replacement = "#define "+permute+"DIM"+std::to_string(i)
		    +"(idx) "+permute+"->getDim("
		    +std::to_string(i)+",idx)";
            Utils::replaceAllString(code,toReplace,replacement);
	}	
    }
    ss <<code;

    for(auto permute: permutes) {
        // Dont delete functions that are
        // part of unknowns.
        if (std::find(unknowns.begin(),unknowns.end(),permute)
                == unknowns.end()) {

            ss << "delete "<< permute << "; \n";
        }
    }

    delete comp;
    return ss.str();
}
// This header is depricated. and Perumte.h is used
std::string CodeSynthesis::GetSupportHeader() {
    std::stringstream ss;
    ss << "#ifndef SYNTH_HEADER\n";
    ss << "#define SYNTH_HEADER\n";
    ss << "#include <functional>\n";
    ss << "#include <algorithm>\n";
    ss << "#include <set>\n";
    ss << "#include <vector>\n";
    ss << "#include <assert.h>\n";
    ss << "#include <iostream>\n";
    ss << "#include <string>\n";
    ss << "#include <sstream>\n";
    ss<<"template <typename T,typename C = std::less<std::vector<T> > >   \n";
    ss<<"class Permutation{   \n";
    ss<<"private:    \n";
    ss<<"    std::set<std::vector<T>,C> d;   \n";
    ss<<"    int tupleSplit = 0;   \n";
    ss<<"public:   \n";
    ss<<"    Permutation(C comp):d(comp), tupleSplit(0){}   \n";
    ss<<"    Permutation(int tupleSplit): tupleSplit(tupleSplit) {}   \n";
    ss<<"    void insert(std::vector<T> tup){   \n";
    ss<<"        d.insert(tup);   \n";
    ss<<"    }   \n";
    ss<<"    int get(std::vector<T> tup){   \n";
    ss<<"        typename std::set<std::vector<T>>::iterator it;   \n";
    ss<<"    	if (tupleSplit == 0){   \n";
    ss<<"	    it = d.find(tup);   \n";
    ss<<"	}else{   \n";
    ss<<"	    it = std::find_if(d.begin(),d.end(),\n";
    ss <<	   "[this,&tup](const std::vector<T> &a){   \n";
    ss<<"			        for(int i=0; i < tupleSplit; i++){   \n";
    ss<<"				    if(a[i] != tup[i]) return false;   \n";
    ss<<"				}   \n";
    ss<<"				return true;   \n";
    ss<<"			    });   \n";
    ss<<"	}   \n";
    ss<<"	if (it == d.end()) {   \n";
    ss<<"	    assert(0 && \"Tuple get tuple does not exist\");   \n";
    ss<<"	}   \n";
    ss<<"	if (tupleSplit == 0) return std::distance(d.begin(),it);   \n";
    ss<<"	else return (*it)[tupleSplit];   \n";
    ss<<"    }   \n";
    ss<<"    std::string toString(){   \n";
    ss<<"	std::stringstream ss;   \n";
    ss<<"	for(int i = 0; i < d.size(); i++){   \n";
    ss<<"	    ss<< \"[\" << i << \"] => {\";   \n";
    ss<<"	    for(int j = 0; j  < d[i].size(); j++){   \n";
    ss<<"	        ss << d[i][j] << \",\";   \n";
    ss<<"	    }   \n";
    ss<<"	    ss << \"}\";   \n";
    ss<<"	}   \n";
    ss<<"	return ss.str();   \n";
    ss<<"    }   \n";
    ss<<"};   \n";
    ss<<"#endif   \n";
    return ss.str();
}



// Function checks if an expression depends on output
// tuple returns true if it does.
// \param arity   arity of the constraint
// \param inArity input arity of the constraint
// \param e       expression
// Params are not adopted
bool CodeSynthesis::dependsOnOutputTuple(int arity, int inArity,
        iegenlib::Exp*e) {

    for(int i = inArity; i < arity; i++) {
        TupleVarTerm tup(i);
        if (e->dependsOn(tup)) {
            return true;
        }
    }
    return false;
}

// This function substitutes direct equalities on output tuple
// if the rhs of the equality is a function of the input tuple.
// Example:
//     R {[n] -> [i,k]: row(n) = i ^ P(i,k) = u }
//     it generates
//     R {[n] -> [i,k]: row(n) = i ^ P(row(n),k) = u
Relation* CodeSynthesis::substituteDirectEqualities(Relation* rel) {
    Relation* res = new Relation(*rel);
    SubMap subMap;
    for(auto conj = res->conjunctionBegin();
            conj!= res->conjunctionEnd(); conj++) {
        auto itE = (*conj)->equalities().begin();
        while(itE != (*conj)->equalities().end()) {
            auto eq = (*itE);
            if(!eq) continue;
            for(int i =res->inArity(); i < res->arity(); i++) {
                TupleVarTerm tup(i);
                if (eq->dependsOn(tup)) {
                    Term* t = tup.clone();
                    Exp* e = eq->solveForFactor(t);
                    if (e && !dependsOnOutputTuple(res->arity(),
                                                   res->inArity(),
                                                   e)) {
                        subMap.insertPair(tup.clone(),e);
                    }
                }
            }
            itE++;
        }
    }
    res->substituteInConstraints(subMap);
    Relation* r = res->Intersect(rel);
    delete res;
    return r;
}

bool CodeSynthesis::IsTupleBoundedByUnknown(TupleVarTerm& t,
        SparseConstraints* sp,
        std::vector<std::string>& unknowns) {
    return false;
}

Exp* CodeSynthesis::getEqualExpr(SparseConstraints* sc, const Term* t) {
    Exp * res = NULL;
    for(auto conj =  sc->conjunctionBegin();
            conj!= sc->conjunctionEnd(); conj++) {
        for( auto exp: (*conj)->equalities()) {
            if (exp->dependsOn(*t)) {
                res = exp->solveForFactor(t->clone());
                break;

            }
        }
    }
    return res;
}


// class ExpComparatorVisitor generates expression
// for comparator.
// F(row(e1)) > H (row(e2))
//
// P->insert(col(n),row(n));
//
// row(n) is at position x[1] in P
// col(n) is at position x[0] in P
//
// x is determined by e1CompName or e2 comp name
//
class ExpComparatorVisitor: public Visitor {
private:
    UFCallTerm* permuteTerm=NULL;
    std::string e1CompName;
    std::string e2CompName;
    std::list<Exp*> equalityExps;
    std::stringstream ss;
public:
    explicit ExpComparatorVisitor(UFCallTerm* permuteTerm,
                                  std::string e1CompName, std::string e2CompName,
                                  std::list<Exp*> eqExps):
        permuteTerm(permuteTerm),e1CompName(e1CompName),
        e2CompName(e2CompName),
        equalityExps(eqExps) {}

    void postVisitTupleVarTerm(iegenlib::TupleVarTerm * t) override;
    void postVisitUFCallTerm(iegenlib::UFCallTerm* t) override;


};

std::string CodeSynthesis::GeneratePermuteConditions(std::string& permute,
        Relation* composeRel,std::vector<UFQuant>& ufQuants) {
    std::stringstream ss;
    // Check if this uf already has a ufquantifier
    auto ufIt = std::find_if(ufQuants.begin(),ufQuants.end(),
    [&](const UFQuant& a) {
        return a.name == permute;
    });
    // This only supports single parameter
    // UFs. should break if this is not a single
    // parameter UF.
    if (ufIt != ufQuants.end()) {
        iegenlib::Set* domain = iegenlib::queryDomainCurrEnv(permute);
        if (domain == NULL) return "";

        if (domain->arity() != 1) {
            delete domain;
            throw iegenlib::assert_exception("arity greater than"
                                             " 1 is not courrently supported");
        }
        delete domain;
        ss << "if(";
        MonotonicType type = ufIt->type;

        if (type == Monotonic_Increasing
                || type == Monotonic_Nondecreasing) {
            ss << "a[0] < b[0]";
        } else if (type == Monotonic_Decreasing
                   || type == Monotonic_Nonincreasing) {
            ss << "a[0] > b[0]";
        }
        ss << ")\n return true ; \n";
    }

    // Discover the range of the permute in
    // the relation. If equivalent to one of the UFs
    // enforce constraints from that UF
    Conjunction* c = *composeRel->conjunctionBegin();
    for(auto e : c->equalities()) {
        UFCallTerm *permTerm = NULL;
        if (permTerm = (UFCallTerm*) findCallTerm(e,permute)) {

            auto ufCase =
                GetUFExpressionSynthCase(e,
                                         permute,composeRel->inArity(),
                                         composeRel->arity());
            // Check what this permTerm is equal to
            // and get the upper and lower bounds
            Term* tClone = permTerm->clone();
            tClone->setCoefficient(1);
            Exp* solvedFor = e->solveForFactor(tClone);
            Term* tupTerm = solvedFor->getTerm();
            // Skip this if it does not only have one term
            if (tupTerm == NULL) continue;

            if (tupTerm->type() != "TupleVarTerm") continue;

            // Get range of term
            std::list<Exp*> lowerBounds =
                c->GetLowerBounds(*((TupleVarTerm*)tupTerm));
            std::list<Exp*> upperBounds =
                c->GetUpperBounds(*((TupleVarTerm*)tupTerm));

            Set * permuteRange = new Set(1);
            Conjunction* permRangeConj =
                *permuteRange->conjunctionBegin();
            // Search for
            TupleVarTerm tup(0);
            // Build permute range with upper and lower.
            // bounds.
            for(Exp* e : upperBounds) {
                Term* tupClone = tup.clone();
                tupClone->setCoefficient(-1);
                e->addTerm(tupClone);
                e->setInequality();
                permRangeConj->addInequality(e);
            }

            for(Exp* e : lowerBounds) {
                Term* tupClone = tup.clone();
                tupClone->setCoefficient(1);
                e->multiplyBy(-1);
                e->addTerm(tupClone);
                e->setInequality();
                permRangeConj->addInequality(e);
            }
            permuteRange->cleanUp();

            // Now we have a permute range, we need to
            // look for which UF satisfies => Domain(UF) = range(Permute)
            for (auto ufQuant : ufQuants) {
                Set* ufQuantDomain = new Set(ufQuant.domain);
                if (*ufQuantDomain == *permuteRange &&
                        ufQuant.rhsProperty!= "") {
                    Set* rhsProperty = new Set(ufQuant.rhsProperty);
                    auto rhsConj = *rhsProperty->conjunctionBegin();
                    // Now we create comparator constraints based on this.
                    ss << "if(";
                    bool first = true;
                    for (Exp* rhsExp : rhsConj->inequalities()) {
                        if (first) {
                            first = false;
                        }
                        else ss << " && ";
                        ss << rhsExp->prettyPrintString(rhsConj->
                                                        getTupleDecl());
                        ss << ">= 0";
                    }
                    for (Exp* rhsExp : rhsConj->equalities()) {
                        ss << " && ";
                        ss << rhsExp->prettyPrintString(rhsConj->
                                                        getTupleDecl());
                        ss << " == 0";

                    }
                    ss << ")\n";
                    // TODO: Replace with actual vectors comming
                    // into permutation.
                    ss << "return true;\n";
                }
            }
        }
    }
    return ss.str();
}


// This function agressively fuses loops with true dependency
// in order.
void CodeSynthesis::ReadReductionFusionOptimization(Computation* comp,
        std::vector<int>& fuseStmts, int level) {
    if (fuseStmts.size() == 0) return;
    // Statements with the same domain and does not write
    int fuseStart = *fuseStmts.begin();
    for (auto it = fuseStmts.begin() + 1; it!=fuseStmts.end() ; it++) {
	comp->fuse(fuseStart,*it,level);
	fuseStart = *it;
    }
}

// Remove statements that write to the same point in memory
void CodeSynthesis::RedundantStatementElimination(Computation* comp) {
    comp->deleteDeadStatements();
}

// Simplify constraints and also optimizes out statements
// involving such constraints.
void CodeSynthesis::ConstraintSimplification(Computation* comp) {

}

std::string CodeSynthesis::GenerateSelfRefPermuteConditions(Exp* e, std::string& permute) {
    UFCallTerm* ut =(UFCallTerm*) findCallTerm(e,permute);
    Term* cloneUT = ut->clone();
    cloneUT->setCoefficient(1);
    Exp* solveE = e->solveForFactor(cloneUT);
    UFCallTerm* ut2 = (UFCallTerm*)findCallTerm(solveE,permute);;
    if (ut2 == NULL) return "";
    std::stringstream ssP;
    for(int k = 0; k < ut->numArgs(); k++) {
        Exp* kLHSExp = ut->getParamExp(k);
        Exp* kRHSExp = ut2->getParamExp(k);
        // Get the difference between lhs and rhs
        // to see which is bigger.
        Exp* diff = new Exp();
        diff->addExp(kLHSExp);
        diff->multiplyBy(-1);
        diff->addExp(kRHSExp);
        int coeff = 0;
        Term * t = diff->getTerm();
        if (t!=NULL) {
            coeff = t->coefficient();
        }
        if (coeff ==0) continue;
        ssP << "if (a[" << k  << "] ";
        if (coeff < 0 ) {
            ssP << (ut->coefficient() < 0 ? "<": ">");
        } else if( coeff > 0) {

            ssP << (ut->coefficient() < 0 ? ">": "<");
        }
        ssP << " b[" << k  << "] )";
        ssP << "    return true;\n";
    }
    return ssP.str();
}


void CodeSynthesis::CreateSortIRComponent(UFCallTerm* currentUF,Computation* comp,
        int executionScheduleIndex) {
    // Add sort statement after
    iegenlib::Set* sortDomain = new iegenlib::Set("{[0]}");

    // Get execution schedule
    iegenlib::Relation* sortExecutionSchedule =
        getExecutionSchedule(
            sortDomain,executionScheduleIndex);
    TupleDecl tdl(currentUF->numArgs());
    std::pair<string,string> dataAccess = {currentUF->name(),"{"+tdl.toString(true)+
                                           "->[0]}"
                                          };
    comp->addStmt(new Stmt(currentUF->name()+ "->sort()",sortDomain->
                           prettyPrintString(),sortExecutionSchedule->
                           prettyPrintString(), {dataAccess}, {dataAccess}));
    delete sortExecutionSchedule;
    delete sortDomain;
}

// Create IR component generates an IR specification
// for an unknown.
void CodeSynthesis::CreateIRComponent(std::string currentUF,
                                      Computation* comp, int executionScheduleIndex,
                                      SynthExpressionCase ufCase, Exp* exp,
                                      std::vector<std::string>& unknowns,
                                      iegenlib::Set* transSet,
				      UFCallTerm* reorderUF) {

    std::string expStmt =
        constraintToStatement(exp,
                              currentUF,transSet->getTupleDecl(),ufCase);

    Set* ufDomain = GetCaseDomain(
                        currentUF,transSet,exp,ufCase);

    // Remove all self referential references to PERMUTE_NAME
    // from the ufDomain. This information is only important for
    // memory allocation. It looks like this code was already
    // called earlier, however projectOut introduces self referential
    // on Permutation everytime so we have to repeat the constraint
    // removal everytime a GetCaseDomain is called.

    for(auto selfRef : selfRefs) {
        // If self referential is on one of
        // the generated permutes, remove constraint
        auto it = std::find(permutes.begin(), permutes.end(),
                            selfRef.first);
        if (it != permutes.end()) {
            RemoveConstraint(ufDomain,selfRef.second);
        }
    }
    // remove constraints involving unknown UFs

    RemoveSymbolicConstraints(unknowns,ufDomain);
     

    // IF this space is uses a reorder function
    // change the iteration spaace to loop through
    // the reorder function
    if (reorderUF){
       Set* inv = GetInverseIterationSpace(ufDomain,reorderUF);
       delete ufDomain;
       ufDomain = inv; 
    }

    // Check if the domain created is valid for code generation,
    // if it isnt; throw exception and do not allow for codegeneration
    // for this set, it is better for the error to be handled here
    // than during late code generation phase.
    if (!IsValidIterationSpace(ufDomain)){
	std::cerr << "Bad: " << ufDomain->prettyPrintString() << "\n";
	throw assert_exception("CreateIRComponent:: Generated Domain"
			" is not feasible");
    }
    // Get reads and writes.
    auto ufWrites =
        GetWrites(currentUF,exp,ufCase,ufDomain->arity());


    auto ufReads =
        GetReads(currentUF,exp,ufCase,ufDomain->arity());

    // add data spaces for reads and writes to
    // IR
    addToDataSpace((*comp),
                   ufReads, "double");

    addToDataSpace((*comp),
                   ufWrites, "double");

    // Get execution schedule
    iegenlib::Relation* ufExecSched =
        getExecutionSchedule(
            ufDomain,executionScheduleIndex);

    comp->addStmt(new Stmt(expStmt,ufDomain->prettyPrintString(),
                           ufExecSched->prettyPrintString(),
                           ufReads,ufWrites));
    delete ufDomain;
    delete ufExecSched;
}

bool CodeSynthesis::IsDomainBoundedByUnknown(Term* term,
        const std::vector<std::string>& unknowns,iegenlib::Relation* rel ) {
    if (term->isUFCall()) {
        UFCallTerm* ufTerm = (UFCallTerm*) term;
        for(int i = 0; i < ufTerm->numArgs(); i++) {
            Exp* argExp = ufTerm->getParamExp(i);
            for(int j =0; j < rel->arity(); j++) {
                TupleVarTerm tuple(1,j);
                if (argExp->solveForFactor(tuple.clone())!=NULL) {
                    auto upperBounds = rel->GetUpperBounds(tuple);
                    auto lowerBounds = rel->GetLowerBounds(tuple);
                    for (auto unknown: unknowns) {
                        for(Exp* exp: upperBounds) {
                            if (findCallTerm(exp,unknown)!=NULL) {
                                return true;
                            }
                        }

                        for(Exp* exp: lowerBounds) {
                            if (findCallTerm(exp,unknown)!=NULL) {
                                return true;
                            }
                        }
                    }
                }
            }
        }
    } else if (term->type() == "TupleVarTerm") {
        TupleVarTerm* tupTerm = (TupleVarTerm*) term;
        auto upperBounds = rel->GetUpperBounds(*tupTerm);
        auto lowerBounds = rel->GetLowerBounds(*tupTerm);
        for (auto unknown: unknowns) {
            for(Exp* exp: upperBounds) {
                if (findCallTerm(exp,unknown)!=NULL) {
                    return true;
                }
            }

            for(Exp* exp: lowerBounds) {
                if (findCallTerm(exp,unknown)!=NULL) {
                    return true;
                }
            }
        }

    }
    return false;
}


// Class returns all aliases for a tuple
// variable in a sparse constraint.
// An alias is when some tv(x) = tv(y)
//

class FindTupleAliases: public Visitor{
private:
    int tv;
    std::vector<int> aliases;
public:
    FindTupleAliases(int tv): tv(tv){}
    void preVisitExp(Exp * e){
        if (e->isEquality()){
	    auto termList = e->getTermList();
	    for(Term* term:termList){
	        TupleVarTerm* t = dynamic_cast<TupleVarTerm*>(term);
		if (t == NULL || t->tvloc() != tv) continue;
		Term* tClone = t->clone();
		tClone->setCoefficient(1);
                Exp* rhs = e->solveForFactor(tClone);
		Term* termA= NULL;
		TupleVarTerm* tupTerm = NULL;
		if((termA = rhs->getTerm()) && 
				(tupTerm=dynamic_cast<TupleVarTerm*>(termA))){
		    aliases.push_back(tupTerm->tvloc());
		}

	    }
	}
    }
    std::vector<int> getAliases() { return aliases;}    
};

// Class finds tuple variable associated with 
// a uf. This class visits a set and returns 
// the expression of the rhs of the first UF
// it finds
class FindTupleVarTermVisitor: public Visitor{
private:
    string ufName;
    Exp* res = NULL;
    int arity;    
    int tvloc;
public:
    FindTupleVarTermVisitor(string ufName,int arity): arity(arity),
	    ufName(ufName){}
    ~FindTupleVarTermVisitor(){ delete res;}
    void preVisitExp(Exp* e){
        if (e->isEquality() && res == NULL){
	    Term* ufTerm = CodeSynthesis::findCallTerm(e,ufName);
	    if (ufTerm != NULL){
	        for(int i = 0; i < arity; i++){
	            TupleVarTerm tv(i);

		    if (e->dependsOn(tv)){
		       tvloc = i;
		       Term* tclone = ufTerm->clone();
		       tclone->setCoefficient(1);
		       res = e->solveForFactor(tclone);  
		    }
	        }
	    }
	} 
    }
    // Pointer is owned by the caller and must
    // be deallocated.
    Exp* getTupleExpression(){ 
	if (res == NULL) throw assert_exception(
	 		"FindTupleVarTermVisitor::getTupleExpression:"
			" expression not found");
	return res->clone();
    }
    int  getFoundTvLoc() { 
	
	if (res == NULL) throw assert_exception(
	 		"FindTupleVarTermVisitor::getTupleExpression:"
			" expression not found");
	return tvloc; 
    }  

};

Set* CodeSynthesis::GetInverseIterationSpace(Set* set, UFCallTerm* domUF){
    // remap tuple vars
    std::vector<int> remap(set->arity());
    // Shift all maps by 2.
    std::iota(remap.begin(), remap.end(),2);
    
    // Stores all iterations used so far
    uint32_t bitIter = 0;

    Set * s1 = new Set(set->arity()+2);
    Conjunction* conj = (*s1->conjunctionBegin());
    TupleDecl td = conj->getTupleDecl();
    for(int i = 0; i < set->arity(); i++){
        td.copyTupleElem(set->getTupleDecl(),i,i+2);
    }
    td.setTupleElem(0,"_no");
    td.setTupleElem(1,"_n");
    s1->setTupleDecl(td);


    // Add the iteration space for the reorder function.
    Exp* lb = new Exp();
    // _no >= 0
    lb->setInequality();
    lb->addTerm(new TupleVarTerm(0));
    bitIter |= (1<<0);
    conj->addInequality(lb);
   
    // _no < P1SIZE
    Exp* ub = new Exp();
    ub->setInequality();
    ub->addTerm(new VarTerm(1,domUF->name()+"SIZE"));
    ub->addTerm(new TupleVarTerm(-1,0));
    ub->addTerm(new Term(-1));
    conj->addInequality(ub);

    // _n = P1MAP(_no)
    UFCallTerm * pMap = new UFCallTerm(-1,domUF->name()+"MAP",1);
    Exp * pMapArg = new Exp();
    pMapArg->addTerm(new TupleVarTerm(0));
    pMap->setParamExp(0,pMapArg);
    Exp*  nOrEq = new Exp();
    nOrEq->addTerm(new TupleVarTerm(1));
    nOrEq->addTerm(pMap);
    nOrEq->setEquality();
    conj->addEquality(nOrEq);
    bitIter |= (1<<1);
    
    for(int i = 0; i < domUF->numArgs(); i++){
        Exp* arg = domUF->getParamExp(i);
        TupleVarTerm* t = dynamic_cast<TupleVarTerm*>(arg->getTerm());
        if (t == NULL) continue;
	TupleVarTerm* tClone = (TupleVarTerm*)t->clone();
        tClone->setCoefficient(1);
        UFCallTerm* p1Dim = new UFCallTerm(-1,
			domUF->name()+ "DIM" +std::to_string(i),1);
        tClone->remapLocation(remap);        
	bitIter|= (1<<tClone->tvloc());
	Exp* argP = new Exp();
	argP->addTerm(new TupleVarTerm(0));
        p1Dim->setParamExp(0,argP);
	Exp* p1DimEx = new Exp();
	p1DimEx->addTerm(tClone);
	p1DimEx->addTerm(p1Dim);
	p1DimEx->setEquality();
        conj->addEquality(p1DimEx);
    }
    //
    // k1 = n  new map
    // Find tuple variable associated with the reordering function
    // in the set.
    FindTupleVarTermVisitor ftv(domUF->name(),set->arity());
    set->acceptVisitor(&ftv);
    Exp* tupExp = ftv.getTupleExpression();
    tupExp->remapTupleVars(remap);
    tupExp->addTerm(new TupleVarTerm(-1,1));
    conj->addEquality(tupExp);
    bitIter|=(1<<(ftv.getFoundTvLoc()+2));
    int n = 0;
    // Find which tuple variable hasn't been 
    // used and assume tuple variable is 
    // original position no. First
    // check if that tuple variable has aliases.
    while(n < s1->arity()){
      if(!(bitIter & (1 << n))){
	   int tuplePos = n-2; // We shifted by 2 before
	                      // so we need to shift back to 
			      // get tuple position in the 
			      // original set.
           // Check if this tuple variable 
	   // is directly asigned to some other tuple
	   // variable. In the Set
	   
	   FindTupleAliases fta(tuplePos);
	   set->acceptVisitor(&fta);
	   std::vector<int> tupleAliases = fta.getAliases();
	   if (tupleAliases.size() > 0 ){
              // If it is add the aliases as part of the 
	      // constraint
	      for(auto tupAlias : tupleAliases){
	          
	          Exp* e = new Exp();
	          e->addTerm(new TupleVarTerm(-1,n));
	          e->addTerm(new TupleVarTerm(tupAlias+2));
	          e->setEquality();
	          conj->addEquality(e); 
	      }
	      n++;
	      continue;
	   }
	   // If there are no aliases 
	   // go ahead and make this the original
	   // position.
	   Exp* e = new Exp();
	   e->addTerm(new TupleVarTerm(-1,n));
	   e->addTerm(new TupleVarTerm(0));
	   e->setEquality();
	   conj->addEquality(e); 
       }
       n++;
    }


    return s1;
}
// Helper class that recursively visits
// a relation and finds output variables
// that are equal to or a function of 
// known UFs or input tuple variables.
class FindResolvedOutputTuple: public Visitor{
private:
    std::vector<int> res;
    int inputArity = -1;
    int outArity = -1;
    std::vector<std::string> unknownUFs;
    std::vector<int> resolvedTuples;
    void addToRes(int tvloc){
        auto it = std::find(res.begin(),res.end(), tvloc);
	if (it == res.end())
	    res.push_back(tvloc);
    }
public:
    FindResolvedOutputTuple(const std::vector<std::string>& unknownUFs,
		    const std::vector<int>& resolvedTuples):
	    unknownUFs(unknownUFs), resolvedTuples(resolvedTuples){
        for(int tup : resolvedTuples){
            addToRes(tup);
	}		
    }
    std::vector<int> getResolvableOutTuples(){
        return res;
    } 
    void preVisitRelation(Relation * rel){ 
        inputArity = rel->inArity();
	outArity   = rel->outArity();
    }

    void preVisitExp(Exp * e){
	if (inputArity == -1 || outArity == -1)
	    throw assert_exception("unknown outArity and"
			    " inArity: Was this run on a relation?");
        if (e->isEquality()){
	    auto termList = e->getTermList();
	    for(Term* term:termList){
	        TupleVarTerm* t = dynamic_cast<TupleVarTerm*>(term);
		if (t == NULL) continue;
                int tv_loc= t->tvloc();
                // We only care about tuple variable location
		// that is an output tuple.
		if (tv_loc < inputArity) continue;
		// Check if this has an unknown UF
	        bool hasUnknownUF = false;
	        for(auto unknown: unknownUFs){
		    if (CodeSynthesis::findCallTerm(e,unknown)!= NULL){
			    hasUnknownUF = true;
			    break;
			}
		    }
                    // Check if it depends on other output tuple
		    bool hasOtherOutputTuple =false;
		    for(int j = inputArity; j < inputArity + outArity; j++){
		        TupleVarTerm t(j);
		        if(j!=tv_loc  && e->dependsOn(t) && 
				std::find(resolvedTuples.begin(),
					resolvedTuples.end(),j) != 
	            			resolvedTuples.end()){
		          hasOtherOutputTuple = true;
			    break;
			}
		    }
		    //std::cout<< e->toString() << "\n";
		    //std::cout <<hasUnknownUF <<"," << hasOtherOutputTuple << "\n";
                     
		    if (!hasUnknownUF  && !hasOtherOutputTuple){
		        addToRes(tv_loc);
		    }

	    }	    

	}	
    }
};
std::vector<int> CodeSynthesis::GetResolvedOutputTuples(Relation* rel,
		const std::vector<std::string>& unknownUFs){
   // Fixed point solution exit once there 
   // are no changes in the resolvable tuple variable
   // results.
   int prevSize = 0;
   std::vector<int> res;
   while(true){
       FindResolvedOutputTuple frv(unknownUFs,res);
       rel->acceptVisitor(&frv);
       res = frv.getResolvableOutTuples();
       if (prevSize == res.size()) break;
       prevSize = res.size();
   }
   return res;
}

bool CodeSynthesis::IsValidIterationSpace(Set* set){
   // TODO: Check if set is satisfiable
	
   // ITerate through all the tuple variables in the set
   TupleDecl tup = set->getTupleDecl();
   for(int i = 0 ; i < tup.size(); i++ ){
      // Skip this tuple variable if it is a constant
      if (tup.elemIsConst(i)) continue; 
      TupleVarTerm tv(i);
      // Check if it has an upper and lower bound
      std::list<Exp*> upperBounds = set->GetUpperBounds(tv);
      std::list<Exp*> lowerBounds = set->GetLowerBounds(tv);
      // In a condition where a tuple variable only has 
      // one and not the other constraint then we return 
      // false because this is gonna fail. It means that the 
      // tuple has an upper bound and not or lower bound 
      // or has a lower bound and not an upper bound
      // which will fail code generation. 
      if(upperBounds.size() > 0 xor lowerBounds.size() > 0){
         //Deallocate content of the list 
	 for(Exp* e : upperBounds) delete e;
	 for(Exp* e : lowerBounds) delete e;
	 return false;
      }
      // This tuple variable has an upper bound and lower
      // bound.
      if(upperBounds.size() > 0 && lowerBounds.size() > 0){
         //Deallocate content of the list 
	 for(Exp* e : upperBounds) delete e;
	 for(Exp* e : lowerBounds) delete e;
         continue;
      }

      // Finally check if there is an equality expression
      // involving this tuple variable.
      bool foundExpression = false;
      for(auto conjIt = set->conjunctionBegin();
		      conjIt != set->conjunctionEnd();
		      conjIt++){
         for(Exp* e: (*conjIt)->equalities()){
	    for(Term* t : e->getTermList()){
	       TupleVarTerm* tvCast = dynamic_cast<TupleVarTerm*>(t);
	       if (tvCast!=nullptr &&
			       tvCast->tvloc() == i){
	          foundExpression = true;
		  break;
	       }
	    }
	 }
      }
      // If we do not find tuple expression 
      // involving tuple variable then 
      // this would noto generate an itertaion space
      if (!foundExpression) return false;
   }
   return true;
}
