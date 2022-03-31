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
using namespace code_synthesis;
using iegenlib::Exp;
/// Function flattens a sparse constraint : set, relation
/// to individual terms
/// \param sparseConstraint
/// \return
std::list<Term*> CodeSynthesis::getTermList(
      SparseConstraints * sparseConstraint){
  TermVisitor termVisitor;
  sparseConstraint->acceptVisitor(&termVisitor);
  return termVisitor.getTerms();
}






std::list<Term *> CodeSynthesis::getDependents(
		Conjunction *conjunction, Term *term) {

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

std::list <Term*> CodeSynthesis:: getParamTermList(
		const UFCallTerm *ufCallTerm) {
  std::list<Term*> ufTerms;
  for(int i = 0; i < ufCallTerm->numArgs();i++ ){
      for(auto t : ufCallTerm->getParamExp(i)->getTermList()){
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
 
    
     std::vector<std::string> unknowns;
    
     iegenlib::StringIterator* iter = destMapR->getSymbolIterator();
     while (iter->hasNext()) {
         std::string symb = iter->next();
    	// Exclude symbols that have already been specified to be 
    	// known.
    	if (std::find(knowns.begin(),knowns.end(),symb) == knowns.end()){
            unknowns.push_back( symb );
        }
     }
     int executionScheduleIndex  = 0;
     for (auto permute : permutes){
         // TODO: Mege permutes with CASE1 and CASE2
	 // CASE1, p0->insert({t1,t2})
	 // CASE2, p0(t1,t2) = t3
	 // Merged Permute 
	 // P0 = Permute(2);
	 // p0->insert({t1,t2,t3}) 
         for(auto e : expList){
             if (findCallTerm(e,permute)!=NULL){
                 auto caseP = 
     		    GetUFExpressionSynthCase(e,permute,
     		    transRel->inArity(),transRel->arity());
    	         if (caseP == SELF_REF) {
    	             selfRefs.push_back({permute,e});
    	         }
    	     } 
         } 
        
         
    
        // Remove all self referential references to permute 
        // from the transRel. This information is only important for 
        // memory allocation.
        for(auto selfRef : selfRefs){
            if (selfRef.first == permute){
    	        RemoveConstraint(transRel,selfRef.second);
    	    }
        }
        
	// Gather all Permute PExpCandidate
        std::vector<std::pair<Exp*,SynthExpressionCase>>
	       	permuteExpCandidate; 
        
        for(auto pExp : expList){
            if (findCallTerm(pExp,permute)!=NULL){
                auto caseP = GetUFExpressionSynthCase(pExp,permute,
     	        	transRel->inArity(),transRel->arity());
                if (caseP == CASE1 || caseP == CASE2)
		   permuteExpCandidate.push_back({pExp,caseP});    
         
	    }
	}
	// If there is just a single candidate. use 
	// the generic case insert.
        if (permuteExpCandidate.size() == 1){
	    auto caseP = permuteExpCandidate[0].second; 	
	    auto pExp = permuteExpCandidate[0].first;
	    std::string pStmt = 
    	           constraintToStatement(pExp,
    		   permute,composeRel->getTupleDecl(),
    		   caseP);
            // Get Domain for P
            iegenlib::Set* pDomain = 
    	        GetCaseDomain(permute,composeSet,pExp,caseP);
        
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
	}else {
            for(auto pCand : permuteExpCandidate) {
		if (pCand.second == CASE1)
		    continue;	
	        auto pExp = pCand.first;
		bool hasUnknown = false;
		for (auto uknwnUF : unknowns){
		    if(findCallTerm(pExp,uknwnUF)!=NULL){
		       hasUnknown = true;
		       break;
		    }
		}
		if (hasUnknown) continue;
	    auto permU =findCallTerm(pExp,permute);
            // add permute to mergedPermutes 
	    // so that Permutation will be initialized
	    // accordingly
            mergedPermutes.push_back({permute,permU->numArgs()});
	    
	    SynthExpressionCase caseP = MERGECASE;
	    std::string pStmt = 
    	           constraintToStatement(pExp,
    		   permute,composeRel->getTupleDecl(),
    		   caseP);
            
	    
	    // Get Domain for P
            iegenlib::Set* pDomain = 
    	        GetCaseDomain(permute,composeSet,pExp,caseP);
        
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
	
	    }
	}

     }
     
     std::vector<iegenlib::Set*> unknownDomain;
    
     for(std::string unknown: unknowns){
         iegenlib::Set* domain = transRel->GetDomain(unknown);
         unknownDomain.push_back(domain);	
     } 
    

    for (auto currentUF : unknowns){
        // skip UF for P since it has already been 
	// synthesized at this point.
        std::list<iegenlib::Exp*> expUfs;
	std::list<std::string> expStmts;
	for(auto e : expList){
	   if(findCallTerm(e,currentUF)!=NULL){


	     // Get case a uf in a constraint falls into.	   
	     auto ufCase = 
		     GetUFExpressionSynthCase(e,
		   currentUF,transRel->inArity(),transRel->arity());
	     if (ufCase ==SELF_REF){ 
	         // Store self referential constraints. Should help
	         // with generating code for montonicity later on.             
	        selfRefs.push_back({currentUF,e->clone()});
	        // skip this loop iteration
		continue;
	     }
	     // IF UF satisifies synthesis case
	     if(ufCase !=UNDEFINED && ufCase != CASE1){
	         std::string expStmt = 
	             constraintToStatement(e,
		       currentUF,transRel->getTupleDecl(),ufCase);
		 Set* ufDomain = GetCaseDomain(
				 currentUF,transSet,e,ufCase);
                 
                // Remove all self referential references to PERMUTE_NAME 
                // from the ufDomain. This information is only important for 
                // memory allocation. It looks like this code was already 
		// called earlier, however projectOut introduces self referential 
		// on Permutation everytime so we have to repeat the constraint 
		// removal everytime a GetCaseDomain is called.
                
		for(auto selfRef : selfRefs){
		    // If self referential is on one of 
		    // the generated permutes, remove constraint
		    auto it = std::find(permutes.begin(), permutes.end(),
				    selfRef.first);
                    if (it != permutes.end()){
	                RemoveConstraint(ufDomain,selfRef.second);
	            }
                }
                 // remove constraints involving unknown UFs
                 
	        RemoveSymbolicConstraints(unknowns,ufDomain);
                 
		 // Get reads and writes.
                 auto ufWrites = 
	                 GetWrites(currentUF,e,ufCase,ufDomain->arity()); 

                 
		 auto ufReads = 
	                 GetReads(currentUF,e,ufCase,ufDomain->arity()); 
                 
		 // add data spaces for reads and writes to 
		 // IR
	         addToDataSpace((*inspector),
				ufReads, "double");
                  
	         addToDataSpace((*inspector),
				ufWrites, "double");
                 
		 // Get execution schedule
                 iegenlib::Relation* ufExecSched = 
	                 getExecutionSchedule(
			    ufDomain,executionScheduleIndex++);

		 inspector->addStmt(new Stmt(expStmt,ufDomain->prettyPrintString(),
					 ufExecSched->prettyPrintString(),
					 ufReads,ufWrites));
                 delete ufDomain;
		 delete ufExecSched;
	     }
	   }
	}
	// Synthesize statements 
    }
    // Generate code to ensure universal constraint
    for(auto uf : unknowns){
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
    
    // CodeGen (RS2->S1(I)) - Data copy Code
    iegenlib::Set* copyDomain = composeRel->ToSet();
    std::string copyStmt = GetCopyStmt(sourceDataName,destDataName,destMapR,
		    sourceMapR);
    iegenlib::Relation* execSchedule = 
	    getExecutionSchedule(
			    copyDomain,executionScheduleIndex++);
    auto copyReads = getCopyReadAccess(sourceMapR,sourceDataName,copyDomain);
    auto copyWrites = getCopyWriteAccess(destMapR,destDataName,copyDomain);
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

UFCallTerm* CodeSynthesis::findCallTerm(Exp* exp, std::string ufName){
   UFCallTerm* res = nullptr;
   for(auto t : exp->getTermList()){
      if (t->isUFCall() && ((UFCallTerm*)t)->name() == ufName){
         res = (UFCallTerm*)t;
	 break;
      }
   }
   return res;
}

std::string CodeSynthesis::
constraintToStatement(Exp* constraint, 
		std::string unknownUF, const TupleDecl& 
		tupDecl, SynthExpressionCase expCase){
   UFCallTerm* ufTerm = findCallTerm(constraint,unknownUF);
   if (ufTerm == NULL){
      throw assert_exception("UFCallTerm must exist in expression");
   }
   int tupleSize = tupDecl.size();
   // Solve for UF term.
   Term* ufClone = ufTerm->clone();
   ufClone->setCoefficient(1);
   Exp* solvedUFConst = constraint->solveForFactor(ufClone);
   
   std::stringstream ss;
   if(expCase == CASE1){
      ss << unknownUF << "->insert({";
      bool firstArg = true;
      for (int i = 0;i <ufTerm->numArgs(); ++i) {
          if (not firstArg) { ss << ", "; }
          if (ufTerm->getParamExp(i)) { 
             ss << ufTerm->getParamExp(i)->prettyPrintString(tupDecl); }
             firstArg = false;
       }
       ss << "})";
    }else if (expCase == CASE2){
       ss << ufTerm->prettyPrintString(tupDecl,true) << "=" 
		    << solvedUFConst->prettyPrintString(tupDecl);
       
    }else if (expCase == CASE3){
	    ss << ufTerm->prettyPrintString(tupDecl,true) << " = " 
	    << "min(" << ufTerm->prettyPrintString(tupDecl,true)<< "," 
	    << solvedUFConst->prettyPrintString(tupDecl)
	    << ")";
       
    }else if (expCase == CASE4){
	
        ss << ufTerm->prettyPrintString(tupDecl,true) << " = " 
		    << "max(" << ufTerm->prettyPrintString(tupDecl,true)<< ","
		    <<solvedUFConst->prettyPrintString(tupDecl)
		    << ")";
    }else if(expCase == MERGECASE){
      ss << unknownUF << "->insert({";
      bool firstArg = true;
      for (int i = 0;i <ufTerm->numArgs(); ++i) {
          if (not firstArg) { ss << ", "; }
          if (ufTerm->getParamExp(i)) { 
             ss << ufTerm->getParamExp(i)->prettyPrintString(tupDecl); }
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
  for(auto l : list){
      if(first){
          ss << l;
          first = false;
      }else{
          ss << "," <<l;
      }

  }
  ss << "]";
  return ss.str();
}

iegenlib::Set* CodeSynthesis::GetCaseDomain(std::string ufName,Set* s,
		Exp * constraint, SynthExpressionCase expCase){
    
    Set * res;
    if (expCase == CASE1){
        res = s->GetDomain(ufName);
    }else {
	// Get the maximum tuple variable 
	// present in this expression.
	
	res = new Set(*s);
	int maxTup = -1;
        for(int i = 0; i < s->arity(); i++){
	    TupleVarTerm t(1,i);
            if(constraint->dependsOn(t)){
	       maxTup =  std::max(i,maxTup);
	    }
	}

	if (maxTup == -1){
		throw assert_exception("GetCaseDomain: no domain"
				" available for expression");
        }
	// Project out tuple variables after 
	// maxTuple.
        while(res->arity() - 1> maxTup ){
	    Set * temp = res->projectOut(res->arity() - 1);
	    delete res;
	    res = temp;
	}
    }
    return res;
}

CodeSynthesis::CodeSynthesis(SparseFormat* source,
		SparseFormat* dest){
    destMapR = new Relation(dest->mapToDense);
    sourceMapR = new Relation(source->mapToDense);
    if(sourceMapR->outArity()!= destMapR->outArity()){
        throw assert_exception("CodeSynthesis:: Format Descriptor map must"
			" have the same output arity");
    }
    auto invDestMap = destMapR->Inverse();  
    // Add the Permutation constraint.
    permutes = AddPermutationConstraint(invDestMap);
    
    composeRel = invDestMap->Compose(sourceMapR);
    
    std::cerr << "Compose Rel: "<< composeRel->prettyPrintString() 
	       << " \n";
    
       
    
    transRel = composeRel->TransitiveClosure();
    
    std::cerr << "Trans Rel: "<< transRel->prettyPrintString()
	    << "\n";

    transRelExpanded = substituteDirectEqualities(transRel);
    

    sourceDataName = source->dataName;
    destDataName = dest->dataName;

    sourceDataConstraint = source->dataConstraint;
    destDataConstraint = source->dataConstraint;
    
    sourceDataAccessMap = new Relation(source->dataAccess);
    destDataAccessMap = new Relation(dest->dataAccess);    
    


    // Setup UFquantifiers
    iegenlib::setCurrEnv();
    for(auto uf: source->ufQuants){
       Set* ufDomain = new Set(uf.domain);
       Set* ufRange = new Set(uf.range);
       iegenlib::appendCurrEnv(uf.name,ufDomain,ufRange,false,
		       uf.type);
    }

    for(auto uf: dest->ufQuants){
       Set* ufDomain = new Set(uf.domain);
       Set* ufRange = new Set(uf.range);
       iegenlib::appendCurrEnv(uf.name,ufDomain,ufRange,false,
		       uf.type);
       ufQuants.push_back(uf);
    }
    
    for(auto known : dest->knowns){
       knowns.push_back(known);
    }

    for(auto known : source->knowns){
       knowns.push_back(known);
    }
}


void CodeSynthesis::RemoveConstraint(SparseConstraints* sc, Exp *e ){
        Conjunction* conj  = *sc->conjunctionBegin();
        if( e->isEquality()){
	
	auto it = std::find_if(conj->equalities().begin(),
			conj->equalities().end(),
			[e](Exp* e1){
			     return e1->toString() == (e)->toString();}
			);
	if (it != conj->equalities().end()){
	    conj->equalities().erase(it);
	}
	}else if (e->isInequality()){
	
	auto it = std::find_if(conj->inequalities().begin(),
			conj->inequalities().end(),
			[e](Exp* e1){
			     return e1->toString() == (e)->toString();}
			);
	if (it != conj->inequalities().end()){
	    conj->inequalities().erase(it);
	}
	}
}

std::vector<std::string> CodeSynthesis::
    AddPermutationConstraint(Relation* rel){
    // TODO: we can't add permutes whose 
    // right hand side equals to some function 
    // in the input tuple.
    std::vector<std::string> permutes;
    for(int i = 0; i < rel->outArity(); i++){
        Exp* e  = new Exp();
	std::string name = PERMUTE_NAME + std::to_string(i);
        UFCallTerm* pUF = new UFCallTerm(1,name,rel->inArity());
        for(int i =0 ; i < rel->inArity(); i++){
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
		    ++it){
            (*it)->addEquality(e);
        }
	permutes.push_back(name);
    }
    return permutes;
}

iegenlib::Relation* CodeSynthesis::solveForOutputTuple(iegenlib::Relation* r){
    assert(r->outArity()==1 && "Output arity must be 1"); 
    // Get Expressions involving output tuple variables
    TupleVarTerm * term= new TupleVarTerm(r->inArity());
    ExpTermVisitor expVisit(term);
    r->acceptVisitor(&expVisit);
    auto exps = expVisit.getExpressions(); 
    
    for(auto &e : exps){
        auto solution = e->solveForFactor(term);
	if( not solution){
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
	       	int pos){
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

    for(int i = 0; i < execTupSize; i++){
	Exp* e = new Exp();
	e->setEquality();
	TupleVarTerm* t2 = new TupleVarTerm(1,size+i+1);
	e->addTerm(t2);
        if ( i %2 == 0){
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
class DataAccessVisitor: public Visitor{
private:
    std::vector<std::pair<std::string,std::string>> dAccess;
    int arity;
public:
    DataAccessVisitor(unsigned int arity): arity(arity){}
    void preVisitUFCallTerm(UFCallTerm * t);
    void preVisitVarTerm(VarTerm* t);
    std::vector<std::pair<std::string,std::string>> getDataAccess(){ 
	    return dAccess;}
    void clearDataAccesses(){dAccess.clear();}
};

void DataAccessVisitor::preVisitVarTerm(VarTerm* t){
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
   std::string relString = rel->prettyPrintString();
   delete rel;
   dAccess.push_back({dataName,relString});
}
void DataAccessVisitor::preVisitUFCallTerm(UFCallTerm* t){
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
    for(int i = 0; i < t->numArgs(); i++){
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
     SynthExpressionCase expCase,int arity){
    std::vector<std::pair<std::string,std::string>> result;
    UFCallTerm* ufTerm = findCallTerm(constraint,uf);
    if (ufTerm == NULL){
       throw assert_exception("UFCallTerm must exist in expression");
    }
    // Solve for UF term.
    Term* ufClone = ufTerm->clone();
    ufClone->setCoefficient(1);
    Exp* solvedUFConst = constraint->solveForFactor(ufClone);
    if (expCase == CASE1 || expCase == MERGECASE){
        // This is for case 1 where 
	// we have an insert abstraction
	TupleDecl tdl(arity);
	result.push_back({uf,"{"+tdl.toString(true)+
            "->[0]}"});
    }else if (expCase ==  CASE2 ||
		   expCase == CASE3 ||
		   expCase == CASE4|| expCase == SELF_REF ){
	DataAccessVisitor dV(arity);
	ufTerm->acceptVisitor(&dV);
	// We know there is only one UF write,
	// so we only pick up data accesses 
	// that involves UF on LHS
	for(auto dAccess: dV.getDataAccess()){
	    if (dAccess.first == uf){
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
	int arity){
    UFCallTerm* ufTerm = findCallTerm(constraint,uf);
    if (ufTerm == NULL){
       throw assert_exception("UFCallTerm must exist in expression");
    }
    
    // Solve for UF term.
    Term* ufClone = ufTerm->clone();
    ufClone->setCoefficient(1);
    Exp* solvedUFConst = constraint->solveForFactor(ufClone);
    DataAccessVisitor dV(arity);
    solvedUFConst->acceptVisitor(&dV);
    auto result = dV.getDataAccess();
    dV.clearDataAccesses();
    // There might be some read accesses inside of the
    // ufTerm. 
    ufTerm->acceptVisitor(&dV);
    auto ufTermDataAccesses = dV.getDataAccess();
    for(auto dataAccess : ufTermDataAccesses){
        if(dataAccess.first != uf){ 
	    result.push_back(dataAccess);
	}
    }
    return result;
}



/// Function returns case of expression as regards a UF. 
SynthExpressionCase CodeSynthesis::GetUFExpressionSynthCase(Exp* constraint,
      std::string unknownUF, int inputArity, int tupleSize){
   SynthExpressionCase caseResult = UNDEFINED;
   UFCallTerm* ufTerm = findCallTerm(constraint,unknownUF);
   if (ufTerm == NULL){
      throw assert_exception("UFCallTerm must exist in expression");
   }
   // Solve for UF term.
   Term* ufClone = ufTerm->clone();
   ufClone->setCoefficient(1);
   Exp* solvedUFConst = constraint->solveForFactor(ufClone);
   
   if (solvedUFConst == NULL)
	 return UNDEFINED;

   // if UF is self referential 
   if (findCallTerm(solvedUFConst,unknownUF)!= NULL){
       return SELF_REF;
   }
   std::stringstream ss;
   if (constraint->isEquality()){
      //Case 1
      // If rhs only has one term and the term is an output tuple
      // var term.
      if (solvedUFConst->getTermList().size()== 1
   		   && solvedUFConst->getTerm()->type() == "TupleVarTerm"
		   && ((TupleVarTerm*)solvedUFConst->getTerm())->
		      tvloc() >= inputArity){
          caseResult = CASE1;
      }else{
         //Case 2
	 //UF(x) = F(x)
	 //solved for must not depend on output term,
	 //Will need the number of tuple declarations here.
         bool dependsOnOutput  = false;
	 for(int i = inputArity; i < tupleSize; i++){
            TupleVarTerm t(1,i);
	    if (solvedUFConst->dependsOn(t)){
	       dependsOnOutput = true;
	       break;
	    }
	 }
	 if (not dependsOnOutput){
		 
             caseResult = CASE2;
           
	 }
      }
   
   }else {
      // All cases in this section 
      // arity(y) > arity(x)
      int x_arity = 0;
      int y_arity = 0;
      bool x_dependsOnTuple = false;
      bool y_dependsOnTuple = false;
      for(int i = 0 ; i < inputArity; i++){
         TupleVarTerm t(1,i);
         for(int k = 0; k < ufTerm->numArgs(); k++){
	    if(ufTerm->getParamExp(k)->
			    dependsOn(t)){
	       x_arity++;
	       x_dependsOnTuple = true;
	       break; 
	    }
	 }
         if(solvedUFConst->dependsOn(t)){
	    y_arity++;
	    y_dependsOnTuple = true;
	 }
      }
      // Relax y_arity >= x_arity constraints
      if (y_arity >= x_arity && x_dependsOnTuple && y_dependsOnTuple){
         // Case 3
         // UF(x) <= F(y) 
         if (ufTerm->coefficient() < 0){
	     caseResult = CASE3; 
	 }else if (ufTerm->coefficient() > 0){
	     caseResult = CASE4;
	 }
      }
      
   }
   delete solvedUFConst;
   return caseResult;
}	

std::string CodeSynthesis::getSupportingMacros (){
   std::stringstream ss;
   ss << "#define min(a,b) a < b ? a : b\n"
	   << "#define max(a,b) a > b ? a: b\n";
   return ss.str();
}

void CodeSynthesis::addToDataSpace(Computation& comp, 
		      std::vector<std::pair<std::string,std::string>> access,std::string baseType){
    for(auto a : access){
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

void CodeSynthesis::RemoveSymbolicConstraints(const std::vector<std::string>& symbNames,
		     SparseConstraints* sc){
   for(auto it = sc->conjunctionBegin();
                   it!= sc->conjunctionEnd();
                   ++it){
      Conjunction* conj = (*it);
      auto itE = conj->equalities().begin();
        while(itE != conj->equalities().end()){
	   bool found = false;   
	   for(auto t: (*itE)->getTermList()){
		 std::string name;
	         if(t->isUFCall()){
	            UFCallTerm* ut = dynamic_cast<UFCallTerm*>(t);
                    name = ut->name();
		 }else{
	            VarTerm* vt = dynamic_cast<VarTerm*>(t);
		    name = vt? vt->symbol():"";
		 }
		 auto itU = 
			std::find(symbNames.begin(),symbNames.end(),name);
                 if(itU != symbNames.end()){
		    delete (*itE);
                    itE = conj->equalities().erase(itE);
		     found = true;
		     break;
		 }
	      }
	   if (!found){ itE++;}
      }

      itE = conj->inequalities().begin();
      while(itE != conj->inequalities().end()){
	   bool found = false;   
	   for(auto t: (*itE)->getTermList()){
		 std::string name;
	         if(t->isUFCall()){
	            UFCallTerm* ut = dynamic_cast<UFCallTerm*>(t);
                    name = ut->name();
		 }else{
	            VarTerm* vt = dynamic_cast<VarTerm*>(t);
		    name = vt? vt->symbol():"";
		 }
		 auto itU = 
			std::find(symbNames.begin(),symbNames.end(),name);
                 if(itU != symbNames.end()){
		     delete (*itE);
                     itE = conj->equalities().erase(itE);
		     found = true;
		     break;
		 }
	      }
	   if (!found){ itE++;}
       }
   }
   

}

// Function returns read accesses for code generated 
// in a montonic statement.
std::vector<std::pair<std::string,std::string>> CodeSynthesis::
getMonotonicReadAccess(std::string uf,MonotonicType type){
    if (type == Monotonic_NONE){
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
		     MonotonicType type){
    if (type == Monotonic_NONE){
        throw assert_exception("getMonotonicStmt: none monotonic type"
			" is not supported.");
    }
    // Only writes to uf(e2) in all cases 
    std::vector<std::pair<std::string,std::string>>res;
    res.push_back({uf,"{[e1,e2] -> [e2]}"});
    return res;
}



std::string CodeSynthesis::getMonotonicStmt(std::string uf,MonotonicType type){
    if (type == Monotonic_NONE){
        throw assert_exception("getMonotonicStmt: none monotonic type"
			" is not supported.");
    }
    std::stringstream ss;
    ss <<  "if (";
    if (type == Monotonic_Increasing){
	// if ( not (uf (e1) < uf(e2)) ){
        ss << " not (" << uf << "(e1) < "<< uf << "(e2))){";
        //   UF(e2) = UF(e1) + 1
        ss << uf << "(e2) = " << uf << "(e1) + 1;";	
    }else  
    if (type == Monotonic_Nondecreasing){
	// if ( not (uf (e1) <= uf(e2)) ){
        ss << " not (" << uf << "(e1) <= "<< uf << "(e2))){";
        //   UF(e2) = UF(e1)
        ss << uf << "(e2) = " << uf << "(e1);";	
    }else 
    if (type == Monotonic_Decreasing){
  	// if ( not (uf (e1) > uf(e2)) ){
        ss << " not (" << uf << "(e1) > "<< uf << "(e2))){";
        //   UF(e2) = UF(e1) - 1
        ss << uf << "(e2) = " << uf << "(e1) - 1;";	
     }else 
     if (type ==Monotonic_Nonincreasing){
  	// if ( not (uf (e1) >= uf(e2)) ){
        ss << " not (" << uf << "(e1) >= "<< uf << "(e2))){";
        //   UF(e2) = UF(e1)
        ss << uf << "(e2) = " << uf << "(e1) - 1;";	
     }
     ss << "}";
    return ss.str();
}




Set* CodeSynthesis::GetMonotonicDomain(std::string uf, MonotonicType type,
       Set* ufDomain){
    if (ufDomain->arity()!=1){
        throw assert_exception("GetMonotonicDomain: Uf domain"
			" must have an arity of 1");
    }
    if (type == Monotonic_NONE){
        throw assert_exception("GetMonotonicDomain: none monotonic type"
			" is not supported.");
    }
    std::stringstream ss; 
    TupleVarTerm tVar(0);

    std::list<Exp*> upperBounds = ufDomain->GetUpperBounds(tVar);
    std::list<Exp*> lowerBounds = ufDomain->GetLowerBounds(tVar);

    ss << "{[e1,e2]: e1 < e2 ";
    for(Exp* e : upperBounds){
        ss << " && e1 <= " << e->toString();
	ss << " && e2 <= " << e->toString(); 
        // DELETE e Since we own e. I do not think this is a good idea
	// I suggest we refactor such that we don't own the 
	// expression. TODO
	delete e; 
    }


    for(Exp* e : lowerBounds){
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



Exp* CodeSynthesis::getMonotonicDiff(std::string uf,Exp* ex){
   // TODO: revisit this.
   UFCallTerm* ufTerm = findCallTerm(ex,uf);
   if (ufTerm == NULL){
      throw assert_exception("getMonotonicDiff:"
		      " UFCallTerm must exist in expression");
   }
   // Solve for UF term.
   Term* ufClone = ufTerm->clone();
   ufClone->setCoefficient(1);
   Exp* solvedUFConst = ex->solveForFactor(ufClone);
   UFCallTerm* ufTerm2 = findCallTerm(solvedUFConst,uf);

   if (ufTerm2 == NULL){
      throw assert_exception("getMonotonicDiff: Expression"
		      " is not self referential");
   }
   
   if(ufTerm2->numArgs() !=  1 || ufTerm->numArgs() != 1){
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
	      Relation* destMap, Relation* sourceMap){
   std::stringstream ss;
   ss << destDataName << "(";
   bool isFirst =true;
   for(int i =0; i < destMap->inArity(); i++){ 
      if(isFirst){ 
         ss << destMap->getTupleDecl().elemVarString(i);
         isFirst = false;
      }else{
         ss << "," <<  destMap->getTupleDecl().elemVarString(i);
      }
   }
   ss << ") = ";
   isFirst = true;
   ss << sourceDataName << "(";
   for(int i =0; i < sourceMap->inArity(); i++){ 
      if(isFirst){ 
         ss << sourceMap->getTupleDecl().elemVarString(i);
         isFirst = false;
      }else{
         ss << "," <<  sourceMap->getTupleDecl().elemVarString(i);
      }
   }
   ss << " )";
   return ss.str();
}

std::vector<std::pair<std::string,std::string>> 
      CodeSynthesis::getCopyWriteAccess(Relation* destMap, std::string destDataName, 
		      Set* domain){
    std::stringstream ss;
    bool isFirst = true;
    ss <<"[";
    for(int i =0; i < destMap->inArity(); i++){ 
       if(isFirst){ 
          ss << destMap->getTupleDecl().elemVarString(i);
          isFirst = false;
       }else{
          ss << "," <<  destMap->getTupleDecl().elemVarString(i);
       }
    }
    ss << "]";
    std::vector<std::pair<std::string,std::string>>res;
    res.push_back({destDataName,"{"+domain->getTupleDecl().
		    toString(true)+" -> "+ss.str()+"}"});
    return res;
}


std::vector<std::pair<std::string,std::string>> 
	      CodeSynthesis::getCopyReadAccess(Relation* sourceMap,
			      std::string sourceDataName, 
			      Set* domain){
    bool isFirst = true;
    std::stringstream ss;
    ss <<"[";
    for(int i =0; i < sourceMap->inArity(); i++){ 
       if(isFirst){ 
          ss << sourceMap->getTupleDecl().elemVarString(i);
          isFirst = false;
       }else{
          ss << "," <<  sourceMap->getTupleDecl().elemVarString(i);
       }
    }
    ss << "]";
    std::vector<std::pair<std::string,std::string>>res;
    res.push_back({sourceDataName,"{"+domain->getTupleDecl().
		    toString(true)+" -> "+ss.str()+"}"});
    return res;
}


std::string CodeSynthesis::generateFullCode(){
    Computation* comp = generateInspectorComputation();
    std::cout << "=======IR======\n";
    comp->printInfo();
    std::cout << "=======IR-END====\n";
    std::stringstream ss;
    ss << getSupportingMacros();
    for(auto permute : permutes ){
        auto mergeIT = std::find_if(mergedPermutes.begin(),
		       mergedPermutes.end(), 
			[&permute](std::pair<std::string,int>& val){
				return val.first == permute;
			});
        if (mergeIT != mergedPermutes.end()){
           ss << "Permutation<int> * " << permute <<
	        " = new Permutation<int>("<<mergeIT->second <<");\n";
	   continue;
	}	
	std::string permInit = "Permutation<int> * "+permute +
	        " = new Permutation<int>();\n";
	auto it = std::find_if(selfRefs.begin(),
			selfRefs.end(),
			[&permute](std::pair<std::string,iegenlib::Exp*>& val){
				return val.first == permute;
			});
	if (it != selfRefs.end()){
	    // This has to be ordered by how much permute
	    // level is different in self referential
            Exp * e = it->second;
	    UFCallTerm* ut = findCallTerm(e,permute);
            Term* cloneUT = ut->clone();
	    cloneUT->setCoefficient(1);
 	    Exp* solveE = e->solveForFactor(cloneUT);
	    UFCallTerm* ut2 = findCallTerm(solveE,permute);;
	    if (ut2 == NULL) continue;
            std::stringstream ssP;
	    ssP << "Permutation<int>* "<< permute
		    << " = new Permutation <int>([]("
		   << " std::vector<int>& a, std::vector<int>& b){\n"; 
            for(int k = 0; k < ut->numArgs(); k++){
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
	       if (t!=NULL){
                   coeff = t->coefficient();
               }
               if (coeff ==0) continue;
	       ssP << "if (a[" << k  << "] ";
	       if (coeff < 0 ){
	          ssP << (ut->coefficient() < 0 ? "<": ">");   	  
	       }else if( coeff > 0){
	       
	          ssP << (ut->coefficient() < 0 ? ">": "<");   	  
	       }
	       ssP << " b[" << k  << "] )";
	       ssP << "    return true;\n";
	    }

	    ssP << "return false;\n";
	    ssP << "});\n";
	    permInit = ssP.str();
	}else {
	    std::string comp = GetPermuteComparator(permute,
			    composeRel,ufQuants);
            permInit = "Permutation<int> * "+permute+
		    " = new Permutation <int>("+comp+");\n";
	}
	ss <<  permInit;
    }
    // Add Datamacros for Source and Destination
    // #define <sourceDataName>(i) <sourceDataName>[i]
    // #define <destDataName>(i) <destDataName>[
    bool isFirst = true;
    std::string srcD1 = sourceDataName+ "(";
    for(int i =0; i < sourceMapR->inArity(); i++){ 
       if(isFirst){ 
          srcD1+= sourceMapR->getTupleDecl().elemVarString(i);
          isFirst = false;
       }else{
          srcD1 += ","  + sourceMapR->getTupleDecl().elemVarString(i);
       }
    }
    srcD1 += ")";
    
    std::string srcD2 = sourceDataName;
    for(int i=0; i < sourceDataAccessMap->outArity(); i++){
       srcD2 += "[" + sourceDataAccessMap->getTupleDecl().
	       elemVarString(i+sourceDataAccessMap->inArity()) + "]";
    }

    ss << "#define " << srcD1 << " "<< srcD2 << "\n";
    isFirst = true;
    std::string destD1 = destDataName+ "(";
    for(int i =0; i < destMapR->inArity(); i++){ 
       if(isFirst){ 
          destD1+= destMapR->getTupleDecl().elemVarString(i);
          isFirst = false;
       }else{
          destD1 += ","  + destMapR->getTupleDecl().elemVarString(i);
       }
    }
    destD1 += ")";
    
    std::string destD2 = destDataName;
    for(int i=0; i < destDataAccessMap->outArity(); i++){
       destD2 += "[" + destDataAccessMap->getTupleDecl().
	       elemVarString(i+destDataAccessMap->inArity()) + "]";
    }
    ss << "#define " << destD1 << " "<< destD2 << "\n";
    std::string code = comp->codeGen();
    
    for(auto permute : permutes){
        // Replace P[][] access to P->get({});
        // #define P(i,j) P[i][j] becomes:
        // #define P(i,j) p->get({i,j})
        std::string p1 = permute+ "(";
        std::string p2 = permute;
        std::string p3 = permute+ "->get({";
        isFirst = true;
        for (int i =0 ; i < destMapR->outArity(); i++){
            if(isFirst){
	        p1+="t" + std::to_string(i);
	        p3+="t" + std::to_string(i);
	        isFirst = false;
            }else{
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
    }
    ss <<code;
    delete comp;
    return ss.str();    
}

std::string CodeSynthesis::GetSupportHeader(){
   std::stringstream ss;
   ss << "#ifndef SYNTH_HEADER\n";
   ss << "#define SYNTH_HEADER\n";
   ss << "#include <functional>\n";
   ss << "#include <algorithm>\n";
   ss << "#include <vector>\n";
   ss << "#include <assert.h>\n";
   ss << "#include <iostream>\n";
   ss << "#include <string>\n";
   ss << "#include <sstream>\n";
   ss << "// Define P Data structure\n";
   ss << "template <typename T>\n";
   ss << "using Comparator = std::function<bool (std::vector<T>&,std::vector<T>&)>;\n";
   ss << "\n";
   ss << "template <typename T>\n";
   ss << "class Permutation{\n";
   ss << "private:\n";
   ss << "    std::vector<std::vector<T>> d;\n";
   ss << "    int tupleSplit = 0;\n";
   ss << "    Comparator<T> sortConstraint;\n";
   ss << "public:\n";
   ss << "    Permutation(Comparator<T> sortConstraint): tupleSplit(tupleSplit),\n";
   ss << "	sortConstraint(sortConstraint){}\n";
   ss << "    Permutation(){\n";
   ss << "       this->sortConstraint = NULL;\n";
   ss << "    }\n";
   ss << "    Permutation(int tupleSplit): tupleSplit(tupleSplit) {}\n";
   ss << "    void insert(std::vector<T> tup){\n";
   ss << "        d.push_back(tup);\n";
   ss << "	if (sortConstraint != NULL){\n";
   ss << "	    std::sort(d.begin(),d.end(),sortConstraint);\n";
   ss << "	}\n";
   ss << "    }\n";
   ss << "    int get(std::vector<T> tup){\n";
   ss << "        typename std::vector<std::vector<T>>::iterator it;\n";
   ss << "    	if (tupleSplit == 0){\n";
   ss << "	    it = std::find(d.begin(),d.end(),tup);\n";
   ss << "	}else{\n";
   ss << "	    it = std::find_if(d.begin(),d.end(),[this,&tup](std::vector<T> &a){\n";
   ss << "			        for(int i=0; i < tupleSplit; i++){\n";
   ss << "				    if(a[i] != tup[i]) return false;\n";
   ss << "				}\n";
   ss << "				return true;\n";
   ss << "			    });\n";
   ss << "	}\n";
   ss << "	if (it == d.end()) {\n";
   ss << "	    std::stringstream ss;\n";
   ss << "	    ss << \"Permutation::get: Tuple {\";\n";
   ss << "\n";
   ss << "	    for(int j = 0; j  < tup.size(); j++){\n";
   ss << "	        ss << tup[j] << ",";\n";
   ss << "	    }\n";
   ss << "	    ss << \"} not found\";\n";
   ss << "	    std::cerr << ss.str();\n";
   ss << "	    assert(0 && ss.str().c_str());\n";
   ss << "	}\n";
   ss << "	if (tupleSplit == 0) return it - d.begin();\n";
   ss << "	else return (*it)[tupleSplit];\n";
   ss << "    }\n";
   ss << "    std::string toString(){\n";
   ss << "	std::stringstream ss;\n";
   ss << "	for(int i = 0; i < d.size(); i++){\n";
   ss << "	    ss<< \"[\" << i << \"] => {\";\n";
   ss << "	    for(int j = 0; j  < d[i].size(); j++){\n";
   ss << "	        ss << d[i][j] << \",\";\n";
   ss << "	    }\n";
   ss << "	    ss << \"}\";\n";
   ss << "	}\n";
   ss << "	return ss.str();\n";
   ss << "    }\n";
   ss << "};\n";
   ss << "\n";
   ss << "#endif\n";
   return ss.str(); 
}



// Function checks if an expression depends on output 
// tuple returns true if it does.
// \param arity   arity of the constraint
// \param inArity input arity of the constraint
// \param e       expression
// Params are not adopted
bool CodeSynthesis::dependsOnOutputTuple(int arity, int inArity,
		      iegenlib::Exp*e){

    for(int i = inArity; i < arity; i++){
           TupleVarTerm tup(i);
           if (e->dependsOn(tup)){
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
Relation* CodeSynthesis::substituteDirectEqualities(Relation* rel){
   Relation* res = new Relation(*rel);
   SubMap subMap;
   for(auto conj = res->conjunctionBegin();
		conj!= res->conjunctionEnd(); conj++){
        auto itE = (*conj)->equalities().begin();
        while(itE != (*conj)->equalities().end()){
	   auto eq = (*itE);
	   if(!eq) continue;
           for(int i =res->inArity(); i < res->arity(); i++){
               TupleVarTerm tup(i);
               if (eq->dependsOn(tup)){
		   Term* t = tup.clone();
	           Exp* e = eq->solveForFactor(t);
		   if (e && !dependsOnOutputTuple(res->arity(),
					   res->inArity(),
					   e)){
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
		      std::vector<std::string>& unknowns){
    return false;
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

std::string CodeSynthesis::GetPermuteComparator(std::string& permute,Relation* composeRel,std::vector<UFQuant>& ufQuants){
    std::stringstream ss;
    ss << "[](std::vector<int>& a,std::vector<int>& b){\n";
    // Discover the range of the permute in 
    // the relation. If equivalent to one of the UFs
    // enforce constraints from that UF 
    Conjunction* c = *composeRel->conjunctionBegin();
    for(auto e : c->equalities()){
        UFCallTerm *permTerm = NULL;
	if (permTerm = findCallTerm(e,permute)){
	  
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
	  for(Exp* e : upperBounds){
	     Term* tupClone = tup.clone();
	     tupClone->setCoefficient(-1);
	     e->addTerm(tupClone);
	     e->setInequality();
	     permRangeConj->addInequality(e);
	  }

	  for(Exp* e : lowerBounds){
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
	  for (auto ufQuant : ufQuants){
	     Set* ufQuantDomain = new Set(ufQuant.domain);
	     if (*ufQuantDomain == *permuteRange && 
			     ufQuant.rhsProperty!= ""){
	        Set* rhsProperty = new Set(ufQuant.rhsProperty);
		auto rhsConj = *rhsProperty->conjunctionBegin();
                // Now we create comparator constraints based on this.
                ss << "if(";
		bool first = true;
		for (Exp* rhsExp : rhsConj->inequalities()){
		   if (first) { first = false;}
		   else ss << " && "; 
		   ss << rhsExp->prettyPrintString(rhsConj->
				   getTupleDecl());
                   ss << ">= 0";   
		}
		for (Exp* rhsExp : rhsConj->equalities()){
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
    ss << "    return false;\n";
    ss << "}";
    return ss.str();
}

