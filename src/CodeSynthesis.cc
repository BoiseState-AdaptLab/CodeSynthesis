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
/// TODO: Change to use TermVisitor
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
  delete mapToNewSpace;
  delete originalSpace;
}

Computation* CodeSynthesis::generateInspectorComputation() {
     // Code is currently been written in code_synthesis_main.cpp
     return NULL;
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
      ss << unknownUF << ".insert(";
      bool firstArg = true;
      for (int i = 0;i <ufTerm->numArgs(); ++i) {
          if (not firstArg) { ss << ", "; }
          if (ufTerm->getParamExp(i)) { 
             ss << ufTerm->getParamExp(i)->prettyPrintString(tupDecl); }
             firstArg = false;
       }
       ss << ")";
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
    if (expCase == CASE1){
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
std::vector<std::pair<std::string,std::string>> getMonotonicReadAccess
  (std::string uf,MonotonicType type,UniQuantRule* rule,Exp* ex){
   return {};
} 
      

// Function returns write accesses for code generated 
// in a montonic statement.
std::vector<std::pair<std::string,std::string>> 
      getMonotonicWriteAccess(std::string uf,MonotonicType type ,
	 UniQuantRule* rule,Exp* ex){
    return {};
}



std::string getMonotonicStmt(std::string uf,MonotonicType type,
                  Exp* montonicDiffExp){
    return "";
}




Set* CodeSynthesis::GetMonotonicDomain(std::string uf, MonotonicType type,
      Exp* monotonicDiff,Set* ufDomain){
   if(ufDomain->arity() !=  1){
      throw assert_exception("GetMonotonicDomain:"
		      " Domain must have a single arity");
   }
    
    Set* res = new Set(*ufDomain);
    TupleVarTerm t(1,0);
    for(auto it = res->conjunctionBegin();
		    it != res->conjunctionEnd(); it++ ){
    
        std::list<Exp*> upperBounds = (*it)->GetUpperBounds(t);
        for(Exp* upperBound : upperBounds){
            //Add a constraints that limits the upper bound
	    //by the monotonic diff
	    // i < upperBound - monotonicDiff
	    // Equiv: -i + upperBound - monotonicDiff >= 0
	    Exp * e1 = new Exp();
	    Term* tupClone = t.clone();
	    tupClone->setCoefficient(-1);
            e1->addTerm(tupClone);
	    e1->addExp(upperBound);
	    Exp* monClone = monotonicDiff->clone();
	    monClone->multiplyBy(-1);
	    e1->addExp(monClone);
	    e1->setInequality();
	    (*it)->addInequality(e1);
        }
	//Simplify constraint
	(*it)->cleanUp();
    }

    return res;
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

