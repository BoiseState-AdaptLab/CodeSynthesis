#include "Visitors.h"

void TermVisitor::preVisitTerm(iegenlib::Term * t){
  if (!t->isConst())
    terms.push_back(t);  
}

void TermVisitor::preVisitUFCallTerm(iegenlib::UFCallTerm * t){ 
  terms.push_back(t);  
}

void TermVisitor::preVisitTupleVarTerm(iegenlib::TupleVarTerm * t){ 
  terms.push_back(t);  
}

void TermVisitor::preVisitVarTerm(iegenlib::VarTerm * t){ 
  terms.push_back(t);  
}
void TermVisitor::preVisitTupleExpTerm(iegenlib::TupleExpTerm * t){
  terms.push_back(t);
}	

void ExpTermVisitor::preVisitTupleVarTerm(iegenlib::TupleVarTerm *t){
    if(t->tvloc() == term->tvloc()){
       auto expression= expStack.top();
       auto it = std::find(exps.begin(),exps.end(), expression);
       if(it == exps.end()){
           exps.push_back(expression);
       }
    }
}

std::vector<iegenlib::Exp*> ExpTermVisitor::getExpressions(){
    return exps;
}	

void ExpTermVisitor::preVisitExp(iegenlib::Exp* e){
    expStack.push(e);
}

void ExpTermVisitor::postVisitExp(iegenlib::Exp* e){
    expStack.pop();
}
