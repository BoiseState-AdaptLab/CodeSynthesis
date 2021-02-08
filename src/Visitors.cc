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

