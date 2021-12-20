#ifndef VISITORS_H
#define VISITORS_H

#include  <iegenlib/set_relation/Visitor.h>
#include  <list>
#include <vector>
#include <algorithm>
#include <stack>
/// TermVisitor class extracts all 
/// the terms in a Set/Relation
///
class TermVisitor : public Visitor {
  private: 
    std::list<iegenlib::Term*> terms;
  public:
    void preVisitTerm(iegenlib::Term * t) override;
    void preVisitUFCallTerm(iegenlib::UFCallTerm * t) override;
    void preVisitTupleVarTerm(iegenlib::TupleVarTerm * t) override;
    void preVisitVarTerm(iegenlib::VarTerm * t) override;
    void preVisitTupleExpTerm(iegenlib::TupleExpTerm * t) override;
    std::list<iegenlib::Term*> getTerms(){return terms;}

};

/// This visitor extracts expressions 
/// involving term
class ExpTermVisitor: public Visitor {
  private:
    std::vector<iegenlib::Exp*> exps;
    iegenlib::TupleVarTerm* term;
    std::stack<iegenlib::Exp*> expStack;
  public:
    explicit ExpTermVisitor(iegenlib::TupleVarTerm* term): term(term){}
    void preVisitTupleVarTerm(iegenlib::TupleVarTerm *t) override;
    void preVisitExp(iegenlib::Exp* e) override;
    void postVisitExp(iegenlib::Exp* e) override;
    std::vector<iegenlib::Exp*> getExpressions(); 

};

#endif
