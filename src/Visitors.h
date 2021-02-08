#ifndef VISITORS_H
#define VISITORS_H

#include  <iegenlib/set_relation/Visitor.h>
#include  <list>
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

#endif
