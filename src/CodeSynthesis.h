//
// Created by Tobi Popoola on 10/7/20.
//

#ifndef CODESYNTHESIS_CODESYNTHESIS_H
#define CODESYNTHESIS_CODESYNTHESIS_H
#include <iegenlib.h>

#include <utility>
#include <vector>
#include <computation/Computation.h>
namespace code_synthesis {
  /// Class contains functionality to generate SPF Computation
  /// from a Relation and a Set. The relation is a mapping from
  /// a previous space to a new space.

  struct TupleLexOrder {
      bool operator()(const TupleVarTerm* a, const TupleVarTerm* b) 
		const {
          return (*a) < (*b);
      }
  };

  struct Stmt {
      TupleDecl tupleDecl;
      Stmt(TupleDecl tupleDecl):tupleDecl(tupleDecl){}
      Exp* rhs;
      Exp* lhs;
      std::string toString() const;

      std::string toPrettyPrintString() const;
  };

  class CodeSynthesis {
  private:
      Relation* mapToNewSpace;
      Set* originalSpace;


  public: 
      CodeSynthesis(std::string mapToNewSpace, std::string originalSpace):
          mapToNewSpace(new Relation(std::move(mapToNewSpace))),
          originalSpace(new Set(std::move(originalSpace))){};

      ~CodeSynthesis();

      /// Generates inspector computation.
      /// \return Computation object, representing the synthesis.
      Computation*  generateInspectorComputation ();

      /// This gets the list of all expressions in a conjunction.
      /// Each expression in this list is newly allocated and
      /// the caller is responsible for deallocating
      /// \param conjunction
      /// \return
      static std::list<Exp*> getExprs(Conjunction* conjunction);



      /// Intersects two lists
      /// \tparam T type of the list
      /// \param a first list
      /// \param b second list
      /// \return returns a list intersection
      template <typename  T>
      static std::list<T*> intersectLists(const std::list<T*>& a,
                                          const std::list<T*>&b);

      /// Compares terms absolutely without regards for
      //  their coefficients.
      /// \param a term a
      /// \param b term b
      /// \return true if abs(a) is equal to abs(b)
      static bool compareAbsTerms(const Term * a, const Term* b);


      /// Function synthesizes list of statements for an unknown
      /// term.
      /// \param unknownTerm term to be synthesized
      /// \return
      std::list<Stmt *> synthesizeStatements(Term *unknownTerm);


      static std::string constraintToStatement(Exp* constraint,
		      std::string unknownUF, int inputArity, int tupleSize);

      static UFCallTerm*
	      findCallTerm(Exp* exp, std::string ufName);



      /// This function uses a transformation relation and a set to
      /// identify unknowns in a transformation. This is an important
      /// step in code synthesis.
      /// \return list of unknown terms in transform relation
      std::list<Term *> evaluateUnknowns();

      /// Function flattens a sparse constraint : set, relation
      /// to individual terms
      /// \param sparseConstraint
      /// \return
      static std::list<Term*>  getTermList(
              SparseConstraints * sparseConstraint);

      /// Function flattens Parameters in a UF
      /// to individual terms
      /// \param ufCallTerm
      /// \return List of terms in ufcall Parameters
      static std::list <Term*>getParamTermList(
			const UFCallTerm *ufCallTerm);


      /// Returns a minimally true expression from
      /// an inequality expression.
      /// \param expr
      /// \return
      static Exp* getMinTrueExpr(Exp* expr);


      /// Check if a term is contained in a list of terms. This
      /// goes recursively deeper into UFCall terms.
      /// \param terms list of terms
      /// \param term  been searched for
      /// \return true if term is in there and returns false otherwise.
      static bool containsTerm (const std::list<Term*>&terms,
                                const Term * term);


      /// Counts number of tuple variables in a term list.
      /// \param terms
      /// \return
      static int getTupleVarCount(std::list<Term*>& terms);


      /// This function gets the domain of an unknown
      /// Term in a relation
      /// \param unknownTerm unkown term currently being investigated
      /// \param unkownTerms unknown terms in the relation
      /// \throw Exception if relation has no constraint
      /// \return
      Set *getDomain(Term *unknownTerm, std::list<Term *> &unkownTerms);

      /// This function extracts dependents to term from a
      /// conjunction, a term x is a dependent to term y
      /// if term x exists in some constraint involving term y
      /// \param conjunction
      /// \param term
      /// \return
      static std::list<Term*> getDependents(Conjunction * conjunction,
			Term* term);

      /// This returns a string to allocate memory for an unknown
      /// term.
      /// \param unknownTerm
      /// \return
      static std::string getAllocationStmt (Term* unknownTerm);

      /// Gets formatted tuple string from a list
      /// \param list
      /// \return a string of formatted tuple strings
      static std::string getFormattedTupleString
                          (const std::list<std::string>& list);
      
       
      /// Function creates an execution schedule from a set
      /// and sets the outermost loop to position pos
      /// Example:
      ///   {[n,k]: C1 ^ C2}, pos=1
      ///   {[n,k] -> [1,n,0,k,0]}
      /// This functionality helps with code synthesis.
      static iegenlib::Relation* getExecutionSchedule(iegenlib::Set*, int pos); 


      /// Function gets expressions involving term.
      static std::vector<iegenlib::Exp*> getExpInvolvingTerm(
          iegenlib::TupleVarTerm*);

      /// This function solves for output tuple in a relation and
      /// gets all constraints involving the output tuple.
      static iegenlib::Relation* solveForOutputTuple(iegenlib::Relation*);



  };
}


#endif //CODESYNTHESIS_CODESYNTHESIS_H
