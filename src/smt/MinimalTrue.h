//
// Created by Tobi Popoola on 8/30/20.
//

#ifndef CODESYNTHESIS_MINIMALTRUE_H
#define CODESYNTHESIS_MINIMALTRUE_H

#include <iegenlib/iegenlib.h>
#include <string>
#include <vector>
namespace code_synthesis {
    namespace smt {
        class Stmt;

        /// This class provides routines for resolving
        /// minimal satisfiablity of a constraint.
        class MinimalSatisfiablity {
        private:
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



        public:


            struct TermCompare {
                bool operator()(const Term*& a, const Term*& b) const {
                    return (*a)==(*b);
                }
            };

            struct TupleLexOrder {
                bool operator()(const TupleVarTerm* a, const TupleVarTerm* b) const {
                    return (*a) < (*b);
                }
            };

            /// This function creates a mapping from a sourceMap
            /// to a destination map. The output arity of sourceMap
            /// and output arity of the destinationMap must be equal.
            /// \param sourceMaptoDense
            /// \param destinationMaptToDense
            /// \return
            static Relation * mapToNewSpace(Relation* sourceMap,
                                            Relation* destinationMap);

            /// Function synthesizes statement from a constraint.
            /// This statement reorders an unknown expression to
            /// be on the lhs of a constraint equality;
            /// \param constraint.
            /// \param known
            /// \param unknownTerm
            /// \return
            static Stmt * synthStmt(Exp *constraint, Term *unknownTerm);

            /// This function uses a transformation relation and a set to
            /// identify unknowns in a transformation. This is an important
            /// step in code synthesis.
            /// \param transformRelation
            /// \param set
            /// \return list of unknown terms in transform relation
            static std::list<Term*> evaluateUnknowns(Relation * transformRelation,Set* set);

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
            static std::list <Term*>
                    getParamTermList(const UFCallTerm *ufCallTerm);


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
            /// \param relation containing domain information
            /// \param unknownTerm unkown term currently being investigated
            /// \param unkownTerms unknown terms in the relation
            /// \throw Exception if relation has no constraint
            /// \return
            static Set* getDomain (Relation* relation, Term *unknownTerm,
                                   std::list<Term*>& unkownTerms);

            /// This function extracts dependents to term from a
            /// conjunction, a term x is a dependent to term y
            /// if term x exists in some constraint involving term y
            /// \param conjunction
            /// \param term
            /// \return
            static std::list<Term*> getDependents(Conjunction * conjunction,Term* term);



        };
        /// Stmt class holds a synthesized statement.
        /// This synthesized statement holds information
        /// about the execution schedule and domain.
        class Stmt {

        public:
            Set * domain;
            Exp* rhs;
            Exp* lhs;
            std::string toString() const;

        };
    }

}



#endif //CODESYNTHESIS_MINIMALTRUE_H
