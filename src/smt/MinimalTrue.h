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
            /// Function recursively visits an expression
            /// and flattens Expressions that are not \param unknownTerm
            /// to be placed on the right hand side.
            /// \param currentExpr
            /// \param unknownTerm
            /// \return
            static std::vector<Exp*> getRHS(
                    Exp* currentExpr, Term *unknownTerm);

            /// Function flattens a sparse constraint : set, relation
            /// to individual terms
            /// \param sparseConstraint
            /// \return
            static std::list<Term*>  getTermList(
                    SparseConstraints * sparseConstraint);
        public:
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
        };

        class Stmt {

        public:
            Exp* rhs;
            Exp* lhs;
            std::string toString() const;

        };
    }

}



#endif //CODESYNTHESIS_MINIMALTRUE_H
