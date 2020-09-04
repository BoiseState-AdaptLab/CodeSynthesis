//
// Created by Tobi Popoola on 8/30/20.
//

#ifndef CODESYNTHESIS_MINIMALSATISFIABLITY_H
#define CODESYNTHESIS_MINIMALSATISFIABLITY_H

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
            static std::vector<Exp*> getRHS(Exp* currentExpr, Term *unknownTerm);

        public:
            /// Function synthesizes statement from a constraint.
            /// This statement reorders an unknown expression to
            /// be on the lhs of a constraint equality;
            /// \param constraint.
            /// \param known
            /// \param unknownTerm
            /// \return
            static Stmt * synthStmt(Exp *constraint, Term *unknownTerm);
        };

        class Stmt {

        public:
            std::unique_ptr<Exp> rhs;
            std::unique_ptr<Exp> lhs;
            std::string toString() const;

        };
    }

}



#endif //CODESYNTHESIS_MINIMALSATISFIABLITY_H
