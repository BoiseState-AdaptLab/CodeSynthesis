//
// Created by Tobi Popoola on 8/30/20.
//

#ifndef CODESYNTHESIS_MINIMALSATISFIABLITY_H
#define CODESYNTHESIS_MINIMALSATISFIABLITY_H

/**
 * This class provides routines for resolving
 * minimal satisfiablity of a constraint.
 */
#include <iegenlib/iegenlib.h>
#include <string>
namespace code_synthesis {
    namespace smt {
        class MinimalSatisfiablity {
        public:
            static bool asserMinimal(Exp* constraint, Exp* minimalConstraint );
        };
    }

}



#endif //CODESYNTHESIS_MINIMALSATISFIABLITY_H
