//
// Created by Tobi Popoola on 10/7/20.
//

#ifndef CODESYNTHESIS_CODESYNTHESIS_H
#define CODESYNTHESIS_CODESYNTHESIS_H
#include <iegenlib.h>
#include <memory>
#include <utility>
#include <vector>
#include <computation/Computation.h>
#include <iegenlib/set_relation/environment.h>
#include <iegenlib/set_relation/UninterpFunc.h>
namespace code_synthesis {
#define PERMUTE_NAME "P"
// Number of times to retry synthesis
// before failing.
#define MAX_TRIES 200
/// Class contains functionality to generate SPF Computation
/// from a Relation and a Set. The relation is a mapping from
/// a previous space to a new space.


/// Synthesizing statements for UFs in expressions results in 4 cases
/// Case 1
/// UF(x) = y
/// If rhs only has one term and the term is an output tuple
/// var term.
/// Case 2
///UF(x) = F(x)
/// RHS of UF expression does not involve an output tuple variable
/// Case 3
/// UF(x) <= F(y)
/// where arity(x) < arity(y) and y is not an output tuple variable.
/// Case 4
/// UF(x) >= F(y)
/// where arity(x) < arity(y) and y is not an output tuple variable.
/// Case 5
/// UF(y) = F(x)
/// arity(y) == 1 , and
/// SELF_REF
/// Special case where the UF
/// References itself
typedef enum {
    CASE1,
    CASE2,
    CASE3,
    CASE4,
    CASE5,
    MERGECASE,
    SELF_REF,
    UNDEFINED
} SynthExpressionCase;


struct UFQuant;

struct SparseFormat {
    std::string mapToDense;
    std::string dataName;
    std::string dataConstraint;
    std::string dataAccess;
    std::vector<UFQuant> ufQuants;
    std::vector<std::string> knowns;
};
typedef std::pair<std::shared_ptr<Set>,std::shared_ptr<Set>>
     PropertyPair;
typedef std::vector<PropertyPair>
     PropertyPairVector;
/*class ReorderTupleInfo holds information on tuple variables
 * present in a map, its range, and its list of ordering constraints
 * this is useful for making decisions to generate permutations
 * or not depending if source and destination formats already have
 * similar ordering.
 * */
struct ReorderTupleInfo{
    int tupleIndex;
    std::unique_ptr<Set> range;
    // Contains all the lhs and rhs properties
    // extracted from uninterprted function 
    // quantifiers. First of the pair has the lhs
    // and second of each pair has the rhs, 
    //
    // Memory is managed using a unique pointer
    PropertyPairVector propertyPairs;
};

//Helper data structure specifying a universal
//quantifier
struct UFQuant {
    std::string domain;
    std::string range;
    std::string name;
    bool isBijective;
    std::string lhsProperty;
    std::string rhsProperty;
    MonotonicType type;
    UFQuant(std::string domain, std::string range, std::string name,
            bool isBijective,MonotonicType type):
        domain(domain), range(range), name(name), isBijective(isBijective),
        type(type) {};

    UFQuant(std::string domain, std::string range, std::string name,
            bool isBijective,MonotonicType type,std::string lhsProperty,
            std::string rhsProperty):
        domain(domain), range(range), name(name), isBijective(isBijective),
        type(type),lhsProperty(lhsProperty),rhsProperty(rhsProperty) {};
};

class CodeSynthesis {
private:
    // Contains permute and other growing functions.
    std::vector<std::string> permutes;

    // Contains unknown functions
    std::vector<std::string> unknowns;
    
    // Permutes exhibiting merge attributes
    // TODO: Mege permutes with CASE1 and CASE2
    // CASE1, p0->insert({t1,t2})
    // CASE2, p0(t1,t2) = t3
    // Merged Permute
    // P0 = Permute(2);
    // p0->insert({t1,t2,t3})
    std::vector<std::pair<std::string,int>> mergedPermutes;

    Relation* sourceMapR;
    Relation* destMapR;

    // Source and destination
    std::string sourceDataName;
    std::string destDataName;
    std::string sourceDataConstraint;
    std::string  destDataConstraint;

    // Contains only ufquantifiers
    // from destination format.
    // These are the quantifiers that
    // has to be satisfied.
    std::vector<UFQuant> ufQuants;

    Relation* sourceDataAccessMap;
    Relation* destDataAccessMap;


    // Compose realtion
    Relation* composeRel;
    // Transitive closure applied on
    // composeRel.
    Relation* transRel;


    // Expanded equalities
    // applied on transRel
    Relation* transRelExpanded;

    // These are ufs or symbolic constants
    // that are already known and need not
    // to be categorized as an unknown for
    // code generation. A good example is
    // NR,NC, NNZ
    std::vector<std::string> knowns;

    /// Generates inspector computation.
    /// \return Computation object, representing the synthesis.
    Computation*  generateInspectorComputation ();

    /// Generates a comparator based on self referential
    /// constraint.
    /// \param e       self refrential expression
    /// \param permute UF been considered
    std::string GenerateSelfRefPermuteConditions(iegenlib::Exp* e,
            std::string& permute);

    // Vector of pairs holds pairs of ufs and self
    // referential constraints. These are interesting constraints
    // for generating monotonicity code later on and getting
    // constraint describing sorting for P. Expressions stored
    // here are expected to be a clone
    std::vector<std::pair<std::string,iegenlib::Exp*>> selfRefs;

    // Create IR component generates an IR specification
    // for an unknown.
    // \param currentUF         Current UF being considered.
    // \param comp              The computationa the IR Component will be added to
    // \param executionSchedule The execution schedule or position of the IR component
    // \param synthCase         Case of the current uf
    // \param exp               expression/constraint involving uf
    // \param unknowns          List of unknowns.
    // \param reorderFunction   Reordering function(permute) if any
    void CreateIRComponent(std::string currentUF,Computation* comp,
                           int executionSchedule, SynthExpressionCase ufCase,
                           iegenlib::Exp* exp,std::vector<std::string>& unknowns,
                           iegenlib::Set* transSet, UFCallTerm* reorderUF=NULL );

    // Create IR component for sort and adds to the computation
    // \param current UF        Current UF being considered.
    // \param comp		The computation the IR component will be added to
    // \param executionSchedule The execution schedule or position of the IR component.
    void CreateSortIRComponent(UFCallTerm* currentUF,Computation* comp,
                               int executionScheduleIndex);

    // Function checks if the domain is bounded by unkowns
    // \param term   term being investigated.
    // \param unknowns List of unknowns
    // \param relation Relation being considered.
    static bool IsDomainBoundedByUnknown(Term* term,
                                         const std::vector<std::string>& unknowns,iegenlib::Relation* rel );
public:
    CodeSynthesis(SparseFormat* source, SparseFormat* dest);

    ~CodeSynthesis();

    std::string generateFullCode( std::vector<int>& fuseStmts,int level);


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

    
    // Function returns the expression the term is equal to.
    // \param  sc sparse constraint. not adopted
    // \param  t  term being serched for. adopted
    // \return expresssion is owned by the caller
    static Exp* getEqualExpr(SparseConstraints* sc, const Term * t);


    static std::string constraintToStatement(Exp* constraint,
            std::string unknownUF, const TupleDecl& tupDecl,
            SynthExpressionCase expCase);

    static Term*
    findCallTerm(Exp* exp, std::string ufName);

    /// Returns supporting macros and
    //  data structures for abstractions
    static std::string getSupportingMacros();

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

    // This function substitutes direct equalities on output tuple
    // if the rhs of the equality is a function of the input tuple.
    // Example:
    //     R {[n] -> [i,k]: row(n) = i ^ P(i,k) = u }
    //     it generates
    //     R {[n] -> [i,k]: row(n) = i ^ P(row(n_
    static iegenlib::Relation* substituteDirectEqualities(iegenlib::Relation* rel);

    /// Function gets expressions involving term.
    static std::vector<iegenlib::Exp*> getExpInvolvingTerm(
        iegenlib::TupleVarTerm*);

    /// Function solves for output tuple in a relation and
    /// gets all constraints involving the output tuple.
    static iegenlib::Relation* solveForOutputTuple(iegenlib::Relation*);


    /// Function returns a list of read access on a UF.
    static std::vector<std::pair<std::string,std::string>> GetReads(
                std::string uf,iegenlib::Exp* constraint,
                SynthExpressionCase expCase,int arity);

    // Function checks if an expression depends on output
    // tuple returns true if it does.
    // \param arity   arity of the constraint
    // \param inArity input arity of the constraint
    // \param e       expression
    // Params are not adopted
    static bool dependsOnOutputTuple(int arity, int inArity,
                                     iegenlib::Exp*e);
    /// Function returns a list of write accesses.
    static std::vector<std::pair<std::string,std::string>> GetWrites(
                std::string uf, iegenlib::Exp* constraint,
                SynthExpressionCase expCase,int arity);

    /// Function returns case of expression as regards a UF.
    /// Case 1
    /// UF(x) = y
    /// If rhs only has one term and the term is an output tuple
    /// var term.
    /// Case 2
    ///UF(x) = F(x)
    /// RHS of UF expression does not involve an output tuple variable
    /// Case 3
    /// UF(x) <= F(y)
    /// where arity(x) < arity(y) and y is not an output tuple variable.
    /// Case 4
    /// UF(x) >= F(y)
    /// where arity(x) < arity(y) and y is not an output tuple variable.
    /// Assumes that the UF exists in the expression.
    //
    //\param constraint constraint under investigation
    //\param unknownUF  associated unknown uf 
    //\param inputArity size of the input arity of mapping
    //\param tupleSize overall tuple size of the mapping
    //\param resolvedOutputTuples output tuple variables that 
    //                            resolve to a function of known ufs
    //                            and tuple variables. 
    static SynthExpressionCase GetUFExpressionSynthCase(Exp* constraint,
            std::string unknownUF, int inputArity, int tupleSize,
	    const std::vector<int> resolvedOutputTuples 
	    = std::vector<int>());

    /// Function removes constraints present in symbNames
    /// from sparse constraints sc
    static void RemoveSymbolicConstraints(const std::vector<std::string>& symbNames,
                                          SparseConstraints* sc);

    /// This returns the domain of a case for synthesis.
    /// Case 1
    /// UF(x) = y
    /// Domain = Const(x) - ( Const(UF(x) = y) U Const(<other unknown ufs>))
    /// where Const means constraints on some vector
    ///
    /// Case 2
    ///UF(x) = F(x)
    /// RHS of UF expression does not involve an output tuple variable
    /// Case 3
    /// UF(x) <= F(y)
    /// where arity(x) < arity(y) and y is not an output tuple variable.
    /// Case 4
    /// UF(x) >= F(y)
    /// where arity(x) < arity(y) and y is not an output tuple variable.
    /// Assumes that the UF exists in the expression.
    static Set* GetCaseDomain(std::string ufName,Set* s,
                              Exp* constraint,SynthExpressionCase expCase);

    // Helper function to add data spaces to a computation
    // object.
    // comp will be modified.
    // baseType the base type of the value of access.
    static void addToDataSpace(Computation& comp,
                               std::vector<std::pair<std::string,std::string>> access,
                               std::string baseType);

    // Function returns read accesses for code generated
    // in a montonic statement.
    static std::vector<std::pair<std::string,std::string>>
            getMonotonicReadAccess(std::string uf,MonotonicType type);


    // Function returns write accesses for code generated
    // in a montonic statement.
    static std::vector<std::pair<std::string,std::string>>
            getMonotonicWriteAccess(std::string uf,MonotonicType type
                                   );

    // This function automatically returns an expression that
    // gives abs(e2 - e1).
    // Example e1 < e2 <=> func(e1) * func(e2), function returns e2 - e1
    // Parameter ex is not adopted.
    // returned expression is owned by caller, caller should deallocate
    static Exp* getMonotonicDiff(std::string uf,Exp* ex);

    // Returns a statement for monotonic type.
    // \param uf current uninterpreted function under consideration
    // \param type monotonic type of uf
    static std::string getMonotonicStmt(std::string uf,iegenlib::MonotonicType type);


    // Function returns the domain of a monotic type
    // \param uf current uf been considered
    // \param type monotonic type
    // \param uQ not adopted
    // Set returned is owned by caller and should be deallocated.
    static Set* GetMonotonicDomain(std::string uf, MonotonicType type,
                                   Set* ufDomain);

    // Function adds permutation constraint to a relation
    // Example
    //     rel = {[i,j] -> [k]}
    //     P(i,j) =  k
    //     rel = {[i,j] -> [k]: P(i,j) = k}
    // For multiple Tuple output variable, permutation becomes
    //     rel = {[i,j] ->[i,k]}
    //     P0(i,j) = i ^ P1(i,j) = k
    // \param rel the relation
    // \param excludeList list of tuple variables
    //                    permutations shouldnt be generated for.
    // Returns a list of permutation names: P0,P1,....,Pn
    //
    static std::vector<std::string> AddPermutationConstraint
    (Relation* rel);


    //Function gets header string for successful compilation
    // of generated code. It includes data structure abstractions
    // for permutation.
    static std::string GetSupportHeader();

    // Function returns statement for data copy
    static std::string GetCopyStmt(std::string sourceDataName, std::string destDataName,
                                   Relation* destMapR, Relation* sourceMap);

    std::vector<std::pair<std::string,std::string>>
            getCopyWriteAccess();


    std::vector<std::pair<std::string,std::string>>
            getCopyReadAccess();
    // Function removes an expression from a constraint
    static void RemoveConstraint(SparseConstraints* sp, Exp* e);

    // Currently unimplemented. This function will determine if
    // a tuple variable is bounded or have constraints
    // involving unknown ufs or symbolic constants. this
    // functionality is usefull in determining what version
    // of permutation we will be synthesizing
    static bool IsTupleBoundedByUnknown(TupleVarTerm& t,
                                        SparseConstraints* sp,
                                        std::vector<std::string>& unknowns);


    // Returns a comparator for a permute. Based on self referential
    // heuristics and user defined functions.
    // \param Permute    permute string been considered.
    // \param composeRel composed relation
    // \param ufQuants   all uf quantifiers
    static std::string GeneratePermuteConditions(std::string& permute,
            iegenlib::Relation* composeRel,
            std::vector<UFQuant>& ufQuants);

    // This function agressively fuses loops with true dependency
    // in order.
    static void ReadReductionFusionOptimization(Computation* comp,
            std::vector<int>& fuseStmts,int level );

    // Remove statements that write to the same point in memory
    static void RedundantStatementElimination(Computation* comp);

    // Simplify constraints and also optimizes out statements
    // involving such constraints.
    static void ConstraintSimplification(Computation* comp);
    
    // Invert Reorder function. Creates an iteration space
    // that iterates throught the inverse of the reorder
    // function.
    static Set* GetInverseIterationSpace(Set* set, UFCallTerm* domUF);
    
    // Function returns a list of output tuple locations that is a 
    // function of input tuple locations and or known UF.
    //
    // Example
    // [i,j]->[k]: k = i
    // Tuple var K is added to the list because
    // k is resolvable by an input tuple variable.
    //
    // This helper function is needed to select candidates for 
    // synthesis that may be overlooked when it involves an output
    // tuple. That way we can be sure that we are making the right 
    // choice. Note that transitive closure circumvents this, as 
    // candidates will still be picked up. However, we will have 
    // a situation where the iteration space is reduced due to the
    // direct assignment which would prevent opportunities for 
    // fusion. Here is an example of such case:
    //
    // TODO: Write an example
    static std::vector<int> GetResolvedOutputTuples(Relation* rel,
		    const std::vector<std::string>& unknownUFs );
    // This function checks to see if this iteration space 
    // is valid for code generation
    // \param set iteration space been considered 
    static bool IsValidIterationSpace(Set* set);
    
    //This function computes information on tuple variables that 
    //has some form of ordering constraint imposed due to universal
    //quantifiers on uninterpreted functions. Function also returns
    //the range of the tuple variable, and the reordering constraint 
    //
    //Memory automatically managed with smart ptrs.
    static std::vector<std::shared_ptr<ReorderTupleInfo>> 
	    ReorderInfoAnalysis(Relation* mapToDense, 
			    std::vector<UFQuant>& ufQuants); 
    
    // This function returns true if the expression e
    // depends on any tuple variable in the range 
    // start and end, start and end inclusive [start,end]
    static bool DependsOnTupleRange(Exp * e, int start , int end);
   
    // This function returns a range that represents the upper and 
    // lower bound of a tuple variable term
    //
    // \param tupTerm (not adopted) 
    // \param conj    conjunction
    //
    // \returns set owned by the caller
    static Set* GetTupleTermRange(TupleVarTerm* tupTerm, Conjunction* conj);
    
    // Function adds new expressions that directly equate tuples
    // in tupleRel together.
    // Throws an exception if the tuple indices in tuple Rels is out
    // of bound.
    static void AddTupleRelationships(SparseConstraints* sc,
		    std::vector<std::pair<int,int>>& tupleRels);

    // Function removes constraints invovling tuple variables
    // and UF from constraint
    static void RemoveConstraintsInvTupleUF(SparseConstraints* sc ,
		    TupleVarTerm& tv, std::string& ufName);
};
}


#endif //CODESYNTHESIS_CODESYNTHESIS_H
