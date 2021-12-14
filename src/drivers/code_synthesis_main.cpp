#include <iostream>
#include <string>
#include <iegenlib/computation/Computation.h>
#include <CodeSynthesis.h>
#include <iegenlib/set_relation/set_relation.h>
int main() {
    Computation inspector;
    // Add universal constraints
    
    iegenlib::Set * colDomain = new iegenlib::Set("{[x]:0 <= x < NNZ}");
    iegenlib::Set * colRange = new iegenlib::Set("{[j]:0 <= j < NC}");
    
    
    iegenlib::Set * rowptrDomain = new iegenlib::Set("{[i]:0 <= i <= NR}");
    iegenlib::Set * rowptrRange = new iegenlib::Set("{[x]:0 <= x < NNZ}");
    
    iegenlib::setCurrEnv();
    iegenlib::appendCurrEnv("col1",colDomain,colRange,
		    false,iegenlib::Monotonic_NONE);
    iegenlib::appendCurrEnv("rowptr",rowptrDomain,rowptrRange,
		    false,iegenlib::Monotonic_Increasing);
    // End Universal constraints
    
    
    int executionScheduleIndex  = 0;

    iegenlib::Relation* map1 = 
	    new iegenlib::Relation("{[i,j]->[k]: i >= 0 and i < NR and"
                              " j >= 0 and j < NC and rowptr(i) <= k < rowptr(i+1)"
			      " and j = col1(k) and P(i,j) = k }");
    iegenlib::Relation* map2 =
	   new iegenlib::Relation("{[n] -> [i,j]:"
                          " row(n) = i and 0 <= n and n < NNZ "
			  "and col2(n) = j and  i >= 0 and "
                          " i < NR and j >= 0 and j < NC}");
    iegenlib::Relation* composeRel = map1->Compose(map2);
    iegenlib::Relation* transRel = composeRel->TransitiveClosure(); 
    std::list<iegenlib::Exp*> expList =
    code_synthesis::CodeSynthesis::getExprs(*transRel->conjunctionBegin());    
    
    std::string pStmt;
    iegenlib::Exp* pExp = NULL ;
    // Solve for P
    for(auto e : expList){
        if (e->isEquality() && code_synthesis::
			CodeSynthesis::findCallTerm(e,"P")!=NULL){
	    pStmt = code_synthesis::CodeSynthesis::
	         constraintToStatement(e,
		   "P",composeRel->inArity(),composeRel->getTupleDecl());
	    pExp = e;
	    break;
	}

	// TODO: Look for constraints that describes sorting.
	// for help in initializing P DS.
        	
    }
    
    assert(pExp && "Synth Failure: No constraints involving P");
     
    std::vector<std::string> unknowns {"rowptr","col1"};
    
    auto expCase = code_synthesis::CodeSynthesis::GetUFExpressionSynthCase(pExp,
		   "P",composeRel->inArity(),composeRel->arity());
    
    // Convert compose relation to set.
    auto composeSet = composeRel->ToSet();
 
    
    // Get Domain for P
    iegenlib::Set* pDomain = code_synthesis::CodeSynthesis::
	   GetCaseDomain("P",composeSet,pExp,expCase);
    
    // remove constraints involving unknown UFs
    code_synthesis::CodeSynthesis::RemoveSymbolicConstraints(unknowns,pDomain);
    
    // Get execution schedule
    iegenlib::Relation* pExecutionSchedule = code_synthesis::
	    CodeSynthesis::getExecutionSchedule(
			    pDomain,executionScheduleIndex++);

    
    // Get reads and writes.
    auto writes = code_synthesis::CodeSynthesis::
	    GetWrites("P",pExp,expCase,pDomain->arity()); 

    auto reads = code_synthesis::CodeSynthesis::
	    GetReads("P",pExp,expCase,pDomain->arity()); 

    code_synthesis::CodeSynthesis::addToDataSpace(inspector,
			reads, "double");
    
    code_synthesis::CodeSynthesis::addToDataSpace(inspector,
				writes, "double");

    inspector.addStmt(new Stmt(pStmt,pDomain->prettyPrintString(),
			    pExecutionSchedule->prettyPrintString(),
			    reads,writes));
     
    std::vector<iegenlib::Set*> unknownDomain;
    
    for(std::string unknown: unknowns){
        iegenlib::Set* domain = transRel->GetDomain(unknown);
        unknownDomain.push_back(domain);	
    } 
    // Convert trans relation to a set
    Set* transSet = transRel->ToSet();
     
    for (auto currentUF : unknowns){
        // skip UF for P since it has already been 
	// synthesized at this point.
	if(currentUF == "P") continue;

        std::list<iegenlib::Exp*> expUfs;
	std::list<std::string> expStmts;
	for(auto e : expList){
	   if(code_synthesis::CodeSynthesis::findCallTerm(e,currentUF)!=NULL){
	      
             auto ufCase = code_synthesis::CodeSynthesis::GetUFExpressionSynthCase(e,
		   currentUF,transRel->inArity(),transRel->arity());
             if(ufCase !=code_synthesis::UNDEFINED){
	         std::string expStmt = code_synthesis::CodeSynthesis::
	             constraintToStatement(e,
		       currentUF,transRel->inArity(),transRel->getTupleDecl());
		 Set* ufDomain = code_synthesis::CodeSynthesis::GetCaseDomain(
				 currentUF,transSet,e,ufCase);
                 
                 // remove constraints involving unknown UFs
                 code_synthesis::CodeSynthesis::
			 RemoveSymbolicConstraints(unknowns,ufDomain);
                 
		 // Get reads and writes.
                 auto ufWrites = code_synthesis::CodeSynthesis::
	                 GetWrites(currentUF,e,ufCase,ufDomain->arity()); 

                 
		 auto ufReads = code_synthesis::CodeSynthesis::
	                 GetReads(currentUF,e,ufCase,ufDomain->arity()); 
                 
		 // add data spaces for reads and writes to 
		 // IR
	         code_synthesis::CodeSynthesis::addToDataSpace(inspector,
				ufReads, "double");
                  
	         code_synthesis::CodeSynthesis::addToDataSpace(inspector,
				ufWrites, "double");
                 
		 // Get execution schedule
                 iegenlib::Relation* ufExecSched = code_synthesis::
	                 CodeSynthesis::getExecutionSchedule(
			    ufDomain,executionScheduleIndex++);

		 inspector.addStmt(new Stmt(expStmt,ufDomain->prettyPrintString(),
					 ufExecSched->prettyPrintString(),
					 reads,writes));
                 delete ufDomain;
		 delete ufExecSched;
	     }
	   }
	}
	// Synthesize statements 
    }
    // CodeGen (RS2->S1(I)) - Data copy Code
    iegenlib::Set* copyDomain = composeRel->ToSet();
    std::cerr << "copy space: "<< copyDomain->prettyPrintString();
    std::string copyStmt = "ACSR(n,k) = ACOO(n,k)";
    iegenlib::Relation* execSchedule = code_synthesis::
	    CodeSynthesis::getExecutionSchedule(
			    copyDomain,executionScheduleIndex);
    inspector.addStmt(new Stmt(copyStmt,copyDomain->prettyPrintString(),
			    execSchedule->prettyPrintString(),
			    {{"ACOO" , "{[n,k] -> [n]}"}}, 
			    {{"ACSR" , "{[n,k] -> [k]}"}}));
    inspector.padExecutionSchedules();
    inspector.printInfo();
    std::cerr << "ToDot \n" << inspector.codeGen();
    return 0;
}
