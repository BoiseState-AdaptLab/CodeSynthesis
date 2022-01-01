#include <iostream>
#include <string>
#include <iegenlib/computation/Computation.h>
#include <CodeSynthesis.h>
#include <iegenlib/set_relation/environment.h>
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
		    false,iegenlib::Monotonic_Nondecreasing);
    // End Universal constraints
    
    
    int executionScheduleIndex  = 0;

    iegenlib::Relation* map1 = 
	    new iegenlib::Relation("{[i,j]->[k]: i >= 0 and i < NR and"
                              " j >= 0 and j < NC and rowptr(i) <= k < rowptr(i+1)"
			      " and j = col2(k) and P(i,j) = k }");
    iegenlib::Relation* map2 =
	   new iegenlib::Relation("{[n] -> [i,j]:"
                          " row1(n) = i and 0 <= n and n < NNZ "
			  "and col1(n) = j and  i >= 0 and "
                          " i < NR and j >= 0 and j < NC}");
    iegenlib::Relation* composeRel = map1->Compose(map2);
    iegenlib::Relation* transRel = composeRel->TransitiveClosure(); 
    std::list<iegenlib::Exp*> expList =
    code_synthesis::CodeSynthesis::getExprs(*transRel->conjunctionBegin());    
    
    // Vector of pairs holds pairs of ufs and self
    // referential constraints. These are interesting constraints
    // for generating monotonicity code later on and getting 
    // constraint describing sorting for P. Expressions stored
    // here are expected to be a clone
    std::vector<std::pair<std::string,iegenlib::Exp*>> selfRefs; 
    
    iegenlib::Exp* pExp = NULL ;
    // Solve for P
    for(auto e : expList){
        if (code_synthesis::
			CodeSynthesis::findCallTerm(e,"P")!=NULL){
            auto caseP = code_synthesis::CodeSynthesis::
		    GetUFExpressionSynthCase(e,"P",
				    transRel->inArity(),transRel->arity());
	    if (caseP == code_synthesis::CASE1) { 
	        pExp = e;
	    }
	    else if (caseP == code_synthesis::SELF_REF) {
	        selfRefs.push_back({"P",e});
	    }
	} 

         	
    }
    
    assert(pExp && "Synth Failure: No constraints involving P");
     
    std::vector<std::string> unknowns {"rowptr","col2"};
    
    std::string pStmt = code_synthesis::CodeSynthesis::
	         constraintToStatement(pExp,
		   "P",composeRel->getTupleDecl(),code_synthesis::CASE1);
    
    
    // Convert compose relation to set.
    auto composeSet = composeRel->ToSet();
 
    
    // Get Domain for P
    iegenlib::Set* pDomain = code_synthesis::CodeSynthesis::
	   GetCaseDomain("P",composeSet,pExp,code_synthesis::CASE1);
    
    // remove constraints involving unknown UFs
    code_synthesis::CodeSynthesis::RemoveSymbolicConstraints(unknowns,pDomain);
    
    // Get execution schedule
    iegenlib::Relation* pExecutionSchedule = code_synthesis::
	    CodeSynthesis::getExecutionSchedule(
			    pDomain,executionScheduleIndex++);

    
    // Get reads and writes.
    auto writes = code_synthesis::CodeSynthesis::
	    GetWrites("P",pExp,code_synthesis::CASE1,pDomain->arity()); 

    auto reads = code_synthesis::CodeSynthesis::
	    GetReads("P",pExp,code_synthesis::CASE1,pDomain->arity()); 
    
    
    code_synthesis::CodeSynthesis::addToDataSpace(inspector,
			reads, "double");
    // Writes to P is considered a single data space
    // but read from P is a 2d data space 
    //code_synthesis::CodeSynthesis::addToDataSpace(inspector,
    //				writes, "double");

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
        std::list<iegenlib::Exp*> expUfs;
	std::list<std::string> expStmts;
	for(auto e : expList){
	   if(code_synthesis::CodeSynthesis::findCallTerm(e,currentUF)!=NULL){


	     // Get case a uf in a constraint falls into.	   
	     auto ufCase = code_synthesis::CodeSynthesis::
		     GetUFExpressionSynthCase(e,
		   currentUF,transRel->inArity(),transRel->arity());
	     if (ufCase ==code_synthesis::SELF_REF){ 
	         // Store self referential constraints. Should help
	         // with generating code for montonicity later on.             
	        selfRefs.push_back({currentUF,e->clone()});
	        // skip this loop iteration
		continue;
	     }
	     // IF UF satisifies synthesis case
	     if(ufCase !=code_synthesis::UNDEFINED){
	         std::string expStmt = code_synthesis::CodeSynthesis::
	             constraintToStatement(e,
		       currentUF,transRel->getTupleDecl(),ufCase);
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
					 ufReads,ufWrites));
                 delete ufDomain;
		 delete ufExecSched;
	     }
	   }
	}
	// Synthesize statements 
    }
    // Generate code to ensure universal constraint
    for(auto uf : unknowns){
	iegenlib::MonotonicType type = iegenlib::queryMonoTypeEnv(uf);
	
	if(type==Monotonic_NONE)  continue;
	iegenlib::Set* domain = iegenlib::queryDomainCurrEnv(uf);
        
	
	UniQuantRule* uQ =  iegenlib::
	    getUQRForFuncDomainRange(uf); 

	// Domain sorrunding statement for montonic
	// synthesis
	iegenlib::Set* stmtDomain = code_synthesis::CodeSynthesis::
	    GetMonotonicDomain(uf,type,domain);


	// Reads and writes data accesses have to calculated 
	// specially but for now we use the default getwrites
	// and getreads

        auto ufWrites = code_synthesis::CodeSynthesis::
	             getMonotonicWriteAccess(uf,type); 

                 
        auto ufReads = code_synthesis::CodeSynthesis::
		   getMonotonicReadAccess(uf,type);

         	
        std::string monStmt = code_synthesis::CodeSynthesis::
	           getMonotonicStmt(uf,type);	
        // Get execution schedule
        iegenlib::Relation* execSched = code_synthesis::
	             CodeSynthesis::getExecutionSchedule(
	   stmtDomain,executionScheduleIndex++);
	    
	inspector.addStmt(new Stmt(monStmt,stmtDomain->prettyPrintString()
				    ,execSched->prettyPrintString()
				    ,ufReads,ufWrites));
        delete stmtDomain;
    	delete execSched;
    }
    
    // CodeGen (RS2->S1(I)) - Data copy Code
    iegenlib::Set* copyDomain = composeRel->ToSet();
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
    std::cout << code_synthesis::CodeSynthesis::getSupportingMacros(); 
    std::cout << inspector.codeGen();
    return 0;
}
