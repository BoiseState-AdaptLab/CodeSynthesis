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
	    new iegenlib::Relation("{[i,j]->[i,k,j]: i >= 0 and i < NR and"
                              " j >= 0 and j < NC and rowptr(i) <= k < rowptr(i+1)"
			      " and j = col1(k) and P(i,j) = k }");
    iegenlib::Relation* map2 =
	   new iegenlib::Relation("{[n,i,j] -> [i,j]:"
                          " row(n) = i and 0 <= n and n < NNZ "
			  "and col2(n) = j and  i >= 0 and "
                          " i < NR and j >= 0 and j < NC}");
    iegenlib::Relation* composeRel = map1->Compose(map2);
    std::cerr << "composed relation: "<< composeRel->prettyPrintString();
     iegenlib::Set* copyDomain = composeRel->ToSet();
    std::cerr << "copy space: "<< copyDomain->prettyPrintString();
    std::string copyStmt = "ACSR(n,k) = ACOO(n,k)";
    iegenlib::Relation* execSchedule = code_synthesis::
	    CodeSynthesis::getExecutionSchedule(
			    copyDomain,executionScheduleIndex);
    inspector.addStmt(new Stmt(copyStmt,copyDomain->toString(),execSchedule->prettyPrintString(),
			    {{"ACOO" , "{[n,k] -> [n]}"}}, 
			    {{"ACSR" , "{[n,k] -> [k]}"}}));

   
    

    iegenlib::Relation* transRel = composeRel->TransitiveClosure(); 
    std::cout <<"Transitive Closure: "<< transRel->prettyPrintString() << "\n";
    // Assuming there is a symbol Iterator 
    std::vector<std::string> unknowns {"rowptr","col1"};

    std::vector<iegenlib::Set*> unknownDomain;
    
    for(std::string unknown: unknowns){
        iegenlib::Set* domain = transRel->GetDomain(unknown);
        unknownDomain.push_back(domain);	
    }
   

    std::list<iegenlib::Exp*> expList =
    code_synthesis::CodeSynthesis::getExprs(*transRel->conjunctionBegin());    
    
    while (unknowns.size() >= 0){
	std::string currentUF = unknowns.front();
        std::list<iegenlib::Exp*> expUfs;
	std::list<std::string> expStmts;
	for(auto e : expList){
	   if(code_synthesis::CodeSynthesis::findCallTerm(e,currentUF)!=NULL){
	      std::string expStmt = code_synthesis::CodeSynthesis::
	         constraintToStatement(e,
		   currentUF,composeRel->inArity(),composeRel->arity());
	      if (expStmt.size() != 0){
	         expUfs.push_back(e);
		 expStmts.push_back(expStmt);
	      }
	   }
	}
	// Synthesize statements 
    }
    // CodeGen (RS2->S1(I)) - Data copy Code
        return 0;
}
