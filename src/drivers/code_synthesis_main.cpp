#include <iostream>
#include <string>
#include <iegenlib/computation/Computation.h>
#include <CodeSynthesis.h>

int main() {
    std::string denseSpace  = "{[i,j]: i >= 0 and i < NR and"
                              " j >= 0 and j < NC and Ad(i,j) > 0}";
    std::string mapFromDenseToCoo = "{[i,j] -> [n]:"
                          " row(n) = i and col(n) = j and  i >= 0 and "
                          " i < NR and j >= 0 and j < NC}";

    code_synthesis::CodeSynthesis* synth = new code_synthesis::
	    CodeSynthesis(mapFromDenseToCoo, denseSpace);
    iegenlib::Computation *comp = synth->generateInspectorComputation();

    comp->printInfo();
    return 0;
}
