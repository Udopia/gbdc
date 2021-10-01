/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Markus Iser, KIT - Karlsruhe Institute of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#include <string>
#include <vector>

#include "util/CNFFormula.h"

void generate_independent_set_problem(std::string filename) {
    CNFFormula F;
    std::vector<std::vector<unsigned>> occ;
    F.readDimacsFromFile(filename.c_str());
    occ.resize(2 * F.nVars() + 2);
    unsigned nNodes = 0;
    unsigned nEdges = 0;
    unsigned node = 0;
    for (Cl* clause : F) {
        nNodes += clause->size();
        nEdges += (clause->size() * (clause->size() - 1)) / 2;
        for (unsigned i = 0; i < clause->size(); i++) {
            unsigned var1 = node + i + 1;
            occ[(*clause)[i]].push_back(var1);
        }
        node += clause->size();
    }
    for (unsigned i = 1; i <= F.nVars(); i++) {
        nEdges += occ[Lit(Var(i), false)].size() * occ[Lit(Var(i), true)].size();
    }
    std::cout << "c satisfiable iff independent set size is " << F.nClauses() << std::endl;
    std::cout << "p edge " << nNodes << " " << nEdges << std::endl;
    node = 0;
    for (Cl* clause : F) {
        for (unsigned i = 0; i < clause->size(); i++) {
            unsigned var1 = node + i + 1;
            for (unsigned j = i; j < clause->size(); j++) {
                unsigned var2 = node + j + 1;
                std::cout << var1 << " " << var2 << " 0" << std::endl;
            }
        }
        node += clause->size();
    }
    for (unsigned i = 1; i <= F.nVars(); i++) {
        for (unsigned node1 : occ[Lit(Var(i), false)]) {
            for (unsigned node2 : occ[Lit(Var(i), true)]) {
                std::cout << node1 << " " << node2 << " 0" << std::endl;
            }
        }
    }
}