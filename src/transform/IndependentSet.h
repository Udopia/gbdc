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
#ifndef SRC_TRANSFORM_INDEPENDENTSET_H_
#define SRC_TRANSFORM_INDEPENDENTSET_H_

#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include "src/util/CNFFormula.h"

void generate_independent_set_problem(std::string filename, std::string output) {
    CNFFormula F;
    std::vector<std::vector<unsigned>> literal2nodes;
    F.readDimacsFromFile(filename.c_str());
    std::shared_ptr<std::ostream> of;
    if (output != "-") {
        of.reset(new std::ofstream(output.c_str(), std::ofstream::out));
    } else {
        of.reset(&std::cout, [](...){});
    }

    literal2nodes.resize(2 * F.nVars() + 2);
    unsigned nNodes = 0;
    unsigned nEdges = 0;
    unsigned nodeId = 1;
    for (Cl* clause : F) {
        nNodes += clause->size();  // one node per literal occurence
        nEdges += (clause->size() * (clause->size() - 1)) / 2;  // number of edges in clique
        for (unsigned i = 0; i < clause->size(); i++) {
            literal2nodes[(*clause)[i]].push_back(nodeId + i);  // remember nodeids of literals
        }
        nodeId += clause->size();
    }
    for (unsigned i = 1; i <= F.nVars(); i++) {
        // count edges between nodes for opposite literals
        nEdges += literal2nodes[Lit(Var(i), false)].size() * literal2nodes[Lit(Var(i), true)].size();
    }
    *of << "c satisfiable iff maximum independent set size is " << F.nClauses() << std::endl;
    *of << "c kis nNodes nEdges k" << std::endl;
    *of << "p kis " << nNodes << " " << nEdges << " " << F.nClauses() << std::endl;
    nodeId = 0;
    // generate cliques
    for (Cl* clause : F) {
        for (unsigned i = 0; i < clause->size(); i++) {
            unsigned var1 = nodeId + i;
            for (unsigned j = i; j < clause->size(); j++) {
                unsigned var2 = nodeId + j;
                *of << var1 << " " << var2 << " 0" << std::endl;
                *of << var2 << " " << var1 << " 0" << std::endl;
            }
        }
        nodeId += clause->size();
    }
    // generate edges between nodes for opposite literals
    for (unsigned i = 1; i <= F.nVars(); i++) {
        for (unsigned node1 : literal2nodes[Lit(Var(i), false)]) {
            for (unsigned node2 : literal2nodes[Lit(Var(i), true)]) {
                *of << node1 << " " << node2 << " 0" << std::endl;
                *of << node2 << " " << node1 << " 0" << std::endl;
            }
        }
    }
}

#endif  // SRC_TRANSFORM_INDEPENDENTSET_H_
