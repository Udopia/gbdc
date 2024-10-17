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

#ifndef SRC_TRANSFORM_BIP_H_
#define SRC_TRANSFORM_BIP_H_

#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include "src/util/CNFFormula.h"

class BipartiteGraphFromCNF {
 private:
    CNFFormula F;

 public:
    explicit BipartiteGraphFromCNF(const char* filename) : F() {
        F.readDimacsFromFile(filename);
    }

    void generate_bipartite_graph(const char* output = nullptr) {
        std::shared_ptr<std::ostream> of;
        if (output != nullptr) {
            of.reset(new std::ofstream(output, std::ofstream::out));
        } else {
            of.reset(&std::cout, [](...){});
        }

        *of << "c directed bipartite graph representation from cnf" << std::endl;
        *of << "p edge " << F.nVars() + F.nClauses() << std::endl;

        unsigned clause_id = F.nVars() + 1;
        for (Cl* clause : F) {
            for (unsigned i = 0; i < clause->size(); i++) {
                if ((*clause)[i].sign()) {
                    *of << "e " << (*clause)[i].var() << " " << clause_id << std::endl;
                } else {
                    *of << "e " << clause_id << " " << (*clause)[i].var() << std::endl;
                }
            }
            clause_id++;
        }
    }
};

#endif  // SRC_TRANSFORM_BIP_H_
