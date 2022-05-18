/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Markus Iser, KIT - Karlsruhe Institute of Technology

Thanks for the original sanitizer code to Marijn Heule

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

#ifndef SRC_TRANSFORM_SANITIZE_H_
#define SRC_TRANSFORM_SANITIZE_H_

#include <stdio.h>
#include <stdlib.h>

#include "src/util/StreamBuffer.h"


void sanitize(const char* filename) {
    StreamBuffer cnf(filename);

    cnf.skipWhitespace();
    while (*cnf == 'c') {
        cnf.skipLine();
    }

    int nv = 0, nc = 0, rv = 0, rc = 0;
    cnf.skipWhitespace();
    if (!cnf.eof() && *cnf == 'p') {
        ++cnf;
        cnf.skipWhitespace();
        cnf.skipString("cnf\0");
        nv = cnf.readInteger();
        nc = cnf.readInteger();
        cnf.skipLine();
        std::cout << "p cnf " << nv << " " << nc << std::endl;
    } else {
        std::cerr << "Error: expected header of the form: 'p cnf [nvars] [nclauses]'" << std::endl;
        return;
    }

    int *mask, stamp = 0;
    mask = (int*) malloc (sizeof (int) * (2*nv + 1));
    mask += nv;

    std::vector<int> clause;
    while (!cnf.eof()) {
        cnf.skipWhitespace();
        if (!cnf.eof() && *cnf != 'c' && *cnf != 'p') {
            bool tautology = false;
            stamp++;

            // read clause without duplicate literals
            for (int plit = cnf.readInteger(); plit != 0; plit = cnf.readInteger()) {
                rv = std::max(abs(plit), rv);
                if (mask[plit] != stamp) {
                    mask[plit] = stamp; 
                    clause.push_back(plit);
                    if (mask[-plit] == stamp) {
                        tautology = true;
                    }
                }
            }

            // print non-tautologic clause
            if (!tautology) {
                rc++;
                for (int lit : clause) {
                    std::cout << lit << " ";
                } 
                std::cout << "0" << std::endl;
            }
            clause.clear();
        }
    }

    if (rc != nc || rv > nv) {
        std::cout << "c header-fix p cnf " << rv << " " << rc << std::endl;
    }

    free(mask);

}

#endif  // SRC_TRANSFORM_SANITIZE_H_