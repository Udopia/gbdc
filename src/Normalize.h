/*************************************************************************************************
GBDHash -- Copyright (c) 2020, Markus Iser, KIT - Karlsruhe Institute of Technology

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

#ifndef Normalize_h
#define Normalize_h

#include "src/StreamBuffer.h"
#include <vector>
#include <algorithm>

void normalize(const char* filename) {
    StreamBuffer in(filename);
    std::vector<int> clause;
    while (!in.eof()) {
        in.skipWhitespace();
        if (in.eof()) {
            break;
        }
        else if (*in == 'c') {
            in.skipLine();
        }
        else if (*in == 'p') {
            ++in;
            in.skipWhitespace();
            in.skipString("cnf\0");
            int nv = in.readInteger();
            int nc = in.readInteger();
            in.skipLine();
            std::cout << "p cnf " << nv << " " << nc << std::endl;
        }
        else {
            clause.clear();
            for (int plit = in.readInteger(); plit != 0; plit = in.readInteger()) {
                clause.push_back(plit);
            }
            std::sort(clause.begin(), clause.end(), [](int l1, int l2) { return abs(l1) < abs(l2); });
            clause.erase(std::unique(clause.begin(), clause.end()), clause.end()); // remove redundant literals
            if (std::adjacent_find(clause.begin(), clause.end(), [](int l1, int l2) { return l1 + l2 == 0; }) != clause.end()) continue; // skip tautology
            for (int plit : clause) {
                std::cout << plit << " ";
            }
            std::cout << "0" << std::endl;
        }
    }
}

#endif