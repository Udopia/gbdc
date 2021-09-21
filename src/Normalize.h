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

#ifndef SRC_NORMALIZE_H_
#define SRC_NORMALIZE_H_

#include <vector>
#include <algorithm>

#include "src/StreamBuffer.h"

void normalize(const char* filename) {
    StreamBuffer in(filename);
    int nv = 0, nc = 0, rv = 0, rc = 0;
    while (!in.eof()) {
        in.skipWhitespace();
        if (in.eof()) {
            break;
        } else if (*in == 'c') {
            if (nv > 0 && nc > 0) {
                in.skipLine();  // skip comments after header
            } else {
                while (!in.eof() && (!isspace(*in) || isblank(*in))) {
                    std::cout << *in;
                    ++in;
                }
                in.skipWhitespace();
                std::cout << std::endl;
            }
        } else if (*in == 'p') {
            ++in;
            in.skipWhitespace();
            in.skipString("cnf\0");
            nv = in.readInteger();
            nc = in.readInteger();
            in.skipLine();
            std::cout << "p cnf " << nv << " " << nc << std::endl;
        } else {
            for (int plit = in.readInteger(); plit != 0; plit = in.readInteger()) {
                std::cout << plit << " ";
                rv = std::max(abs(plit), rv);
            }
            std::cout << "0" << std::endl;
            rc++;
        }
    }
    if (rc != nc || rv > nv) {
        std::cerr << "c Warning (" << filename << "), header is: p cnf " << nv << " " << nc
            << ", but correct would be: p cnf " << rv << " " << rc << std::endl;
    }
}

#endif  // SRC_NORMALIZE_H_
