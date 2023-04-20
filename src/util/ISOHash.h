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

#ifndef SRC_UTIL_ISOHASH_H_
#define SRC_UTIL_ISOHASH_H_

#include <vector>
#include <algorithm>
#include <stdio.h>

#include "lib/md5/md5.h"

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

std::string iso_hash_from_dimacs(const char* filename) {
    StreamBuffer in(filename);
    struct Node { unsigned in; unsigned out; };
    std::vector<Node> degrees;
    // count in- and out-degress per node in VIG
    while (!in.eof()) {
        in.skipWhitespace();
        if (in.eof()) {
            break;
        }
        if (*in == 'p' || *in == 'c') {
            in.skipLine();
        } else {
            for (int plit = in.readInteger(); plit != 0; plit = in.readInteger()) {
                if (2 * (unsigned)abs(plit) > degrees.size()) {
                    degrees.resize(2 * abs(plit) + 1);
                }
                if (plit < 0) {
                    ++degrees[abs(plit)].in;
                } else {
                    ++degrees[abs(plit)].out;
                }
            }
        }
    }
    // get invariant w.r.t. polarity flips
    for (Node& degree : degrees) {
        if (degree.out < degree.in) std::swap(degree.out, degree.in);
    }
    // sort lexicographically by degree
    std::sort(degrees.begin(), degrees.end(), [](const Node& one, const Node& two) { return one.in != two.in ? one.in < two.in : one.out < two.out; } );
    // hash
    MD5 md5;
    char buffer[64];
    for (Node node : degrees) {
        int n = snprintf(buffer, sizeof(buffer), "%u %u ", node.in, node.out);
        md5.consume(buffer, n);
    }
    return md5.produce();
}

#endif  // SRC_UTIL_ISOHASH_H_
