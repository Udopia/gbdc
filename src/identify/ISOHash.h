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

#ifndef ISOHASH_H_
#define ISOHASH_H_

#include <vector>
#include <algorithm>
#include <stdio.h>

#include "lib/md5/md5.h"

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"


namespace CNF {
    /**
     * @brief Hashsum of ordered degree sequence of literal incidence graph
     * - literal nodes are grouped pairwise and sorted lexicographically
     * - edge weight of 1/n ==> literal node degree = occurence count
     * @param filename benchmark instance
     * @return std::string isohash
     */
    std::string isohash(const char* filename) {
        StreamBuffer in(filename);
        struct Node { unsigned neg; unsigned pos; };
        std::vector<Node> degrees;
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                int plit;
                while (in.readInteger(&plit)) {
                    if (abs(plit) >= degrees.size()) degrees.resize(abs(plit));
                    if (plit == 0) break;
                    else if (plit < 0) ++degrees[abs(plit) - 1].neg;
                    else ++degrees[abs(plit) - 1].pos;
                }
            }
        }
        // get invariant w.r.t. polarity flips
        for (Node& degree : degrees) {
            if (degree.pos < degree.neg) std::swap(degree.pos, degree.neg);
        }
        // sort lexicographically by degree
        std::sort(degrees.begin(), degrees.end(), [](const Node& one, const Node& two) { 
            return one.neg != two.neg ? one.neg < two.neg : one.pos < two.pos; 
        } );
        // hash
        MD5 md5;
        char buffer[64];
        for (Node node : degrees) {
            int n = snprintf(buffer, sizeof(buffer), "%u %u ", node.neg, node.pos);
            md5.consume(buffer, n);
        }
        return md5.produce();
    }
} // namespace CNF

#endif  // ISOHASH_H_
