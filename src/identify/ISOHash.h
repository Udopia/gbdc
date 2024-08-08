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
                    if (static_cast<size_t>(abs(plit)) >= degrees.size()) degrees.resize(abs(plit));
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
            if (node.neg == 0 && node.pos == 0) continue;  // get invariant against variable gaps
            int n = snprintf(buffer, sizeof(buffer), "%u %u ", node.neg, node.pos);
            md5.consume(buffer, n);
        }
        return md5.produce();
    }
} // namespace CNF

namespace WCNF {
    std::string isohash(const char* filename) {
        StreamBuffer in(filename);
        struct Node { uint64_t neg; uint64_t pos; };
        std::vector<Node> hard_degrees;
        std::vector<Node> soft_degrees;
        uint64_t top = 0; // if top is 0, parsing new file format
        while (in.skipWhitespace()) {
            if (*in == 'c') {
                if (!in.skipLine()) break;
            } else if (*in == 'p') {
                // old format: extract top
                in.skip();
                in.skipWhitespace();
                in.skipString("wcnf");
                // skip vars
                in.skipNumber();
                // skip clauses
                in.skipNumber();
                // extract top
                in.readUInt64(&top);
                in.skipLine();
            } else if (*in == 'h') {
                assert(top == 0); // should not have top in new format
                in.skip();
                int plit;
                while (in.readInteger(&plit)) {
                    if (static_cast<size_t>(abs(plit)) >= hard_degrees.size()) hard_degrees.resize(abs(plit));
                    if (plit == 0) break;
                    else if (plit < 0) ++hard_degrees[abs(plit) - 1].neg;
                    else ++hard_degrees[abs(plit) - 1].pos;
                }
            } else {
                uint64_t weight;
                in.readUInt64(&weight);
                if (top != 0 && weight >= top) {
                    // old format hard clause
                    int plit;
                    while (in.readInteger(&plit)) {
                        if (static_cast<size_t>(abs(plit)) >= hard_degrees.size()) hard_degrees.resize(abs(plit));
                        if (plit == 0) break;
                        else if (plit < 0) ++hard_degrees[abs(plit) - 1].neg;
                        else ++hard_degrees[abs(plit) - 1].pos;
                    }
                } else {
                    // soft clause
                    int plit;
                    while (in.readInteger(&plit)) {
                        if (static_cast<size_t>(abs(plit)) >= soft_degrees.size()) soft_degrees.resize(abs(plit));
                        if (plit == 0) break;
                        else if (plit < 0) ++soft_degrees[abs(plit) - 1].neg += weight;
                        else ++soft_degrees[abs(plit) - 1].pos += weight;
                    }
                }
            }
        }
        // add hard degrees to soft degrees
        if (soft_degrees.size() < hard_degrees.size()) soft_degrees.resize(hard_degrees.size());
        std::transform(hard_degrees.begin(), hard_degrees.end(), soft_degrees.begin(), soft_degrees.begin(), [](Node a, Node b) { return Node { .neg = a.neg + b.neg, .pos = a.pos + b.pos }; });
        // get invariant w.r.t. polarity flips
        for (Node& degree : hard_degrees) {
            if (degree.pos < degree.neg) std::swap(degree.pos, degree.neg);
        }
        for (Node& degree : soft_degrees) {
            if (degree.pos < degree.neg) std::swap(degree.pos, degree.neg);
        }
        // sort lexicographically by degree
        auto lex_smaller = [](const Node& one, const Node& two) { return one.neg != two.neg ? one.neg < two.neg : one.pos < two.pos; };
        std::sort(hard_degrees.begin(), hard_degrees.end(), lex_smaller);
        std::sort(soft_degrees.begin(), soft_degrees.end(), lex_smaller);
        // hash
        MD5 md5;
        char buffer[64];
        for (Node node : hard_degrees) {
            if (node.neg == 0 && node.pos == 0) continue;  // get invariant against variable gaps
            int n = snprintf(buffer, sizeof(buffer), "%lu %lu ", node.neg, node.pos);
            md5.consume(buffer, n);
        }
        md5.consume("softs ", 6);
        for (Node node : soft_degrees) {
            if (node.neg == 0 && node.pos == 0) continue;  // get invariant against variable gaps
            int n = snprintf(buffer, sizeof(buffer), "%lu %lu ", node.neg, node.pos);
            md5.consume(buffer, n);
        }
        return md5.produce();
    }
} // namespace CNF

#endif  // ISOHASH_H_
