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

#ifndef GBDHASH_H_
#define GBDHASH_H_

#include <string>
#include <sstream>

#include "src/external/md5/md5.h"
#include "src/util/StreamBuffer.h"

namespace CNF {
    std::string gbdhash(const char* filename) {
        MD5 md5;
        StreamBuffer in(filename);
        bool notfirst = false;
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                if (notfirst) md5.consume(" ", 1);
                std::string plit;
                while (in.readNumber(&plit)) {
                    if (plit == "0") break;
                    md5.consume(plit.c_str(), plit.length());
                    md5.consume(" ", 1);
                }
                md5.consume("0", 1);
                notfirst = true;
            }
        }
        return md5.produce();
    }
} // namespace CNF 

namespace PQBF {
    std::string gbdhash(const char* filename) {
        MD5 md5;
        StreamBuffer in(filename);
        bool notfirst = false;
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                if (notfirst) md5.consume(" ", 1);
                if (*in == 'e' || *in == 'a') {
                    md5.consume(*in == 'e' ? "e " : "a ", 2);
                    in.skip();
                    in.skipWhitespace();
                }
                std::string plit;
                while (in.readNumber(&plit)) {
                    if (plit == "0") break;
                    md5.consume(plit.c_str(), plit.length());
                    md5.consume(" ", 1);
                }
                md5.consume("0", 1);
                notfirst = true;
            }
        }
        return md5.produce();
    }
} // namespace PQBF

namespace OPB {
    std::string gbdhash(const char* filename) {
        MD5 md5;
        StreamBuffer in(filename);
        std::string num;
        while (in.skipWhitespace()) {
            if (*in == '*') {
                if (!in.skipLine()) break;
            } 
            else if (*in == 'm') {
                md5.consume("min:", 4);
                in.skipString("min:");
                in.skipWhitespace();
                while (*in != ';') {
                    if (*in == 'x') {
                        md5.consume(" x", 2);
                        in.skip();
                    } else if (*in == '~') {
                        md5.consume(" ~x", 3);
                        in.skip();
                        in.skipWhitespace();
                        in.skip();
                    } else {
                        md5.consume(" ", 1);
                    }
                    in.readNumber(&num);
                    md5.consume(num.c_str(), num.length());
                    in.skipWhitespace();
                }
                md5.consume(";", 1);
            }
            else {
                while (*in != '>' && *in != '<' && *in != '=') {
                    if (*in == 'x') {
                        md5.consume("x", 1);
                        in.skip();
                    } else if (*in == '~') {
                        md5.consume("~x", 2);
                        in.skip();
                        in.skipWhitespace();
                        in.skip();
                    }
                    in.readNumber(&num);
                    md5.consume(num.c_str(), num.length());
                    md5.consume(" ", 1);
                    in.skipWhitespace();
                }
                while (*in == '>' || *in == '<' || *in == '=') {
                    const char c = *in;
                    md5.consume(&c, 1);
                    in.skip();
                }
                in.readNumber(&num);
                md5.consume(" ", 1);
                md5.consume(num.c_str(), num.length());
                md5.consume(";", 1);
                in.skipWhitespace();
            }
            if (*in == ';') in.skip();
        }
        return md5.produce();
    }
} // namespace OPB

namespace WCNF {
    std::string gbdhash(const char* filename) {
        MD5 md5;
        StreamBuffer in(filename);
        uint64_t top = 0; // if top is 0, parsing new file format
        bool notfirst = false;
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
                assert(top == 0);  // should not have top in new format
                in.skip();
                if (notfirst) md5.consume(" ", 1);
                md5.consume("h ", 2);
                std::string plit;
                while (in.readNumber(&plit)) {
                    if (plit == "0") break;
                    md5.consume(plit.c_str(), plit.length());
                    md5.consume(" ", 1);
                }
                md5.consume("0", 1);
            } else {
                if (notfirst) md5.consume(" ", 1);
                if (top > 0) {
                    // parse old format clause
                    uint64_t nbr;
                    in.readUInt64(&nbr);
                    if (nbr >= top) {
                        // hard clause
                        md5.consume("h ", 2);
                    } else {
                        // soft clause
                        std::string weight = std::to_string(nbr);
                        md5.consume(weight.c_str(), weight.length());
                        md5.consume(" ", 1);
                    }
                }
                std::string plit;
                while (in.readNumber(&plit)) {
                    if (plit == "0") break;
                    md5.consume(plit.c_str(), plit.length());
                    md5.consume(" ", 1);
                }
                md5.consume("0", 1);
                notfirst = true;
            }
        }
        return md5.produce();
    }
} // namespace WCNF

#endif  // GBDHASH_H_
