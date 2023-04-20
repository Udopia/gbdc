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

#ifndef SRC_UTIL_GBDHASH_H_
#define SRC_UTIL_GBDHASH_H_

#include <string>
#include <sstream>

#include "lib/md5/md5.h"

#include "src/util/StreamBuffer.h"

std::string gbd_hash_from_dimacs(const char* filename) {
    MD5 md5;
    StreamBuffer in(filename);
    std::string clause("");
    while (!in.eof()) {
        in.skipWhitespace();
        if (in.eof()) {
            break;
        }
        if (*in == 'p' || *in == 'c') {
            in.skipLine();
        } else {
            for (int plit = in.readInteger(); plit != 0; plit = in.readInteger()) {
                clause.append(std::to_string(plit));
                clause.append(" ");
            }
            clause.append("0");
            //std::cout << "Hashing " << clause << std::endl;
            md5.consume(clause.c_str(), clause.length());
            clause.assign(" ");
        }
    }
    return md5.produce();
}


std::string opb_hash(const char* filename) {
    MD5 md5;
    StreamBuffer in(filename);
    std::stringstream norm;
    while (!in.eof()) {
        in.skipWhitespace();
        if (in.eof()) {
            break;
        }
        if (*in == '*') {
            in.skipLine();
        } 
        else if (*in == 'm') {
            in.skipString("min:");
            norm << "min:";
            for (; *in != ';'; in.skipWhitespace()) {
                long long weight = in.readLongLong();
                in.skipWhitespace();
                in.skipString("x");
                long long variable = in.readLongLong();
                norm << " " << weight << " x" << variable;
            }            
            norm << ";";
        }
        else {
            for (; *in != '>' && *in != '<' && *in != '='; in.skipWhitespace()) {
                long long weight = in.readLongLong();
                in.skipWhitespace();
                in.skipString("x");
                long long variable = in.readLongLong();
                norm << weight << " x" << variable << " ";
            }
            while (*in == '>' || *in == '<' || *in == '=') {
                norm << *in;
                ++in;
            }
            long long rhs = in.readLongLong();
            norm << " " << rhs << ";";
        }
        std::string str = norm.str();
        md5.consume(str.c_str(), str.length());
        norm.str(std::string());
        norm.clear();
        in.skipLine();
    }
    return md5.produce();
}


#endif  // SRC_UTIL_GBDHASH_H_
