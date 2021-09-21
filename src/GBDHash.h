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

#ifndef SRC_GBDHASH_H_
#define SRC_GBDHASH_H_

#include <string>

#include "src/StreamBuffer.h"
#include "lib/md5/md5.h"

std::string gbd_hash_from_dimacs(const char* filename) {
    unsigned char sig[MD5_SIZE];
    char str[MD5_STRING_SIZE];
    md5::md5_t md5;
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
            md5.process(clause.c_str(), clause.length());
            clause.assign(" ");
        }
    }
    md5.finish(sig);
    md5::sig_to_string(sig, str, sizeof(str));
    return std::string(str);
}

#endif  // SRC_GBDHASH_H_
