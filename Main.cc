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

#include <iostream>
#include "src/GBDHash.h"
#include "src/Normalize.h"
#include "lib/cxxopts/cxxopts.hpp"

int main(int argc, char** argv) {
    cxxopts::Options options("cnftools", "Toolbox for DIMACS CNF");

    options.add_options()
    ("t,tool", "Select Tool: gbdhash, normalize", cxxopts::value<std::string>()->default_value("gbdhash"))
    ("f,file", "Filename", cxxopts::value<std::string>(), "Filename")
    ("h,help", "Usage Info", cxxopts::value<bool>()->default_value("false"))
    ;

    options.parse_positional("file");
    options.positional_help("Filename");
    
    cxxopts::ParseResult result;
    try {
        result = options.parse(argc, argv);
    } 
    catch (cxxopts::OptionParseException e) {
        std::cout << e.what() << std::endl;
    }

    if (result.count("help") || !result.count("file")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    std::string filename = result["file"].as<std::string>();
    std::string toolname = result["tool"].as<std::string>();

    if (toolname == "gbdhash") {
        std::cout << gbd_hash_from_dimacs(filename.c_str()) << std::endl;
    }
    else if (toolname == "normalize") {
        std::cerr << "normalizing " << filename << std::endl;
        normalize(filename.c_str());
    }

    return 0;
}