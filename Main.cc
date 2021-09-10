/*************************************************************************************************
CNFTools -- Copyright (c) 2020, Markus Iser, KIT - Karlsruhe Institute of Technology

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
#include "src/ipasir.h"
#include "src/GBDHash.h"
#include "src/Normalize.h"
#include "src/util/CNFFormula.h"
#include "src/util/SolverTypes.h"
#include "src/util/Runtime.h"
#include "src/gates/GateAnalyzer.h"
#include "src/IndependentSet.h"
#include "lib/cxxopts/cxxopts.hpp"

#include <array>

int main(int argc, char** argv) {
    cxxopts::Options options("cnftools", "Toolbox for DIMACS CNF");

    options.add_options()
    ("t,tool", "Select Tool: gbdhash, normalize, gates, isp", cxxopts::value<std::string>()->default_value("gbdhash"))
    ("f,file", "Filename", cxxopts::value<std::string>(), "Filename")
    ("n,number", "Magic Number", cxxopts::value<unsigned>()->default_value("0"), "Magic Number")
    ("h,help", "Usage Info", cxxopts::value<bool>()->default_value("false"))
// "Gates" Tool Options:
    ("p,gates-patterns", "Patterns", cxxopts::value<bool>()->default_value("true"))
    ("s,gates-semantic", "Semantic", cxxopts::value<bool>()->default_value("true"))
    ("r,gates-repeat", "Repeated Root Selection", cxxopts::value<unsigned>()->default_value("1"), "Number of Root Selections")
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

    std::cerr << "Running " << toolname << " with " << filename << std::endl;

    if (toolname == "gbdhash") {
        std::cout << gbd_hash_from_dimacs(filename.c_str()) << std::endl;
    }
    else if (toolname == "normalize") {
        std::cerr << "Normalizing " << filename << std::endl;
        normalize(filename.c_str());
    }
    else if (toolname == "isp") {
        std::cerr << "Generating Independent Set Problem " << filename << std::endl;
        generate_independent_set_problem(filename);
    }
    else if (toolname == "gates") {
        bool patterns = result["gates-patterns"].as<bool>();
        bool semantic = result["gates-semantic"].as<bool>();
        unsigned repeat = result["gates-repeat"].as<unsigned>();
        
        CNFFormula F;
        F.readDimacsFromFile(filename.c_str());
        std::cerr << "c Read instance of " << F.nVars() << " variables and " << F.nClauses() << " clauses" << std::endl;
        
        std::cerr << "c Recognizing Gates " << filename << std::endl;
        GateAnalyzer<> A (F, patterns, semantic, repeat);
        A.analyze();
        A.getGateFormula().printGates();
    }
    else if (toolname == "gates2") {
        bool patterns = result["gates-patterns"].as<bool>();
        bool semantic = result["gates-semantic"].as<bool>();
        unsigned repeat = result["gates-repeat"].as<unsigned>();
        
        CNFFormula F;
        F.readDimacsFromFile(filename.c_str());
        std::cerr << "c Read instance of " << F.nVars() << " variables and " << F.nClauses() << " clauses" << std::endl;
        
        std::cerr << "c Recognizing Gates " << filename << std::endl;
        GateAnalyzer<BlockList> A (F, patterns, semantic, repeat);
        A.analyze();
        A.getGateFormula().printGates();
    }
    else if (toolname == "solve") {
        CNFFormula F;
        F.readDimacsFromFile(filename.c_str());
        std::cerr << "c Read instance of " << F.nVars() << " variables and " << F.nClauses() << " clauses" << std::endl;
        
        std::cerr << "c Intitializing Solver " << ipasir_signature() << std::endl;
        void* S = ipasir_init();
        for (Cl* clause : F) {
            for (Lit lit : *clause) {
                ipasir_add(S, lit.toDimacs());
            }
            ipasir_add(S, 0);
        }
        int res = ipasir_solve(S);
        std::cout << ((res == 10) ? "s SATISFIABLE" : "s UNSATISFIABLE") << std::endl;
        ipasir_release(S);
    }

    return 0;
}