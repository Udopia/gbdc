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
#include <array>

#include "src/ipasir.h"
#include "src/GBDHash.h"
#include "src/Normalize.h"
#include "src/util/CNFFormula.h"
#include "src/util/CNFStats.h"
#include "src/util/SolverTypes.h"
#include "src/util/Runtime.h"
#include "src/gates/GateAnalyzer.h"
#include "src/IndependentSet.h"
#include "lib/cxxopts/cxxopts.hpp"


int main(int argc, char** argv) {
    cxxopts::Options options("cnftools", "Toolbox for DIMACS CNF");

    options.add_options()
    ("t,tool", "Select Tool: gbdhash, normalize, gates, isp", cxxopts::value<std::string>()->default_value("gbdhash"))
    ("f,file", "Filename", cxxopts::value<std::string>(), "Filename")
    ("h,help", "Usage Info", cxxopts::value<bool>()->default_value("false"))
// "Gates" Tool Options:
    ("p,gates-patterns", "Patterns", cxxopts::value<bool>()->default_value("true"))
    ("s,gates-semantic", "Semantic", cxxopts::value<bool>()->default_value("true"))
    ("r,gates-repeat", "Nof Root Selections", cxxopts::value<unsigned>()->default_value("1"), "Nof Root Selections");

    options.parse_positional("file");
    options.positional_help("Filename");

    cxxopts::ParseResult result;
    try {
        result = options.parse(argc, argv);
    } catch (cxxopts::OptionParseException& e) {
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
    } else if (toolname == "normalize") {
        std::cerr << "Normalizing " << filename << std::endl;
        normalize(filename.c_str());
    } else if (toolname == "isp") {
        std::cerr << "Generating Independent Set Problem " << filename << std::endl;
        generate_independent_set_problem(filename);
    } else if (toolname == "extract") {
        std::cerr << "Extracing features of " << filename << std::endl;
        CNFFormula F;
        F.readDimacsFromFile(filename.c_str());
        Runtime runtime;
        runtime.start();
        CNFStats stats(F);
        stats.analyze();
        runtime.stop();
        std::vector<float> record = stats.SatzillaFeatures();
        std::cout << "clauses="                 << record[0] << std::endl;
        std::cout << "variables="               << record[1] << std::endl;
        std::cout << "vcg_vdegrees_mean="       << record[2] << std::endl;
        std::cout << "vcg_vdegrees_variance="   << record[3] << std::endl;
        std::cout << "vcg_vdegrees_min="        << record[4] << std::endl;
        std::cout << "vcg_vdegrees_max="        << record[5] << std::endl;
        std::cout << "vcg_vdegrees_entropy="    << record[6] << std::endl;
        std::cout << "vcg_cdegrees_mean="       << record[7] << std::endl;
        std::cout << "vcg_cdegrees_variance="   << record[8] << std::endl;
        std::cout << "vcg_cdegrees_min="        << record[9] << std::endl;
        std::cout << "vcg_cdegrees_max="        << record[10] << std::endl;
        std::cout << "vcg_cdegrees_entropy="    << record[11] << std::endl;
        std::cout << "vg_degrees_mean="         << record[12] << std::endl;
        std::cout << "vg_degrees_variance="     << record[13] << std::endl;
        std::cout << "vg_degrees_min="          << record[14] << std::endl;
        std::cout << "vg_degrees_max="          << record[15] << std::endl;
        std::cout << "vg_degrees_entropy="      << record[16] << std::endl;
        std::cout << "vg_jwdegrees_mean="       << record[17] << std::endl;
        std::cout << "vg_jwdegrees_variance="   << record[18] << std::endl;
        std::cout << "vg_jwdegrees_min="        << record[19] << std::endl;
        std::cout << "vg_jwdegrees_max="        << record[20] << std::endl;
        std::cout << "vg_jwdegrees_entropy="    << record[21] << std::endl;
        std::cout << "cg_degrees_mean="         << record[22] << std::endl;
        std::cout << "cg_degrees_variance="     << record[23] << std::endl;
        std::cout << "cg_degrees_min="          << record[24] << std::endl;
        std::cout << "cg_degrees_max="          << record[25] << std::endl;
        std::cout << "cg_degrees_entropy="      << record[26] << std::endl;
        std::cout << "balance_clause_mean="     << record[27] << std::endl;
        std::cout << "balance_clause_variance=" << record[28] << std::endl;
        std::cout << "balance_clause_min="      << record[29] << std::endl;
        std::cout << "balance_clause_max="      << record[30] << std::endl;
        std::cout << "balance_clause_entropy="  << record[31] << std::endl;
        std::cout << "balance_vars_mean="       << record[32] << std::endl;
        std::cout << "balance_vars_variance="   << record[33] << std::endl;
        std::cout << "balance_vars_min="        << record[34] << std::endl;
        std::cout << "balance_vars_max="        << record[35] << std::endl;
        std::cout << "balance_vars_entropy="    << record[36] << std::endl;
        std::cout << "clause_size_1="           << record[37] << std::endl;
        std::cout << "clause_size_2="           << record[38] << std::endl;
        std::cout << "clause_size_3="           << record[39] << std::endl;
        std::cout << "clause_size_4="           << record[40] << std::endl;
        std::cout << "clause_size_5="           << record[41] << std::endl;
        std::cout << "clause_size_6="           << record[42] << std::endl;
        std::cout << "clause_size_7="           << record[43] << std::endl;
        std::cout << "clause_size_8="           << record[44] << std::endl;
        std::cout << "clause_size_9="           << record[45] << std::endl;
        std::cout << "horn_clauses="            << record[46] << std::endl;
        std::cout << "horn_vars_mean="          << record[47] << std::endl;
        std::cout << "horn_vars_variance="      << record[48] << std::endl;
        std::cout << "horn_vars_min="           << record[49] << std::endl;
        std::cout << "horn_vars_max="           << record[50] << std::endl;
        std::cout << "horn_vars_entropy="       << record[51] << std::endl;
        std::cout << "inv_horn_clauses="        << record[52] << std::endl;
        std::cout << "inv_horn_vars_mean="      << record[53] << std::endl;
        std::cout << "inv_horn_vars_variance="  << record[54] << std::endl;
        std::cout << "inv_horn_vars_min="       << record[55] << std::endl;
        std::cout << "inv_horn_vars_max="       << record[56] << std::endl;
        std::cout << "inv_horn_vars_entropy="   << record[57] << std::endl;
        std::cout << "feature_extraction_time=" << runtime.get() << std::endl;
    } else if (toolname == "gates") {
        bool patterns = result["gates-patterns"].as<bool>();
        bool semantic = result["gates-semantic"].as<bool>();
        unsigned repeat = result["gates-repeat"].as<unsigned>();

        CNFFormula F;
        F.readDimacsFromFile(filename.c_str());
        std::cerr << "c Read instance of " << F.nVars() << " variables and " << F.nClauses() << " clauses" << std::endl;

        std::cerr << "c Recognizing Gates " << filename << std::endl;
        GateAnalyzer<> A(F, patterns, semantic, repeat);
        A.analyze();
        A.getGateFormula().printGates();
    } else if (toolname == "solve") {
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
