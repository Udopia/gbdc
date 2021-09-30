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

#include "lib/argparse/argparse.hpp"

#include "src/ipasir.h"
#include "src/GBDHash.h"
#include "src/Normalize.h"
#include "src/util/CNFFormula.h"
#include "src/util/SolverTypes.h"
#include "src/IndependentSet.h"

#include "src/features/GateStats.h"
#include "src/features/CNFStats.h"


int main(int argc, char** argv) {
    argparse::ArgumentParser program("CNF Tools");

    program.add_argument("tool").help("Select Tool: solve, gbdhash, normalize, isp, extract, gates")
        .default_value("gbdhash")
        .action([](const std::string& value) {
            static const std::vector<std::string> choices = { "solve", "gbdhash", "normalize", "isp", "extract", "gates" };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            return std::string{ "gbdhash" };
        });

    program.add_argument("file").help("Give Path");

    program.add_argument("-r", "--repeat")
        .help("Give number of root selections for gate recognition")
        .default_value(1)
        .scan<'i', int>();

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cout << err.what() << std::endl;
        std::cout << program;
        exit(0);
    }

    std::string filename = program.get("file");
    std::string toolname = program.get("tool");
    unsigned repeat = program.get<int>("repeat");

    std::cerr << "c Running: '" << toolname << " " << filename << std::endl;

    if (toolname == "gbdhash") {
        std::cout << gbd_hash_from_dimacs(filename.c_str()) << std::endl;
    } else if (toolname == "normalize") {
        std::cerr << "Normalizing " << filename << std::endl;
        normalize(filename.c_str());
    } else if (toolname == "isp") {
        std::cerr << "Generating Independent Set Problem " << filename << std::endl;
        generate_independent_set_problem(filename);
    } else if (toolname == "extract") {
        CNFFormula F;
        F.readDimacsFromFile(filename.c_str());

        CNFStats stats(F);
        stats.analyze();
        std::vector<float> record = stats.BaseFeatures();
        std::vector<std::string> names = CNFStats::BaseFeatureNames();
        for (unsigned i = 0; i < record.size(); i++) {
            std::cout << names[i] << "=" << record[i] << std::endl;
        }
    } else if (toolname == "gates") {
        CNFFormula F;
        F.readDimacsFromFile(filename.c_str());

        GateStats stats(F);
        stats.analyze(repeat);
        std::vector<float> record = stats.GateFeatures();
        std::vector<std::string> names = GateStats::GateFeatureNames();
        for (unsigned i = 0; i < record.size(); i++) {
            std::cout << names[i] << "=" << record[i] << std::endl;
        }
    } else if (toolname == "solve") {
        CNFFormula F;
        F.readDimacsFromFile(filename.c_str());

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
