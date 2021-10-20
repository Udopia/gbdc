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
#include "lib/ipasir.h"

#include "src/util/GBDHash.h"
#include "src/util/CNFFormula.h"
#include "src/util/SolverTypes.h"

#include "src/transform/IndependentSet.h"
#include "src/transform/Normalize.h"

#include "src/features/GateStats.h"
#include "src/features/CNFStats.h"


int main(int argc, char** argv) {
    argparse::ArgumentParser argparse("CNF Tools");

    argparse.add_argument("tool").help("Select Tool: solve, gbdhash, normalize, isp, extract, gates")
        .default_value("gbdhash")
        .action([](const std::string& value) {
            static const std::vector<std::string> choices = { "solve", "gbdhash", "normalize", "isp", "extract", "gates" };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            return std::string{ "gbdhash" };
        });

    argparse.add_argument("file").help("Give Path");

    argparse.add_argument("-t", "--timeout")
        .help("Timeout in seconds (default: 0, disabled)")
        .default_value(0)
        .scan<'i', int>();

    argparse.add_argument("-m", "--memout")
        .help("Memout in megabytes (default: 0, disabled)")
        .default_value(0)
        .scan<'i', int>();

    argparse.add_argument("-v", "--verbose")
        .help("Verbosity level (default: 0, disabled)")
        .default_value(0)
        .scan<'i', int>();

    argparse.add_argument("-r", "--repeat")
        .help("Give number of root selections for gate recognition")
        .default_value(1)
        .scan<'i', int>();

    try {
        argparse.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cout << err.what() << std::endl;
        std::cout << argparse;
        exit(0);
    }

    std::string filename = argparse.get("file");
    std::string toolname = argparse.get("tool");
    int repeat = argparse.get<int>("repeat");
    ResourceLimits limits(argparse.get<int>("timeout"), argparse.get<int>("memout"));
    int verbose = argparse.get<int>("verbose");

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
        CNFFormula formula;
        formula.readDimacsFromFile(filename.c_str());

        CNFStats stats(formula, limits);
        stats.analyze();
        std::vector<float> record = stats.BaseFeatures();
        std::vector<std::string> names = CNFStats::BaseFeatureNames();
        for (unsigned i = 0; i < record.size(); i++) {
            std::cout << names[i] << "=" << record[i] << std::endl;
        }
    } else if (toolname == "gates") {
        CNFFormula formula;
        formula.readDimacsFromFile(filename.c_str());
        std::cout << "Finished Reading " << std::endl;
        GateStats stats(formula, limits);
        stats.analyze(repeat, verbose);
        std::vector<float> record = stats.GateFeatures();
        std::vector<std::string> names = GateStats::GateFeatureNames();
        for (unsigned i = 0; i < record.size(); i++) {
            std::cout << names[i] << "=" << record[i] << std::endl;
        }
    }

    return 0;
}
