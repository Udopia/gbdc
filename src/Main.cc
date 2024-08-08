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
#include <cstdio>
#include <filesystem>

#include "lib/argparse/argparse.hpp"
#include "lib/ipasir.h"

#include "src/identify/GBDHash.h"
#include "src/identify/ISOHash.h"

#include "src/util/SolverTypes.h"

#include "src/transform/IndependentSet.h"
#include "src/transform/Normalize.h"
#include "src/util/ResourceLimits.h"
#include "src/transform/cnf2bip.h"

#include "src/extract/CNFGateFeatures.h"
#include "src/extract/CNFBaseFeatures.h"
#include "src/extract/WCNFBaseFeatures.h"
#include "src/extract/OPBBaseFeatures.h"

#include "src/util/StreamCompressor.h"

int main(int argc, char** argv) {
    argparse::ArgumentParser argparse("CNF Tools");

    argparse.add_argument("tool").help("Select Tool: solve, id|identify (gbdhash, opbhash, pqbfhash), isohash, normalize, sanitize, checksani, cnf2kis, cnf2bip, extract, gates")
        .default_value("identify")
        .action([](const std::string& value) {
            static const std::vector<std::string> choices = { "solve", "id", "identify", "gbdhash", "opbhash", "pqbfhash", "isohash", "normalize", "sanitize", "checksani", "cnf2kis", "cnf2bip", "extract", "gates", "test" };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            return std::string{ "identify" };
        });

    argparse.add_argument("file").help("Path to Input File");
    argparse.add_argument("-o", "--output").help("Path to Output File (used by cnf2kis if set, default is stdout)").default_value(std::string("-"));

    argparse.add_argument("-t", "--timeout")
        .help("Timeout in seconds (default: 0, disabled)")
        .default_value(0)
        .scan<'i', int>();

    argparse.add_argument("-m", "--memout")
        .help("Memout in megabytes (default: 0, disabled)")
        .default_value(0)
        .scan<'i', int>();

    argparse.add_argument("-f", "--fileout")
        .help("Maximum generated file size in megabytes (default: 0, disabled)")
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
    std::string output = argparse.get("output");
    int verbose = argparse.get<int>("verbose");
    int repeat = argparse.get<int>("repeat");

    ResourceLimits limits(argparse.get<int>("timeout"), argparse.get<int>("memout"), argparse.get<int>("fileout"));
    limits.set_rlimits();
    std::cerr << "c Running: " << toolname << " " << filename << std::endl;

    try {
        if (toolname == "id" || toolname == "identify") {
            std::string ext = std::filesystem::path(filename).extension();
            if (ext == ".xz" || ext == ".lzma" || ext == ".bz2" || ext == ".gz") {
                ext = std::filesystem::path(filename).stem().extension();
            }
            if (ext == ".cnf" || ext == ".wecnf") {
                std::cerr << "Detected CNF, using CNF hash" << std::endl;
                std::cout << CNF::gbdhash(filename.c_str()) << std::endl;
            }
            else if (ext == ".opb") {
                std::cerr << "Detected OPB, using OPB hash" << std::endl;
                std::cout << OPB::gbdhash(filename.c_str()) << std::endl;
            }
            else if (ext == ".qcnf" || ext == ".qdimacs") {
                std::cerr << "Detected QBF, using QBF hash" << std::endl;
                std::cout << PQBF::gbdhash(filename.c_str()) << std::endl;
            }
            else if (ext == ".wcnf") {
                std::cerr << "Detected WCNF, using WCNF hash" << std::endl;
                std::cout << WCNF::gbdhash(filename.c_str()) << std::endl;
            }
        } else if (toolname == "gbdhash") {
            std::cout << CNF::gbdhash(filename.c_str()) << std::endl;
        } else if (toolname == "isohash") {
            std::string ext = std::filesystem::path(filename).extension();
            if (ext == ".xz" || ext == ".lzma" || ext == ".bz2" || ext == ".gz") {
                ext = std::filesystem::path(filename).stem().extension();
            }
            if (ext == ".cnf") {
                std::cerr << "Detected CNF, using CNF isohash" << std::endl;
                std::cout << CNF::isohash(filename.c_str()) << std::endl;
            } else if (ext == ".wcnf") {
                std::cerr << "Detected WCNF, using WCNF isohash" << std::endl;
                std::cout << WCNF::isohash(filename.c_str()) << std::endl;
            }
        } else if (toolname == "opbhash") {
            std::cout << OPB::gbdhash(filename.c_str()) << std::endl;
        } else if (toolname == "pqbfhash") {
            std::cout << PQBF::gbdhash(filename.c_str()) << std::endl;
        } else if (toolname == "normalize") {
            std::cerr << "Normalizing " << filename << std::endl;
            normalize(filename.c_str());
        } else if (toolname == "checksani") {
            if (!check_sanitized(filename.c_str())) {
                std::cerr << filename << " needs sanitization" << std::endl;
            }
        } else if (toolname == "sanitize") {
            sanitize(filename.c_str());
        } else if (toolname == "cnf2kis") {
            std::cerr << "Generating Independent Set Problem " << filename << std::endl;
            IndependentSetFromCNF gen(filename.c_str());
            gen.generate_independent_set_problem(output == "-" ? nullptr : output.c_str());
        } else if (toolname == "cnf2bip") {
            std::cerr << "Generating Bipartite Graph " << filename << std::endl;
            BipartiteGraphFromCNF gen(filename.c_str());
            gen.generate_bipartite_graph(output == "-" ? nullptr : output.c_str());
        } else if (toolname == "extract") {
            std::string ext = std::filesystem::path(filename).extension();
            if (ext == ".xz" || ext == ".lzma" || ext == ".bz2" || ext == ".gz") {
                ext = std::filesystem::path(filename).stem().extension();
            }
            if (ext == ".cnf") {
                std::cerr << "Detected CNF, extracting CNF base features" << std::endl;
                CNF::BaseFeatures stats(filename.c_str());
                stats.extract();
                std::vector<double> record = stats.getFeatures();
                std::vector<std::string> names = stats.getNames();
                for (unsigned i = 0; i < record.size(); i++) {
                    std::cout << names[i] << "=" << record[i] << std::endl;
                }
            } else if (ext == ".wcnf") {
                std::cerr << "Detected WCNF, extracting WCNF base features" << std::endl;
                WCNF::BaseFeatures stats(filename.c_str());
                stats.extract();
                std::vector<double> record = stats.getFeatures();
                std::vector<std::string> names = stats.getNames();
                for (unsigned i = 0; i < record.size(); i++) {
                    std::cout << names[i] << "=" << record[i] << std::endl;
                }
            } else if (ext == ".opb") {
                std::cerr << "Detected OPB, extracting OPB base features" << std::endl;
                OPB::BaseFeatures stats(filename.c_str());
                stats.extract();
                std::vector<double> record = stats.getFeatures();
                std::vector<std::string> names = stats.getNames();
                for (unsigned i = 0; i < record.size(); i++) {
                    std::cout << names[i] << "=" << record[i] << std::endl;
                }
            }
        } else if (toolname == "gates") {
            CNFGateFeatures stats(filename.c_str());
            stats.extract();
            std::vector<double> record = stats.getFeatures();
            std::vector<std::string> names = stats.getNames();
            for (unsigned i = 0; i < record.size(); i++) {
                std::cout << names[i] << "=" << record[i] << std::endl;
            }
        } else if (toolname == "test") {
            std::cout << "Testing something ... " << std::endl;
            StreamCompressor cmpr(filename.c_str(), 100);
            for (int i = 0; i < 10; i++)
                cmpr.write("0123456789", 10);
        }
    }
    catch (std::bad_alloc& e) {
        std::cerr << "Memory Limit Exceeded" << std::endl;
        return 1;
    }
    catch (MemoryLimitExceeded& e) {
        std::cerr << "Memory Limit Exceeded" << std::endl;
        return 1;
    }
    catch (TimeLimitExceeded& e) {
        std::cerr << "Time Limit Exceeded" << std::endl;
        return 1;
    }
    catch (FileSizeLimitExceeded& e) {
        std::remove(output.c_str());
        std::cerr << "File Size Limit Exceeded" << std::endl;
        return 1;
    }
    return 0;
}
