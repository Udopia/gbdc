/*************************************************************************************************
CNFTools -- Copyright (c) 2020, Ashlin Iser, KIT - Karlsruhe Institute of Technology

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


#include "src/external/argparse/argparse.h"
#include "src/external/ipasir.h"

#include "src/identify/GBDHash.h"
#include "src/identify/ISOHash.h"
#include "src/identify/ISOHash2.h"

#include "src/util/SolverTypes.h"

#include "src/util/ResourceLimits.h"
#include "src/transform/cnf2bip.h"
#include "src/transform/cnf2kis.h"
#include "src/transform/cnf2cnf.h"

#include "src/extract/CNFSaniCheck.h"
#include "src/extract/CNFBaseFeatures.h"
#include "src/extract/CNFGateFeatures.h"
#include "src/extract/WCNFBaseFeatures.h"
#include "src/extract/OPBBaseFeatures.h"

#include "src/util/StreamCompressor.h"

int main(int argc, char** argv) {
    argparse::ArgumentParser argparse("CNF Tools");

    argparse.add_argument("tool").help("Select Tool: id, isohash, isohash2, normalize, sanitize, checksani, cnf2kis, cnf2bip, extract, gates")
        .default_value("identify")
        .action([](const std::string& value) {
            static const std::vector<std::string> choices = { "id", "isohash", "isohash2", "normalize", "sanitize", "checksani", "cnf2kis", "cnf2bip", "extract", "gates" };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            return std::string{ "identify" };
        });

    argparse.add_argument("file").help("Path to Input File");
    argparse.add_argument("-o", "--output").default_value(std::string("-")).help("Path to Output File (used by cnf2* transformers, default is stdout)");
    argparse.add_argument("-t", "--timeout").default_value(0).scan<'i', int>().help("Time limit in seconds");
    argparse.add_argument("-m", "--memout").default_value(0).scan<'i', int>().help("Memory limit in MB");
    argparse.add_argument("-f", "--fileout").default_value(0).scan<'i', int>().help("File size limit in MB"); 
    argparse.add_argument("--max-iters").default_value(6).scan<'i', int>().help("Maximum Isohash2 iterations before stopping");

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

    ResourceLimits limits(argparse.get<int>("timeout"), argparse.get<int>("memout"), argparse.get<int>("fileout"));
    limits.set_rlimits();
    std::cerr << "c Running: " << toolname << " " << filename << std::endl;

    std::string ext = std::filesystem::path(filename).extension();
    if (ext == ".xz" || ext == ".lzma" || ext == ".bz2" || ext == ".gz") {
        ext = std::filesystem::path(filename).stem().extension();
    }

    if (ext == ".cnf" || ext == ".wecnf") {
        std::cerr << "Detected CNF" << std::endl;
    } else if (ext == ".opb") {
        std::cerr << "Detected OPB" << std::endl;
    } else if (ext == ".qcnf" || ext == ".qdimacs") {
        std::cerr << "Detected QBF" << std::endl;
    } else if (ext == ".wcnf") {
        std::cerr << "Detected WCNF" << std::endl;
    }

    try {
        if (toolname == "id") {
            if (ext == ".cnf" || ext == ".wecnf") {
                std::cout << CNF::gbdhash(filename.c_str()) << std::endl;
            }
            else if (ext == ".opb") {
                std::cout << OPB::gbdhash(filename.c_str()) << std::endl;
            }
            else if (ext == ".qcnf" || ext == ".qdimacs") {
                std::cout << PQBF::gbdhash(filename.c_str()) << std::endl;
            }
            else if (ext == ".wcnf") {
                std::cout << WCNF::gbdhash(filename.c_str()) << std::endl;
            }
        } else if (toolname == "isohash") {
            if (ext == ".cnf") {
                std::cout << CNF::isohash(filename.c_str()) << std::endl;
            }
            else if (ext == ".wcnf") {
                std::cout << WCNF::isohash(filename.c_str()) << std::endl;
            }
        } 
        else if (toolname == "isohash2") {
            if (ext == ".cnf") {
                CNF::IsoHash2Settings config;
                config.max_iterations = argparse.get<int>("--max-iters");
                std::cout << CNF::isohash2(filename.c_str(), config) << std::endl;
            }
        } 
        else if (toolname == "normalize") {
            std::cerr << "Normalizing " << filename << std::endl;
            CNF::Normaliser norm(filename.c_str(), output == "-" ? nullptr : output.c_str());
            norm.run();
        } 
        else if (toolname == "sanitize") {
            CNF::Sanitiser sani(filename.c_str(), output == "-" ? nullptr : output.c_str());
            sani.run();
        }
        else if (toolname == "checksani") {
            CNF::SaniCheck ana(filename.c_str(), true);
            ana.run();
            bool sani = true;
            std::cout << "hash" << " " << CNF::gbdhash(filename.c_str()) << std::endl;
            std::cout << "filename " << filename << std::endl;
            sani = ana.getFeature("head_vars") == ana.getFeature("norm_vars") && ana.getFeature("head_clauses") == ana.getFeature("norm_clauses");
            // std::cout << ana.getFeature("head_vars") << " " << ana.getFeature("norm_vars") << std::endl;
            // std::cout << ana.getFeature("head_clauses") << " " << ana.getFeature("norm_clauses") << std::endl;
            std::cout << "header_consistent " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("whitespace_normalised") == 1.0;
            std::cout << "whitespace_normalised " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("has_comment") == 0.0;
            std::cout << "no_comment " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("has_tautological_clause") == 0.0;
            std::cout << "no_tautological_clause " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("has_duplicate_literals") == 0.0;
            std::cout << "no_duplicate_literals " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("has_empty_clause") == 0.0;
            std::cout << "no_empty_clause " << (sani ? "yes" : "no") << std::endl;
        }
        else if (toolname == "cnf2kis") {
            std::cerr << "Generating Independent Set Problem " << filename << std::endl;
            IndependentSetFromCNF gen(filename.c_str());
            gen.generate_independent_set_problem(output == "-" ? nullptr : output.c_str());
        }
        else if (toolname == "cnf2bip") {
            std::cerr << "Generating Bipartite Graph " << filename << std::endl;
            CNF::cnf2bip gen(filename.c_str(), output == "-" ? nullptr : output.c_str());
            gen.run();
        }
        else if (toolname == "extract" || toolname == "gates") {
            IExtractor *stats;
            if (toolname == "extract") {
                if (ext == ".cnf") {
                    stats = new CNF::BaseFeatures(filename.c_str());
                }
                else if (ext == ".wcnf") {
                    stats = new WCNF::BaseFeatures(filename.c_str());
                }
                else if (ext == ".opb") {
                    stats = new OPB::BaseFeatures(filename.c_str());
                }
                else {
                    std::cerr << "Format " << ext << " not supported by extract" << std::endl;
                    return 1;
                }
            }
            else if (toolname == "gates") {
                stats = new CNF::GateFeatures(filename.c_str());
                if (ext != ".cnf") {
                    std::cerr << "Format " << ext << " not supported by extract" << std::endl;
                    return 1;
                }
            }
            stats->run();
            for (std::string name : stats->getNames()) {
                std::cout << name << " ";
            }
            std::cout << std::endl;
            for (double feature : stats->getFeatures()) {
                std::cout << feature << " ";
            }
            std::cout << std::endl;
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
