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

/*
 * gbdc standalone tools.
 *
 * This binary implements the gbd external-tool contract (see gbd/doc/external-tools-contract.md).
 * A single binary is dispatched to a concrete tool either by its invocation name (argv[0], e.g.
 * "gbd-extract-base") or, when invoked as "gbdc"/"gbdctool", by the first positional argument
 * (legacy multi-tool mode: "gbdc base file").
 *
 * Machine-readable (--gbd) contract:
 *   - stdout carries a stream of "<feature> <value>" lines (one per line).
 *   - extractors emit the features of the input instance (the input hash is attached by gbd).
 *   - transformers emit features of the produced instance, including "local", "hash" and links.
 *   - two reserved lines convey the outcome: "status <success|timeout|memout>" and "runtime <sec>".
 *   - the produced instance of a transformer goes to -o (optionally xz-compressed) or, without -o,
 *     to stderr; stdout stays reserved for the metadata stream.
 *   - "--feature-names" prints "<feature> [default]" per line; a default marks a unique (1:1)
 *     feature, its absence marks a non-unique (1:n) feature.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "src/external/argparse/argparse.h"

#include "src/identify/GBDHash.h"
#include "src/identify/ISOHash.h"
#include "src/identify/ISOHash2.h"

#include "src/util/ResourceLimits.h"
#include "src/util/StreamCompressor.h"

#include "src/transform/cnf2bip.h"
#include "src/transform/cnf2kis.h"
#include "src/transform/cnf2cnf.h"

#include "src/extract/CNFSaniCheck.h"
#include "src/extract/CNFBaseFeatures.h"
#include "src/extract/CNFGateFeatures.h"
#include "src/extract/WCNFBaseFeatures.h"
#include "src/extract/OPBBaseFeatures.h"


namespace {

enum class Mode { HUMAN, GBD };

/* Extract the base name of a path. */
std::string basename_of(const std::string& path) {
    return std::filesystem::path(path).filename().string();
}

/* Map an invocation name to an internal tool id, following the "gbd-<tool>" convention.
 * Returns an empty string when the name does not encode a tool (e.g. plain "gbdc"). */
std::string tool_from_invocation(const std::string& argv0) {
    const std::string name = basename_of(argv0);
    static const std::vector<std::pair<std::string, std::string>> map = {
        {"gbd-extract-base", "base"},
        {"gbd-extract-gate", "gate"},
        {"gbd-extract-wcnf", "wcnfbase"},
        {"gbd-extract-opb", "opbbase"},
        {"gbd-checksani", "checksani"},
        {"gbd-isohash2", "isohash2"},
        {"gbd-isohash", "isohash"},
        {"gbd-identify", "identify"},
        {"gbd-cnf2kis", "cnf2kis"},
        {"gbd-cnf2bip", "cnf2bip"},
        {"gbd-sanitize", "sanitize"},
        {"gbd-normalize", "normalize"},
    };
    for (const auto& [invocation, tool] : map) {
        if (name == invocation) return tool;
    }
    return "";
}

/* Normalise a few legacy subcommand spellings to the canonical tool ids. */
std::string canonical_tool(const std::string& tool) {
    if (tool == "id" || tool == "hash") return "identify";
    if (tool == "extract") return "base";
    if (tool == "gates") return "gate";
    return tool;
}

/* Detect the input format from the file extension, ignoring a compression suffix. */
std::string detect_extension(const std::string& filename) {
    std::filesystem::path p(filename);
    std::string ext = p.extension().string();
    if (ext == ".xz" || ext == ".lzma" || ext == ".bz2" || ext == ".gz") {
        ext = p.stem().extension().string();
    }
    return ext;
}

/* Format a feature value: integral values are printed without a fractional part or scientific
 * notation so that large counts survive the round-trip through gbd. */
std::string format_value(double value) {
    if (std::isfinite(value) && std::floor(value) == value && std::fabs(value) < 9e15) {
        return std::to_string(static_cast<long long>(value));
    }
    std::ostringstream oss;
    oss.precision(std::numeric_limits<double>::max_digits10);
    oss << value;
    return oss.str();
}


/* --- Extractors ---------------------------------------------------------------------------- */

/* Instantiate the extractor matching the tool id and the input format. */
IExtractor* make_extractor(const std::string& tool, const std::string& ext, const std::string& filename) {
    if (tool == "base") {
        if (ext == ".cnf") return new CNF::BaseFeatures(filename.c_str());
        throw std::runtime_error("base extractor requires a .cnf file");
    }
    if (tool == "gate") {
        if (ext == ".cnf") return new CNF::GateFeatures(filename.c_str());
        throw std::runtime_error("gate extractor requires a .cnf file");
    }
    if (tool == "wcnfbase") {
        if (ext == ".wcnf") return new WCNF::BaseFeatures(filename.c_str());
        throw std::runtime_error("wcnf extractor requires a .wcnf file");
    }
    if (tool == "opbbase") {
        if (ext == ".opb") return new OPB::BaseFeatures(filename.c_str());
        throw std::runtime_error("opb extractor requires a .opb file");
    }
    throw std::runtime_error("unknown extractor: " + tool);
}

/* Names of the features produced by an extractor tool (used by --feature-names). */
std::vector<std::string> extractor_feature_names(const std::string& tool) {
    if (tool == "base") return CNF::BaseFeatures("").getNames();
    if (tool == "gate") return CNF::GateFeatures("").getNames();
    if (tool == "wcnfbase") return WCNF::BaseFeatures("").getNames();
    if (tool == "opbbase") return OPB::BaseFeatures("").getNames();
    throw std::runtime_error("unknown extractor: " + tool);
}

int run_extractor(const std::string& tool, const std::string& filename, const std::string& ext,
                  ResourceLimits& limits, Mode mode) {
    IExtractor* extractor = nullptr;
    try {
        extractor = make_extractor(tool, ext, filename);
        extractor->run();
    } catch (TimeLimitExceeded&) {
        delete extractor;
        if (mode == Mode::GBD) {
            std::cout << "status timeout" << std::endl;
            std::cout << "runtime " << limits.get_runtime() << std::endl;
            return 0;
        }
        std::cerr << "Time Limit Exceeded" << std::endl;
        return 1;
    } catch (MemoryLimitExceeded&) {
        delete extractor;
        if (mode == Mode::GBD) {
            std::cout << "status memout" << std::endl;
            std::cout << "runtime " << limits.get_runtime() << std::endl;
            return 0;
        }
        std::cerr << "Memory Limit Exceeded" << std::endl;
        return 1;
    }

    const std::vector<std::string> names = extractor->getNames();
    const std::vector<double> features = extractor->getFeatures();

    if (mode == Mode::GBD) {
        for (size_t i = 0; i < names.size(); ++i) {
            std::cout << names[i] << " " << format_value(features[i]) << std::endl;
        }
        std::cout << "status success" << std::endl;
        std::cout << "runtime " << limits.get_runtime() << std::endl;
    } else {
        for (const std::string& name : names) std::cout << name << " ";
        std::cout << std::endl;
        for (double feature : features) std::cout << feature << " ";
        std::cout << std::endl;
    }
    delete extractor;
    return 0;
}

std::vector<std::string> checksani_feature_names() {
    return {"header_consistent", "whitespace_normalised", "no_comment",
            "no_tautological_clause", "no_duplicate_literals", "no_empty_clause"};
}

/* checksani reports normalisation/sanitation flags rather than numeric features. */
int run_checksani(const std::string& filename, Mode mode) {
    CNF::SaniCheck ana(filename.c_str(), true);
    ana.run();
    struct Flag { const char* name; bool value; };
    const std::vector<Flag> flags = {
        {"header_consistent", ana.getFeature("head_vars") == ana.getFeature("norm_vars") &&
                              ana.getFeature("head_clauses") == ana.getFeature("norm_clauses")},
        {"whitespace_normalised", ana.getFeature("whitespace_normalised") == 1.0},
        {"no_comment", ana.getFeature("has_comment") == 0.0},
        {"no_tautological_clause", ana.getFeature("has_tautological_clause") == 0.0},
        {"no_duplicate_literals", ana.getFeature("has_duplicate_literals") == 0.0},
        {"no_empty_clause", ana.getFeature("has_empty_clause") == 0.0},
    };
    if (mode == Mode::GBD) {
        for (const Flag& f : flags) std::cout << f.name << " " << (f.value ? "yes" : "no") << std::endl;
        std::cout << "status success" << std::endl;
    } else {
        std::cout << "hash " << CNF::gbdhash(filename.c_str()) << std::endl;
        std::cout << "filename " << filename << std::endl;
        for (const Flag& f : flags) std::cout << f.name << " " << (f.value ? "yes" : "no") << std::endl;
    }
    return 0;
}


/* --- Identifiers --------------------------------------------------------------------------- */

int run_identify(const std::string& filename, const std::string& ext) {
    std::string hash;
    if (ext == ".cnf" || ext == ".wecnf") hash = CNF::gbdhash(filename.c_str());
    else if (ext == ".opb") hash = OPB::gbdhash(filename.c_str());
    else if (ext == ".qcnf" || ext == ".qdimacs") hash = PQBF::gbdhash(filename.c_str());
    else if (ext == ".wcnf") hash = WCNF::gbdhash(filename.c_str());
    else throw std::runtime_error("identify: unsupported format " + ext);
    std::cout << hash << std::endl;
    return 0;
}

int run_isohash(const std::string& filename, const std::string& ext, Mode mode) {
    std::string value;
    if (ext == ".wcnf") value = WCNF::isohash(filename.c_str());
    else if (ext == ".cnf") value = CNF::isohash(filename.c_str());
    else throw std::runtime_error("isohash: unsupported format " + ext);
    if (mode == Mode::GBD) std::cout << "isohash " << value << std::endl;
    else std::cout << value << std::endl;
    return 0;
}

int run_isohash2(const std::string& filename, const std::string& ext, argparse::ArgumentParser& args, Mode mode) {
    if (ext != ".cnf") throw std::runtime_error("isohash2: unsupported format " + ext);
    CNF::IsoHash2Settings config;
    if (auto max_iters = args.present<int>("--max-iters")) config.max_iterations = *max_iters;
    const std::string value = CNF::isohash2(filename.c_str(), config);
    if (mode == Mode::GBD) std::cout << "isohash2 " << value << std::endl;
    else std::cout << value << std::endl;
    return 0;
}


/* --- Transformers -------------------------------------------------------------------------- */

/* Parse the -z/--compress value into a libarchive compression format. */
CompressionFormat compression_format(const std::string& name) {
    if (name == "xz") return CompressionFormat::XZ;
    if (name == "gz") return CompressionFormat::GZIP;
    if (name == "bz2") return CompressionFormat::BZIP2;
    throw std::runtime_error("unknown compression format: " + name + " (expected none, xz, gz, or bz2)");
}

/* Run a transformer. The transformer classes emit the produced instance to std::cout; the driver
 * points std::cout's buffer at the chosen destination so the instance streams there directly,
 * without buffering the whole payload:
 *   - human/CLI mode without -o: the instance is the primary output and streams to stdout;
 *   - -o (plain): the instance streams to the output file;
 *   - -o with -z <xz|gz|bz2>: the instance streams through the matching libarchive compressor.
 * In --gbd mode stdout instead carries the feature/metadata stream, so -o is required (and gbd
 * always passes it). */
int run_transformer(const std::string& tool, const std::string& filename, const std::string& output,
                    const std::string& compress, ResourceLimits& limits, Mode mode) {
    const bool has_output = !(output.empty() || output == "-");
    if (mode == Mode::GBD && !has_output) {
        throw std::runtime_error("transformer requires -o/--output in --gbd mode");
    }

    // Set up the destination stream for the produced instance.
    std::string local;
    std::ofstream file_stream;
    std::unique_ptr<StreamCompressor> compressor;
    std::unique_ptr<CompressorStreamBuf> compressor_buf;
    std::streambuf* const real_cout = std::cout.rdbuf();
    std::streambuf* sink = real_cout;  // CLI without -o: stream the instance to stdout

    if (has_output) {
        if (compress == "none") {
            local = output;
            file_stream.open(output, std::ofstream::out);
            if (!file_stream) throw std::runtime_error("Could not open output file: " + output);
            sink = file_stream.rdbuf();
        } else {
            const CompressionFormat format = compression_format(compress);
            const std::string suffix = compression_suffix(format);
            local = output;
            if (local.size() < suffix.size() || local.substr(local.size() - suffix.size()) != suffix) {
                local += suffix;
            }
            compressor = std::make_unique<StreamCompressor>(local.c_str(), 0, format);
            compressor_buf = std::make_unique<CompressorStreamBuf>(*compressor);
            sink = compressor_buf.get();
        }
    }

    std::vector<std::pair<std::string, std::string>> derived;
    std::cout.rdbuf(sink);
    try {
        if (tool == "cnf2kis") {
            IndependentSetFromCNF gen(filename.c_str());
            derived.emplace_back("nodes", format_value(gen.numNodes()));
            derived.emplace_back("edges", format_value(gen.numEdges()));
            derived.emplace_back("k", format_value(gen.minK()));
            gen.generate_independent_set_problem(nullptr);
        } else if (tool == "sanitize") {
            CNF::Sanitiser(filename.c_str(), nullptr).run();
        } else if (tool == "normalize") {
            CNF::Normaliser(filename.c_str(), nullptr).run();
        } else if (tool == "cnf2bip") {
            CNF::cnf2bip gen(filename.c_str(), "");
            derived.emplace_back("nodes", format_value(gen.getFeature("nodes")));
            derived.emplace_back("edges", format_value(gen.getFeature("edges")));
            gen.run();
        } else {
            std::cout.rdbuf(real_cout);
            throw std::runtime_error("unknown transformer: " + tool);
        }
    } catch (...) {
        std::cout.rdbuf(real_cout);
        throw;
    }
    std::cout.flush();
    std::cout.rdbuf(real_cout);

    // Finalise the destination (order matters: flush the compressor buffer before closing it).
    if (compressor_buf) compressor_buf->pubsync();
    if (compressor) compressor->close();
    if (file_stream.is_open()) file_stream.close();

    if (!has_output) return 0;  // CLI: the instance was streamed to stdout

    const std::string hash = CNF::gbdhash(local.c_str());
    if (mode == Mode::GBD) {
        std::cout << "local " << local << std::endl;
        std::cout << "hash " << hash << std::endl;
        for (const auto& [name, value] : derived) std::cout << name << " " << value << std::endl;
        if (tool == "cnf2kis" || tool == "sanitize") {
            std::cout << "to_cnf " << CNF::gbdhash(filename.c_str()) << std::endl;
        }
        std::cout << "status success" << std::endl;
        std::cout << "runtime " << limits.get_runtime() << std::endl;
    } else {
        std::cerr << "Produced " << local << " with hash " << hash << std::endl;
    }
    return 0;
}

std::vector<std::pair<std::string, std::string>> transformer_feature_names(const std::string& tool) {
    /* the second entry is the default; an empty default marks a non-unique (1:n) feature. */
    if (tool == "cnf2kis") {
        return {{"local", ""}, {"to_cnf", ""}, {"nodes", "empty"}, {"edges", "empty"}, {"k", "empty"}};
    }
    if (tool == "sanitize") return {{"local", ""}, {"to_cnf", ""}};
    if (tool == "normalize") return {{"local", ""}};
    if (tool == "cnf2bip") return {{"local", ""}, {"nodes", "empty"}, {"edges", "empty"}};
    throw std::runtime_error("unknown transformer: " + tool);
}


/* --- Dispatch helpers ---------------------------------------------------------------------- */

bool is_extractor(const std::string& tool) {
    return tool == "base" || tool == "gate" || tool == "wcnfbase" || tool == "opbbase";
}

bool is_transformer(const std::string& tool) {
    return tool == "cnf2kis" || tool == "sanitize" || tool == "normalize" || tool == "cnf2bip";
}

int print_feature_names(const std::string& tool, Mode mode) {
    if (is_extractor(tool)) {
        for (const std::string& name : extractor_feature_names(tool)) {
            std::cout << name << (mode == Mode::GBD ? " empty" : "") << std::endl;
        }
        return 0;
    }
    if (tool == "checksani") {
        for (const std::string& name : checksani_feature_names()) {
            std::cout << name << (mode == Mode::GBD ? " empty" : "") << std::endl;
        }
        return 0;
    }
    if (tool == "isohash") { std::cout << "isohash" << (mode == Mode::GBD ? " empty" : "") << std::endl; return 0; }
    if (tool == "isohash2") { std::cout << "isohash2" << (mode == Mode::GBD ? " empty" : "") << std::endl; return 0; }
    if (is_transformer(tool)) {
        for (const auto& [name, def] : transformer_feature_names(tool)) {
            if (mode == Mode::GBD && !def.empty()) std::cout << name << " " << def << std::endl;
            else std::cout << name << std::endl;
        }
        return 0;
    }
    throw std::runtime_error("--feature-names not supported for tool: " + tool);
}

}  // namespace


int main(int argc, char** argv) {
    const std::string invocation_tool = tool_from_invocation(argc > 0 ? argv[0] : "");

    argparse::ArgumentParser program("gbdc");

    /* In legacy mode (invoked as "gbdc"/"gbdctool") the tool is the first positional argument. */
    if (invocation_tool.empty()) {
        program.add_argument("tool").help(
            "Tool: identify, isohash, isohash2, normalize, sanitize, checksani, "
            "cnf2kis, cnf2bip, base, gate, wcnfbase, opbbase");
    }
    program.add_argument("file").remaining().help("Path to input file");
    program.add_argument("-o", "--output").default_value(std::string("-"))
        .help("Output file for transformers (default: stderr)");
    program.add_argument("-z", "--compress").default_value(std::string("none"))
        .help("Compression for -o output: none, xz, gz, or bz2");
    program.add_argument("-t", "--tlim").default_value(0).scan<'i', int>().help("Time limit in seconds");
    program.add_argument("-m", "--mlim").default_value(0).scan<'i', int>().help("Memory limit in MB");
    program.add_argument("-f", "--flim").default_value(0).scan<'i', int>().help("Output file size limit in MB");
    program.add_argument("--max-iters").scan<'i', int>().help("Maximum isohash2 iterations");
    program.add_argument("--gbd").default_value(false).implicit_value(true)
        .help("Emit machine-readable output for gbd");
    program.add_argument("--feature-names").default_value(false).implicit_value(true)
        .help("Print the features this tool produces and exit");

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    std::string tool = invocation_tool.empty() ? program.get("tool") : invocation_tool;
    tool = canonical_tool(tool);
    const Mode mode = program.get<bool>("--gbd") ? Mode::GBD : Mode::HUMAN;

    if (program.get<bool>("--feature-names")) {
        try {
            return print_feature_names(tool, mode);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }

    const auto files = program.present<std::vector<std::string>>("file");
    if (!files || files->empty()) {
        std::cerr << "No input file given" << std::endl;
        std::cerr << program;
        return 1;
    }
    const std::string filename = files->front();
    const std::string output = program.get("output");
    const std::string compress = program.get("compress");

    ResourceLimits limits(program.get<int>("tlim"), program.get<int>("mlim"), program.get<int>("flim"));
    limits.set_rlimits();

    const std::string ext = detect_extension(filename);
    std::cerr << "c Running: " << tool << " " << filename << std::endl;

    try {
        if (is_extractor(tool)) return run_extractor(tool, filename, ext, limits, mode);
        if (tool == "checksani") return run_checksani(filename, mode);
        if (tool == "identify") return run_identify(filename, ext);
        if (tool == "isohash") return run_isohash(filename, ext, mode);
        if (tool == "isohash2") return run_isohash2(filename, ext, program, mode);
        if (is_transformer(tool)) return run_transformer(tool, filename, output, compress, limits, mode);
        std::cerr << "Unknown tool: " << tool << std::endl;
        return 1;
    } catch (std::bad_alloc&) {
        std::cerr << "Memory Limit Exceeded" << std::endl;
        return 1;
    } catch (MemoryLimitExceeded&) {
        std::cerr << "Memory Limit Exceeded" << std::endl;
        return 1;
    } catch (TimeLimitExceeded&) {
        std::cerr << "Time Limit Exceeded" << std::endl;
        return 1;
    } catch (FileSizeLimitExceeded&) {
        if (!output.empty() && output != "-") std::remove(output.c_str());
        std::cerr << "File Size Limit Exceeded" << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
