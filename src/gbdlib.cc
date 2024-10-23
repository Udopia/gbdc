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

#include <cstdio>
#include <cstdint>
#include <string>
#include <chrono>
#include <regex>
#include <fstream>
#include <memory>
#include <sstream>
#include <future>
#include <unordered_map>
#include <variant>

#include "src/identify/GBDHash.h"
#include "src/identify/ISOHash.h"

#include "src/extract/CNFBaseFeatures.h"
#include "src/extract/CNFGateFeatures.h"
#include "src/extract/WCNFBaseFeatures.h"
#include "src/extract/OPBBaseFeatures.h"

#include "src/transform/IndependentSet.h"
#include "src/transform/Normalize.h"
#include "src/util/ResourceLimits.h"

// #include "src/util/pybind11/include/pybind11/pybind11.h"
// #include "src/util/pybind11/include/pybind11/stl.h"
// #include "src/util/pybind11/include/pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"

namespace py = pybind11;

std::string version() {
    std::ifstream file("setup.py");
    if (!file.is_open()) {
        return "Error: Could not open setup.py";
    }

    std::string line;
    std::regex version_regex(R"(version\s*=\s*['"]([^'"]+)['"])");
    std::smatch match;

    while (std::getline(file, line)) {
        if (std::regex_search(line, match, version_regex)) {
            return match[1].str();
        }
    }

    return "Error: Version not found in setup.py";
}

template <typename Extractor>
auto feature_names() {
    auto ex = Extractor("");
    auto names = ex.getNames();
    names.push_back(ex.getRuntimeDesc());
    return names;
}

py::dict cnf2kis(const std::string filename, const std::string output) {
    py::dict dict;

    IndependentSetFromCNF gen(filename.c_str());

    dict[py::str("nodes")] = gen.numNodes();
    dict[py::str("edges")] = gen.numEdges();
    dict[py::str("k")] = gen.minK();

    gen.generate_independent_set_problem(output.c_str());
    dict[py::str("local")] = output;
    dict[py::str("hash")] = CNF::gbdhash(output.c_str());
    return dict;
}

template <typename Extractor>
py::dict extract_features(const std::string filepath, const size_t rlim, const size_t mlim) {
    py::dict dict;
    Extractor stats(filepath.c_str());
    ResourceLimits limits(rlim, mlim);
    limits.set_rlimits();
    try {
        stats.extract();
        dict[py::str(stats.getRuntimeDesc())] = (double)limits.get_runtime();
        const auto names = stats.getNames();
        const auto features = stats.getFeatures();
        for (size_t i = 0; i < features.size(); ++i) {
            dict[py::str(names[i])] = features[i];
        }
    }
    catch (TimeLimitExceeded &e) {
        dict[py::str(stats.getRuntimeDesc())] = "timeout";
    }
    catch (MemoryLimitExceeded &e) {
        dict[py::str(stats.getRuntimeDesc())] = "memout";
    }
    return dict;
}

PYBIND11_MODULE(gbdc, m) {
    m.doc() = "GBDC Python Bindings";
    m.def("extract_base_features", &extract_features<CNF::BaseFeatures>, "Extract cnf base features", py::arg("filepath"), py::arg("rlim"), py::arg("mlim"));
    m.def("extract_gate_features", &extract_features<CNF::GateFeatures>, "Extract cnf gate features", py::arg("filepath"), py::arg("rlim"), py::arg("mlim"));
    m.def("extract_wcnf_base_features", &extract_features<WCNF::BaseFeatures>, "Extract wcnf base features", py::arg("filepath"), py::arg("rlim"), py::arg("mlim"));
    m.def("extract_opb_base_features", &extract_features<OPB::BaseFeatures>, "Extract opb base features", py::arg("filepath"), py::arg("rlim"), py::arg("mlim"));
    m.def("version", &version, "Return current version of gbdc.");
    m.def("cnf2kis", &cnf2kis, "Create k-ISP Instance from given CNF Instance.", py::arg("filename"), py::arg("output"));
    m.def("sanitize", &sanitize, "Print sanitized, i.e., no duplicate literals in clauses and no tautologic clauses, CNF to stdout.", py::arg("filename"));
    m.def("base_feature_names", &feature_names<CNF::BaseFeatures>, "Get Base Feature Names");
    m.def("gate_feature_names", &feature_names<CNF::GateFeatures>, "Get Gate Feature Names");
    m.def("wcnf_base_feature_names", &feature_names<WCNF::BaseFeatures>, "Get WCNF Base Feature Names");
    m.def("opb_base_feature_names", &feature_names<OPB::BaseFeatures>, "Get OPB Base Feature Names");
    m.def("gbdhash", &CNF::gbdhash, "Calculates GBD-Hash (md5 of normalized file) of given DIMACS CNF file.", py::arg("filename"));
    m.def("isohash", &CNF::isohash, "Calculates ISO-Hash (md5 of sorted degree sequence) of given DIMACS CNF file.", py::arg("filename"));
    m.def("opbhash", &OPB::gbdhash, "Calculates OPB-Hash (md5 of normalized file) of given OPB file.", py::arg("filename"));
    m.def("pqbfhash", &PQBF::gbdhash, "Calculates PQBF-Hash (md5 of normalized file) of given PQBF file.", py::arg("filename"));
    m.def("wcnfhash", &WCNF::gbdhash, "Calculates WCNF-Hash (md5 of normalized file) of given WCNF file.", py::arg("filename"));
    m.def("wcnfisohash", &WCNF::isohash, "Calculates WCNF ISO-Hash of given WCNF file.", py::arg("filename"));
}
