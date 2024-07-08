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

#include "src/util/threadpool/ThreadPool.h"

#include "src/extract/CNFBaseFeatures.h"
#include "src/extract/CNFGateFeatures.h"
#include "src/extract/WCNFBaseFeatures.h"
#include "src/extract/OPBBaseFeatures.h"

#include "src/transform/IndependentSet.h"
#include "src/transform/Normalize.h"
#include "src/util/ResourceLimits.h"

#include "src/util/pybind11/include/pybind11/pybind11.h"
#include "src/util/pybind11/include/pybind11/stl.h"
#include "src/util/pybind11/include/pybind11/functional.h"
#include "gbdlib.h"

namespace py = pybind11;
namespace tp = threadpool;

template <typename Ret>
static void bind_MPSCQueue(py::module &m, const std::string &type_name)
{
    using q = MPSCQueue<tp::result_t<Ret>>;
    py::class_<q, std::shared_ptr<q>>(m, type_name.c_str())
        .def(py::init<std::uint16_t>())
        .def("pop", &q::pop, "Pop an element from the queue in a synchronized manner.")
        .def("empty", &q::empty, "Returns true if queue is empty, false otherwise.")
        .def("done", &q::done, "Returns true if producers are done, false otherwise.");
}

py::dict cnf2kis(const std::string filename, const std::string output, const size_t maxEdges, const size_t maxNodes)
{
    py::dict dict;

    IndependentSetFromCNF gen(filename.c_str());
    unsigned nNodes = gen.numNodes();
    unsigned nEdges = gen.numEdges();
    unsigned minK = gen.minK();

    dict[py::str("nodes")] = nNodes;
    dict[py::str("edges")] = nEdges;
    dict[py::str("k")] = minK;

    if ((maxEdges > 0 && nEdges > maxEdges) || (maxNodes > 0 && nNodes > maxNodes))
    {
        dict[py::str("hash")] = "fileout";
        return dict;
    }

    gen.generate_independent_set_problem(output.c_str());
    dict[py::str("local")] = output;

    dict[py::str("hash")] = CNF::gbdhash(output.c_str());
    return dict;
}

template <typename Extractor, int type_id = 0>
py::dict extract_features(const std::string filepath, const size_t rlim, const size_t mlim)
{
    constexpr auto d = get_runtime_desc<type_id>();
    py::dict emergency;
    ResourceLimits limits(rlim, mlim);
    limits.set_rlimits();
    try
    {
        Extractor stats(filepath.c_str());
        stats.extract();
        py::dict dict;
        const auto names = stats.getNames();
        const auto features = stats.getFeatures();
        dict[py::str(d)] = (double)limits.get_runtime();
        for (size_t i = 0; i < features.size(); ++i)
        {
            dict[py::str(names[i])] = features[i];
        }
        return dict;
    }
    catch (TimeLimitExceeded &e)
    {
        emergency[py::str(d)] = "timeout";
        return emergency;
    }
    catch (MemoryLimitExceeded &e)
    {
        emergency[py::str(d)] = "memout";
        return emergency;
    }
}

PYBIND11_MODULE(gbdc, m)
{
    bind_MPSCQueue<tp::extract_ret_t>(m, "MPSCQueue_Extract");
    m.def("tp_extract_base_features", &tp_extract<CNF::BaseFeatures>, "Extract base features", py::arg("mlim"), py::arg("jobs"), py::arg("args"));
    m.def("extract_base_features", &extract_features<CNF::BaseFeatures<>>, "Extract base features", py::arg("filepath"), py::arg("rlim"), py::arg("mlim"));
    m.def("tp_extract_gate_features", &tp_extract<CNFGateFeatures>, "Extract gate features", py::arg("mlim"), py::arg("jobs"), py::arg("args"));
    m.def("extract_gate_features", &extract_features<CNFGateFeatures<>, 1>, "Extract gate features", py::arg("filepath"), py::arg("rlim"), py::arg("mlim"));
    m.def("tp_extract_wcnf_base_features", &tp_extract<WCNF::BaseFeatures>, "Extract wcnf base features", py::arg("mlim"), py::arg("jobs"), py::arg("args"));
    m.def("extract_wcnf_base_features", &extract_features<WCNF::BaseFeatures<>>, "Extract wcnf base features", py::arg("filepath"), py::arg("rlim"), py::arg("mlim"));
    m.def("tp_extract_opb_base_features", &tp_extract<OPB::BaseFeatures>, "Extract opb base features", py::arg("mlim"), py::arg("jobs"), py::arg("args"));
    m.def("extract_opb_base_features", &extract_features<OPB::BaseFeatures<>>, "Extract opb base features", py::arg("filepath"), py::arg("rlim"), py::arg("mlim"));
    m.def("version", &version, "Return current version of gbdc.");
    m.def("cnf2kis", &cnf2kis, "Create k-ISP Instance from given CNF Instance.", py::arg("filename"), py::arg("output"), py::arg("maxEdges"), py::arg("maxNodes"));
    m.def("sanitize", &sanitize, "Print sanitized, i.e., no duplicate literals in clauses and no tautologic clauses, CNF to stdout.", py::arg("filename"));
    m.def("base_feature_names", &feature_names<CNF::BaseFeatures<>>, "Get Base Feature Names");
    m.def("gate_feature_names", &feature_names<CNFGateFeatures<>, 1>, "Get Gate Feature Names");
    m.def("wcnf_base_feature_names", &feature_names<WCNF::BaseFeatures<>>, "Get WCNF Base Feature Names");
    m.def("opb_base_feature_names", &feature_names<OPB::BaseFeatures<>>, "Get OPB Base Feature Names");
    m.def("gbdhash", &CNF::gbdhash, "Calculates GBD-Hash (md5 of normalized file) of given DIMACS CNF file.", py::arg("filename"));
    m.def("isohash", &WCNF::isohash, "Calculates ISO-Hash (md5 of sorted degree sequence) of given DIMACS CNF file.", py::arg("filename"));
    m.def("opbhash", &OPB::gbdhash, "Calculates OPB-Hash (md5 of normalized file) of given OPB file.", py::arg("filename"));
    m.def("pqbfhash", &PQBF::gbdhash, "Calculates PQBF-Hash (md5 of normalized file) of given PQBF file.", py::arg("filename"));
    m.def("wcnfhash", &WCNF::gbdhash, "Calculates WCNF-Hash (md5 of normalized file) of given WCNF file.", py::arg("filename"));
    m.def("wcnfisohash", &WCNF::isohash, "Calculates WCNF ISO-Hash of given WCNF file.", py::arg("filename"));
}
