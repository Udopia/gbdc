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

#include "Python.h"

#include "src/identify/GBDHash.h"
#include "src/identify/ISOHash.h"

#include "src/util/CNFFormula.h"
#include "src/util/ResourceLimits.h"
#include "src/util/py_util.h"

#include "src/extract/CNFBaseFeatures.h"
#include "src/extract/CNFGateFeatures.h"

#include "src/transform/IndependentSet.h"
#include "src/transform/Normalize.h"

static PyObject* version(PyObject* self) {
    return pytype(1);
}

static PyObject* gbdhash(PyObject* self, PyObject* arg) {
    const char* filename;
    PyArg_ParseTuple(arg, "s", &filename);
    std::string result = CNF::gbdhash(filename);
    return pytype(result.c_str());
}

static PyObject* isohash(PyObject* self, PyObject* arg) {
    const char* filename;
    PyArg_ParseTuple(arg, "s", &filename);
    std::string result = CNF::isohash(filename);
    return pytype(result.c_str());
}

static PyObject* opbhash(PyObject* self, PyObject* arg) {
    const char* filename;
    PyArg_ParseTuple(arg, "s", &filename);
    std::string result = OPB::gbdhash(filename);
    return pytype(result.c_str());
}

static PyObject* pqbfhash(PyObject* self, PyObject* arg) {
    const char* filename;
    PyArg_ParseTuple(arg, "s", &filename);
    std::string result = PQBF::gbdhash(filename);
    return pytype(result.c_str());
}

static PyObject* wcnfhash(PyObject* self, PyObject* arg) {
    const char* filename;
    PyArg_ParseTuple(arg, "s", &filename);
    std::string result = WCNF::gbdhash(filename);
    return pytype(result.c_str());
}


static PyObject* extract_base_features(PyObject* self, PyObject* arg) {
    const char* filename;
    unsigned rlim = 0, mlim = 0;
    PyArg_ParseTuple(arg, "s|II", &filename, &rlim, &mlim);

    PyObject *emergency = pydict();
    pydict(emergency, "base_features_runtime", "memout");

    ResourceLimits limits(rlim, mlim);
    limits.set_rlimits();
    try {
        CNFBaseFeatures stats(filename);
        stats.extract();
        std::vector<double> record = stats.getFeatures();
        std::vector<std::string> names = stats.getNames();
        PyObject *dict = pydict();
        pydict(dict, "base_features_runtime", limits.get_runtime());
        for (unsigned int i = 0; i < record.size(); i++) {
            pydict(dict, names[i].c_str(), record[i]);
        }
        return dict;
    } catch (TimeLimitExceeded& e) {
        pydict(emergency, "base_features_runtime", "timeout");
        return emergency;
    } catch (MemoryLimitExceeded& e) {
        return emergency;
    }
}


static PyObject* extract_gate_features(PyObject* self, PyObject* arg) {
    const char* filename;
    unsigned rlim = 0, mlim = 0;
    PyArg_ParseTuple(arg, "s|II", &filename, &rlim, &mlim);

    PyObject *emergency = pydict();
    pydict(emergency, "gate_features_runtime", "memout");

    ResourceLimits limits(rlim, mlim);
    limits.set_rlimits();
    try {
        CNFGateFeatures stats(filename);
        stats.extract();
        std::vector<double> record = stats.getFeatures();
        std::vector<std::string> names = stats.getNames();
        PyObject *dict = pydict();
        for (unsigned int i = 0; i < record.size(); i++) {
            pydict(dict, names[i].c_str(), record[i]);
        }
        pydict(dict, "gate_features_runtime", limits.get_runtime());
        return dict;
    } catch (TimeLimitExceeded& e) {
        pydict(emergency, "gate_features_runtime", "timeout");
        return emergency;
    } catch (MemoryLimitExceeded& e) {
        return emergency;
    }
}

static PyObject* base_feature_names(PyObject* self) {
    PyObject *list = PyList_New(0);
    PyList_Append(list, pytype("base_features_runtime"));
    CNFBaseFeatures stats("");
    std::vector<std::string> names = stats.getNames();
    for (unsigned int i = 0; i < names.size(); i++) {
        PyList_Append(list, pytype(names[i].c_str()));
    }
    return list;
}

static PyObject* gate_feature_names(PyObject* self) {
    PyObject *list = PyList_New(0);
    PyList_Append(list, pytype("gate_features_runtime"));
    CNFGateFeatures stats("");
    std::vector<std::string> names = stats.getNames();
    for (unsigned int i = 0; i < names.size(); i++) {
        PyList_Append(list, pytype(names[i].c_str()));
    }
    return list;
}


static PyObject* cnf2kis(PyObject* self, PyObject* arg) {
    const char* filename;
    const char* output;
    unsigned maxEdges, maxNodes;
    unsigned rlim = 0, mlim = 0, flim = 0;
    PyArg_ParseTuple(arg, "ssII|III", &filename, &output, &maxEdges, &maxNodes, &rlim, &mlim, &flim);

    PyObject *dict = pydict();
    pydict(dict, "nodes", 0);
    pydict(dict, "edges", 0);
    pydict(dict, "k", 0);

    ResourceLimits limits(rlim, mlim, flim);
    limits.set_rlimits();
    try {
        IndependentSetFromCNF gen(filename);
        unsigned nNodes = gen.numNodes();
        unsigned nEdges = gen.numEdges();
        unsigned minK = gen.minK();

        pydict(dict, "nodes", nNodes);
        pydict(dict, "edges", nEdges);
        pydict(dict, "k", minK);

        if ((maxEdges > 0 && nEdges > maxEdges) || (maxNodes > 0 && nNodes > maxNodes)) {
            pydict(dict, "hash", "fileout");
            return dict;
        }

        gen.generate_independent_set_problem(output);
        pydict(dict, "local", output);

        std::string hash = CNF::gbdhash(output);
        pydict(dict, "hash", hash.c_str());

        return dict;
    } catch (TimeLimitExceeded& e) {
        std::remove(output);
        pydict(dict, "hash", "timeout");
        return dict;
    } catch (MemoryLimitExceeded& e) {
        std::remove(output);
        pydict(dict, "hash", "memout");
        return dict;
    } catch (FileSizeLimitExceeded& e) {
        std::remove(output);
        pydict(dict, "hash", "fileout");
        return dict;
    }
}

static PyObject* print_sanitized(PyObject* self, PyObject* arg) {
    const char* filename;
    unsigned rlim = 0, mlim = 0;
    PyArg_ParseTuple(arg, "s|II", &filename, &rlim, &mlim);

    ResourceLimits limits(rlim, mlim);
    limits.set_rlimits();
    try {
        sanitize(filename);
        Py_RETURN_TRUE;
    } catch (TimeLimitExceeded& e) {
        Py_RETURN_FALSE;
    } catch (MemoryLimitExceeded& e) {
        Py_RETURN_FALSE;
    }
}

static PyMethodDef myMethods[] = {
    {"extract_gate_features", extract_gate_features, METH_VARARGS, "Extract Gate Features."},
    {"extract_base_features", extract_base_features, METH_VARARGS, "Extract Base Features."},
    {"base_feature_names", (PyCFunction)base_feature_names, METH_NOARGS, "Get Base Feature Names."},
    {"gate_feature_names", (PyCFunction)gate_feature_names, METH_NOARGS, "Get Gate Feature Names."},
    {"sanitize", print_sanitized, METH_VARARGS, "Print sanitized, i.e., no duplicate literals in clauses and no tautologic clauses, CNF to stdout."},
    {"cnf2kis", cnf2kis, METH_VARARGS, "Create k-ISP Instance from given CNF Instance."},
    {"gbdhash", gbdhash, METH_VARARGS, "Calculates GBD-Hash (md5 of normalized file) of given DIMACS CNF file."},
    {"isohash", isohash, METH_VARARGS, "Calculates ISO-Hash (md5 of sorted degree sequence) of given DIMACS CNF file."},
    {"opbhash", opbhash, METH_VARARGS, "Calculates OPB-Hash (md5 of normalized file) of given OPB file."},
    {"pqbfhash", pqbfhash, METH_VARARGS, "Calculates PQBF-Hash (md5 of normalized file) of given PQBF file."},
    {"wcnfhash", wcnfhash, METH_VARARGS, "Calculates WCNF-Hash (md5 of normalized file) of given WCNF file."},
    {"version", (PyCFunction)version, METH_NOARGS, "Returns Version"},
    {nullptr, nullptr, 0, nullptr}
};

static struct PyModuleDef myModule = {
    PyModuleDef_HEAD_INIT,
    "gbdc",
    "GBDC Accelerator Module",
    -1,
    myMethods
};

PyMODINIT_FUNC PyInit_gbdc(void) {
    return PyModule_Create(&myModule);
}
