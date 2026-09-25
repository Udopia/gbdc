#include "src/util/OutputWrapper.h"

#include "cnf2lig.h"

CNF::cnf2lig::cnf2lig(const char* filename, const char* output) : F(), filename_(filename), output_(output) { 
    F.readDimacsFromFile(filename);
    setFeature("nodes", 2 * F.nVars() + F.nClauses());
    setFeature("edges", F.nVars() + F.nLits());
}

CNF::cnf2lig::~cnf2lig() { }

void CNF::cnf2lig::run() {
    std::string outputStr(output_);
    OutputWrapper out(&outputStr);

    size_t n = F.nVars();

    out << "p edge " << 2 * n + F.nClauses() << " " << n + F.nLits() << std::endl;

    if (n > 0) {
        out << "n 1 0" << std::endl;
    }

    for (size_t i = 0; i < F.nClauses(); i++) {
        out << "n " << 2 * n + i + 1 << " 1" << std::endl;
    }

    for (size_t v = 1; v <= n; v++) {
        out << "e " << 2 * v - 1 << " " << 2 * v << std::endl;
    }

    size_t clause_id = 2 * n + 1;
    for (Cl* clause : F) {
        for (size_t i = 0; i < clause->size(); i++) {
            if ((*clause)[i].sign()) {
                out << "e " << 2 * (*clause)[i].var() << " " << clause_id << std::endl;
            } else {
                out << "e " << 2 * (*clause)[i].var() - 1 << " " << clause_id << std::endl;
            }
        }
        clause_id++;
    }
}