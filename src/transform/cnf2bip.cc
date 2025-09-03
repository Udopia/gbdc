/**
 * MIT License
 * Copyright (c) 2025 Ashlin Iser 
 */

#include "src/util/OutputWrapper.h"

#include "cnf2bip.h"

CNF::cnf2bip::cnf2bip(const char* filename, const char* output) : F(), filename_(filename), output_(output) { 
    F.readDimacsFromFile(filename);
    setFeature("nodes", F.nVars() + F.nClauses());
    setFeature("edges", F.nLits());
}

CNF::cnf2bip::~cnf2bip() { }

void CNF::cnf2bip::run() {
    std::string outputStr(output_);
    OutputWrapper out(&outputStr);

    out << "c directed bipartite graph representation from cnf" << std::endl;
    out << "p edge " << F.nVars() + F.nClauses() << " " << F.nLits() << std::endl;

    size_t clause_id = F.nVars() + 1;
    for (Cl* clause : F) {
        for (size_t i = 0; i < clause->size(); i++) {
            if ((*clause)[i].sign()) {
                out << "e " << (*clause)[i].var() << " " << clause_id << std::endl;
            } else {
                out << "e " << clause_id << " " << (*clause)[i].var() << std::endl;
            }
        }
        clause_id++;
    }
}