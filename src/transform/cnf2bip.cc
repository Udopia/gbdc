/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

 #include <iostream>
 #include <fstream>
 #include <memory>

#include "cnf2bip.h"

CNF::cnf2bip::cnf2bip(const char* filename, const char* output) : F(), filename_(filename), output_(output), features(), names() { 
    F.readDimacsFromFile(filename);
    names.push_back("nodes");
    names.push_back("edges");
    features.push_back(F.nVars() + F.nClauses());
    features.push_back(0);
}

CNF::cnf2bip::~cnf2bip() {

}

std::vector<double> CNF::cnf2bip::getFeatures() const {
    return features;
}

std::vector<std::string> CNF::cnf2bip::getNames() const {
    return names;
}

void CNF::cnf2bip::run() {
    std::shared_ptr<std::ostream> of;
    if (output_ != nullptr) {
        of.reset(new std::ofstream(output_, std::ofstream::out), [](std::ostream* p){ delete p; });
    } else {
        of.reset(&std::cout, [](...){});
    }

    *of << "c directed bipartite graph representation from cnf" << std::endl;
    *of << "p edge " << F.nVars() + F.nClauses() << std::endl;

    unsigned clause_id = F.nVars() + 1;
    for (Cl* clause : F) {
        for (unsigned i = 0; i < clause->size(); i++) {
            if ((*clause)[i].sign()) {
                *of << "e " << (*clause)[i].var() << " " << clause_id << std::endl;
            } else {
                *of << "e " << clause_id << " " << (*clause)[i].var() << std::endl;
            }
            features[1]++;
        }
        clause_id++;
    }
}