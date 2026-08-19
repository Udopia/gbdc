/**
 * Standalone gate feature extractor — part of extras/gate.
 * Build with: cmake -B build && cmake --build build
 * Usage: gbdc-gate <file.cnf>
 */

#include <iostream>
#include <string>

#include "CNFGateFeatures.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <file.cnf>" << std::endl;
        return 1;
    }
    const std::string filename = argv[1];
    CNF::GateFeatures stats(filename.c_str());
    stats.run();
    for (const std::string& name : stats.getNames()) {
        std::cout << name << " ";
    }
    std::cout << std::endl;
    for (double v : stats.getFeatures()) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    return 0;
}
