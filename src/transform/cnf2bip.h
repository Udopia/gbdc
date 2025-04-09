/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#pragma once

#include <vector>

#include "src/extract/IExtractor.h"
#include "src/util/CNFFormula.h"

namespace CNF {

class cnf2bip : public IExtractor {
 private:
    CNFFormula F;
    const char* filename_;
    const char* output_;
    std::vector<double> features;
    std::vector<std::string> names;

 public:
    cnf2bip(const char* filename, const char* output = nullptr);
    virtual ~cnf2bip();
    virtual void run();
    virtual std::vector<double> getFeatures() const;
    virtual std::vector<std::string> getNames() const;
};

}  // namespace CNF
