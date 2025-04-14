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

 public:
    cnf2bip(const char* filename, const char* output = nullptr);
    virtual ~cnf2bip();
    virtual void run();
};

}  // namespace CNF
