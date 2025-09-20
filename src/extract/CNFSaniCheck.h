/**
 * MIT License
 * Copyright (c) 2025 Ashlin Iser 
 */

#pragma once

#include "IExtractor.h"
#include "src/util/CNFFormula.h"

namespace CNF {

class SaniCheck : public IExtractor {
    const char* filename_;
    bool sanicheck;
    void checkNormalised();
    void checkSanitised();

public:
SaniCheck(const char* filename, bool sanicheck = false) 
        : filename_(filename), sanicheck(sanicheck) { }
    virtual ~SaniCheck() {}
    virtual void run();
};

} // namespace CNF