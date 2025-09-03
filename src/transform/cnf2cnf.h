/**
 * MIT License
 * Copyright (c) 2025 Ashlin Iser 
 */

#pragma once

#include "src/extract/IExtractor.h"
#include "src/util/CNFFormula.h"
 
namespace CNF {
 
class Normaliser : public IExtractor {
    const char* filename_;
    const char* output_;

public:
    Normaliser(const char* filename, const char* output = nullptr) 
        : filename_(filename), output_(output) { }
    virtual ~Normaliser() {}
    virtual void run();
};

class Sanitiser : public IExtractor {
    const char* filename_;
    const char* output_;

public:
    Sanitiser(const char* filename, const char* output = nullptr) 
        : filename_(filename), output_(output) { }
    virtual ~Sanitiser() {}
    virtual void run();
};
 
}  // namespace CNF
