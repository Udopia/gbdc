/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#pragma once

#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>

#include "src/extract/IExtractor.h"

namespace CNF {

class GateFeatures : public IExtractor {
    const char *filename_;
    std::vector<double> features;
    std::vector<std::string> names;

    unsigned n_vars = 0, n_gates = 0, n_roots = 0;
    unsigned n_none = 0, n_generic = 0, n_mono = 0;
    unsigned n_and = 0, n_or = 0, n_triv = 0, n_equiv = 0, n_full = 0;

    std::vector<unsigned> levels, levels_none, levels_generic, levels_mono;
    std::vector<unsigned> levels_and, levels_or, levels_triv;
    std::vector<unsigned> levels_equiv, levels_full;

    void load_feature_records();

public:
    GateFeatures(const char* filename);
    virtual ~GateFeatures();
    virtual void extract();
    virtual std::vector<double> getFeatures() const;
    virtual std::vector<std::string> getNames() const;
    virtual std::string getRuntimeDesc() const;
};

} // namespace CNF