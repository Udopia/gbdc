/**
 * MIT License
 * Copyright (c) 2025 Ashlin Iser, Christoph Jabs
 */

#pragma once

#include "IExtractor.h"
#include "src/util/StreamBuffer.h"
#include "src/util/CaptureDistribution.h"
#include <array>

namespace WCNF {

class BaseFeatures1 : public IExtractor {
    const char* filename_;

    unsigned n_vars = 0, n_hard_clauses = 0, n_soft_clauses = 0;
    uint64_t weight_sum = 0;
    // count occurences of hard clauses of small size
    std::array<unsigned, 11> hard_clause_sizes;
    // count occurences of soft clauses of small size
    std::array<unsigned, 11> soft_clause_sizes;
    // numbers of (inverted) horn clauses
    unsigned horn = 0, inv_horn = 0;
    // number of positive and negative clauses
    unsigned positive = 0, negative = 0;
    // occurrence counts in horn clauses (per variable)
    std::vector<unsigned> variable_horn, variable_inv_horn;
    // pos-neg literal balance (per clause)
    std::vector<double> balance_clause;
    // pos-neg literal balance (per variable)
    std::vector<double> balance_variable;
    // Literal Occurrences
    std::vector<unsigned> literal_occurrences;    
    // Soft clause weights
    std::vector<uint64_t> weights;

    void load_feature_record();

  public:
    BaseFeatures1(const char* filename);
    virtual ~BaseFeatures1();
    virtual void run();
};

class BaseFeatures2 : public IExtractor {
    const char* filename_;

    unsigned n_vars = 0;
    // VCG Degree Distribution
    std::vector<unsigned> vcg_cdegree; // clause sizes
    std::vector<unsigned> vcg_vdegree; // occurence counts
    // VIG Degree Distribution
    std::vector<unsigned> vg_degree;
    // CG Degree Distribution
    std::vector<unsigned> clause_degree;

    void load_feature_records();

  public:
    BaseFeatures2(const char* filename);
    virtual ~BaseFeatures2();
    virtual void run();
};

class BaseFeatures : public IExtractor {
    const char* filename_;

    void extractBaseFeatures1();
    void extractBaseFeatures2();

  public:
    BaseFeatures(const char* filename);
    virtual ~BaseFeatures();
    virtual void run();
};

}; // namespace WCNF