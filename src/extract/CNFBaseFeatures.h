/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#pragma once

#include "IExtractor.h"
#include <array>

namespace CNF {

class BaseFeatures : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;

    void extractBaseFeatures1();
    void extractBaseFeatures2();

  public:
    BaseFeatures(const char* filename);
    virtual ~BaseFeatures();
    virtual void extract();
    virtual std::vector<double> getFeatures() const;
    virtual std::vector<std::string> getNames() const;
};

class BaseFeatures1 : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;
    unsigned n_vars = 0, n_clauses = 0, bytes = 0, ccs = 0;
    // count occurences of clauses of small size
    std::array<unsigned, 11> clause_sizes;
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

    void load_feature_record();

  public:
    BaseFeatures1(const char* filename);
    virtual ~BaseFeatures1();
    virtual void extract();
    virtual std::vector<double> getFeatures() const;
    virtual std::vector<std::string> getNames() const;
};

class BaseFeatures2 : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;
    unsigned n_vars = 0, n_clauses = 0;
    // VCG Degree Distribution:
    std::vector<unsigned> vcg_cdegree; // clause sizes
    std::vector<unsigned> vcg_vdegree; // occurence counts
    // VIG Degree Distribution:
    std::vector<unsigned> vg_degree;
    // CG Degree Distribution:
    std::vector<unsigned> clause_degree;

    void load_feature_records();

  public:
    BaseFeatures2(const char* filename);
    virtual ~BaseFeatures2();
    virtual void extract();
    virtual std::vector<double> getFeatures() const;
    virtual std::vector<std::string> getNames() const;
};

}; // namespace CNF