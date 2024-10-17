/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#pragma once

#include "IExtractor.h"

#include "src/util/StreamBuffer.h"

namespace OPB {

class Constr;
class BaseFeatures;

class TermSum {
    friend Constr;
    friend BaseFeatures;

    std::vector<double> coeffs{};
    double max = 0;
    double min = 0;
    double abs_min_coeff = std::numeric_limits<double>::max();
    Var max_var{0};

  public:
    TermSum(StreamBuffer &in);
    
    inline size_t nTerms();
    inline double const minVal();
    inline double const maxVal();
    inline Var const maxVar();
    inline double const minCoeff();
};

class Constr {
    friend BaseFeatures;
    
  public:
    enum Rel { GE, EQ };
    struct Analysis {
        bool tautology : 1;
        bool unsat : 1;
        bool assignment : 1;
        bool clause : 1;
        bool card : 1;
    };

  private:
    TermSum terms;
    Rel rel;
    std::string strbound;
    double bound;

  public:
    Constr(StreamBuffer &in);
    inline Analysis analyse();
    inline Var maxVar();
};

class BaseFeatures : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;

    unsigned n_vars = 0, n_constraints = 0;
    unsigned n_pbs_ge = 0, n_pbs_eq = 0;
    unsigned n_cards_ge = 0, n_cards_eq = 0;
    unsigned n_clauses = 0, n_assignments = 0;
    bool trivially_unsat = false;
    
    unsigned obj_terms = 0;
    double obj_max_val = 0, obj_min_val = 0;
    std::vector<double> obj_coeffs{};

    void load_feature_record();

  public:
    BaseFeatures(const char* filename);
    virtual ~BaseFeatures();
    virtual void extract();
    virtual std::vector<double> getFeatures() const;
    virtual std::vector<std::string> getNames() const;
};

}; // namespace OPB