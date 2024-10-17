/**
 * MIT License
 * Copyright (c) 2024 Markus Iser, Christoph Jabs
 */

#include "WCNFBaseFeatures.h"

#include <cassert>
#include <algorithm>

WCNF::BaseFeatures1::BaseFeatures1(const char* filename) : filename_(filename), features(), names() { 
    hard_clause_sizes.fill(0);
    soft_clause_sizes.fill(0);
    names.insert(names.end(), { "h_clauses", "variables" });
    names.insert(names.end(), { "h_cls1", "h_cls2", "h_cls3", "h_cls4", "h_cls5", "h_cls6", "h_cls7", "h_cls8", "h_cls9", "h_cls10p" });
    names.insert(names.end(), { "h_horn", "h_invhorn", "h_positive", "h_negative" });
    names.insert(names.end(), { "h_hornvars_mean", "h_hornvars_variance", "h_hornvars_min", "h_hornvars_max", "h_hornvars_entropy" });
    names.insert(names.end(), { "h_invhornvars_mean", "h_invhornvars_variance", "h_invhornvars_min", "h_invhornvars_max", "h_invhornvars_entropy" });
    names.insert(names.end(), { "h_balancecls_mean", "h_balancecls_variance", "h_balancecls_min", "h_balancecls_max", "h_balancecls_entropy" });
    names.insert(names.end(), { "h_balancevars_mean", "h_balancevars_variance", "h_balancevars_min", "h_balancevars_max", "h_balancevars_entropy" });
    names.insert(names.end(), { "s_clauses", "s_weight_sum" });
    names.insert(names.end(), { "s_cls1", "s_cls2", "s_cls3", "s_cls4", "s_cls5", "s_cls6", "s_cls7", "s_cls8", "s_cls9", "s_cls10p" });
    names.insert(names.end(), { "s_weight_mean", "s_weight_variance", "s_weight_min", "s_weight_max", "s_weight_entropy" });
}

WCNF::BaseFeatures1::~BaseFeatures1() { }

void WCNF::BaseFeatures1::extract() {
    StreamBuffer in(filename_);

    Cl clause;
    uint64_t top = 0; // if top is 0, parsing new file format
    uint64_t weight = 0; // if weight is 0, parsing hard clause
    while (in.skipWhitespace()) {
        if (*in == 'c') {
            if (!in.skipLine()) break;
            continue;
        } else if (*in == 'p') {
            // old format: extract top
            in.skip();
            in.skipWhitespace();
            in.skipString("wcnf");
            // skip vars
            in.skipNumber();
            // skip clauses
            in.skipNumber();
            // extract top
            in.readUInt64(&top);
            in.skipLine();
            continue;
        } else if (*in == 'h') {
            assert(top == 0);
            weight = 0;
            
            in.skip();
            in.readClause(clause);
        } else {
            in.readUInt64(&weight);
            if (top > 0 && weight >= top) {
                // old hard clause
                weight = 0;
            }
            in.readClause(clause);
        }
        
        for (Lit lit : clause) {
            if (static_cast<unsigned>(lit.var()) > n_vars) {
                n_vars = lit.var();
                variable_horn.resize(n_vars + 1);
                variable_inv_horn.resize(n_vars + 1);
                literal_occurrences.resize(2 * n_vars + 2);
            }
        }

        // record statistics
        if (!weight) {
            ++n_hard_clauses;
            
            if (clause.size() < 10) {
                ++hard_clause_sizes[clause.size()];
            } else {
                ++hard_clause_sizes[10];
            }

            unsigned n_neg = 0;
            for (Lit lit : clause) {
                // count negative literals
                if (lit.sign()) ++n_neg;
                ++literal_occurrences[lit];
            }

            // horn statistics
            unsigned n_pos = clause.size() - n_neg;
            if (n_neg <= 1) {
                if (n_neg == 0) ++positive;
                ++horn;
                for (Lit lit : clause) {
                    ++variable_horn[lit.var()];
                }
            }
            if (n_pos <= 1) {
                if (n_pos == 0) ++negative;
                ++inv_horn;
                for (Lit lit : clause) {
                    ++variable_inv_horn[lit.var()];
                }
            }

            // balance of positive and negative literals per clause
            if (clause.size() > 0) {
                balance_clause.push_back((double)std::min(n_pos, n_neg) / (double)std::max(n_pos, n_neg));
            }
        } else {
            ++n_soft_clauses;
            weight_sum += weight;

            if (clause.size() < 10) {
                ++soft_clause_sizes[clause.size()];
            } else {
                ++soft_clause_sizes[10];
            }

            weights.push_back(weight);
        }
    }

    // balance of positive and negative literals per variable
    for (unsigned v = 0; v < n_vars; v++) {
        double pos = (double)literal_occurrences[Lit(v, false)];
        double neg = (double)literal_occurrences[Lit(v, true)];
        if (std::max(pos, neg) > 0) {
            balance_variable.push_back(std::min(pos, neg) / std::max(pos, neg));
        }
    }

    load_feature_record();
}

void WCNF::BaseFeatures1::load_feature_record() {
    features.insert(features.end(), { (double)n_hard_clauses, (double)n_vars });
    for (unsigned i = 1; i < 11; ++i) {
        features.push_back((double)hard_clause_sizes[i]);
    }
    features.insert(features.end(), { (double)horn, (double)inv_horn, (double)positive, (double)negative });
    push_distribution(features, variable_horn);
    push_distribution(features, variable_inv_horn);
    push_distribution(features, balance_clause);
    push_distribution(features, balance_variable);
    features.insert(features.end(), { (double)n_soft_clauses, (double)weight_sum });
    for (unsigned i = 1; i < 11; ++i) {
        features.push_back((double)soft_clause_sizes[i]);
    }
    push_distribution(features, weights);
}

std::vector<double> WCNF::BaseFeatures1::getFeatures() const {
    return features;
}

std::vector<std::string> WCNF::BaseFeatures1::getNames() const {
    return names;
}

WCNF::BaseFeatures2::BaseFeatures2(const char* filename) : filename_(filename), features(), names() { 
    names.insert(names.end(), { "h_vcg_vdegree_mean", "h_vcg_vdegree_variance", "h_vcg_vdegree_min", "h_vcg_vdegree_max", "h_vcg_vdegree_entropy" });
    names.insert(names.end(), { "h_vcg_cdegree_mean", "h_vcg_cdegree_variance", "h_vcg_cdegree_min", "h_vcg_cdegree_max", "h_vcg_cdegree_entropy" });
    names.insert(names.end(), { "h_vg_degree_mean", "h_vg_degree_variance", "h_vg_degree_min", "h_vg_degree_max", "h_vg_degree_entropy" });
    names.insert(names.end(), { "h_cg_degree_mean", "h_cg_degree_variance", "h_cg_degree_min", "h_cg_degree_max", "h_cg_degree_entropy" });
}

WCNF::BaseFeatures2::~BaseFeatures2() { }

void WCNF::BaseFeatures2::extract() {
    StreamBuffer in(filename_);

    Cl clause;
    uint64_t top = 0; // if top is 0, parsing new file format
    uint64_t weight;
    while (in.skipWhitespace()) {
        if (*in == 'c') {
            if (!in.skipLine()) break;
            continue;
        } else if (*in == 'p') {
            // old format: extract top
            in.skip();
            in.skipWhitespace();
            in.skipString("wcnf");
            // skip vars
            in.skipNumber();
            // skip clauses
            in.skipNumber();
            // extract top
            in.readUInt64(&top);
            in.skipLine();
            continue;
        } else if (*in == 'h') {
            assert(top == 0);
            in.skip();
            in.readClause(clause);
            weight = 0;
        } else {
            in.readUInt64(&weight);
            in.readClause(clause);
            // don't skip soft clause here since we need the true variable count
        }
        
        vcg_cdegree.push_back(clause.size());

        for (Lit lit : clause) {
            // resize vectors if necessary
            if (static_cast<unsigned>(lit.var()) > n_vars) {
                n_vars = lit.var();
                vcg_vdegree.resize(n_vars + 1);
                vg_degree.resize(n_vars + 1);
            }

            // count variable occurrences (only for hard clauses)
            if ((!top && !weight) || weight >= top) {
                ++vcg_vdegree[lit.var()];
                vg_degree[lit.var()] += clause.size();
            }
        }
    }

    // clause graph features
    StreamBuffer in2(filename_);
    while (in2.skipWhitespace()) {
        if (*in2 == 'c' || *in2 == 'p') {
            if (!in2.skipLine()) break;
            continue;
        } else if (*in2 == 'h') {
            assert(top == 0);
            in2.skip();
            in2.readClause(clause);
        } else {
            in2.readUInt64(&weight);
            in2.readClause(clause);
            // skip soft clauses
            if (!top || weight < top) continue;
        }

        unsigned degree = 0;
        for (Lit lit : clause) {
            degree += vcg_vdegree[lit.var()];
        }
        clause_degree.push_back(degree);
    }

    load_feature_records();
}

void WCNF::BaseFeatures2::load_feature_records() {
    push_distribution(features, vcg_vdegree);
    push_distribution(features, vcg_cdegree);
    push_distribution(features, vg_degree);
    push_distribution(features, clause_degree);
}

std::vector<double> WCNF::BaseFeatures2::getFeatures() const {
    return features;
}

std::vector<std::string> WCNF::BaseFeatures2::getNames() const {
    return names;
}

WCNF::BaseFeatures::BaseFeatures(const char* filename) : filename_(filename), features(), names() { 
    BaseFeatures1 baseFeatures1(filename_);
    auto names1 = baseFeatures1.getNames();
    names.insert(names.end(), names1.begin(), names1.end());
    BaseFeatures2 baseFeatures2(filename_);
    auto names2 = baseFeatures2.getNames();
    names.insert(names.end(), names2.begin(), names2.end());
}

WCNF::BaseFeatures::~BaseFeatures() { }

void WCNF::BaseFeatures::extract() {
    extractBaseFeatures1();
    extractBaseFeatures2();
}

void WCNF::BaseFeatures::extractBaseFeatures1() {
    BaseFeatures1 baseFeatures1(filename_);
    baseFeatures1.extract();
    auto feat = baseFeatures1.getFeatures();
    features.insert(features.end(), feat.begin(), feat.end());
}

void WCNF::BaseFeatures::extractBaseFeatures2() {
    BaseFeatures2 baseFeatures2(filename_);
    baseFeatures2.extract();
    auto feat = baseFeatures2.getFeatures();
    features.insert(features.end(), feat.begin(), feat.end());
}

std::vector<double> WCNF::BaseFeatures::getFeatures() const {
    return features;
}

std::vector<std::string> WCNF::BaseFeatures::getNames() const {
    return names;
}