/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#include <cmath>

#include "src/extract/CNFBaseFeatures.h"

#include "src/util/StreamBuffer.h"
#include "src/util/CaptureDistribution.h"
#include "src/util/UnionFind.h"

CNF::BaseFeatures1::BaseFeatures1(const char* filename) : filename_(filename) { 
    clause_sizes.fill(0);
    initFeatures({ "clauses", "variables", "bytes", "ccs" });
    initFeatures({ "cls1", "cls2", "cls3", "cls4", "cls5", "cls6", "cls7", "cls8", "cls9", "cls10p" });
    initFeatures({ "horn", "invhorn", "positive", "negative" });
    initFeatures({ "hornvars_mean", "hornvars_variance", "hornvars_min", "hornvars_max", "hornvars_entropy" });
    initFeatures({ "invhornvars_mean", "invhornvars_variance", "invhornvars_min", "invhornvars_max", "invhornvars_entropy" });
    initFeatures({ "balancecls_mean", "balancecls_variance", "balancecls_min", "balancecls_max", "balancecls_entropy" });
    initFeatures({ "balancevars_mean", "balancevars_variance", "balancevars_min", "balancevars_max", "balancevars_entropy" });
}

CNF::BaseFeatures1::~BaseFeatures1() { }

void CNF::BaseFeatures1::run() {
    StreamBuffer in(filename_);
    UnionFind uf;
    Cl clause;
    while (in.readClause(clause)) {
        ++n_clauses;            
        ++clause_sizes[std::min(clause.size(), 10UL)];
        bytes += 2;

        uf.insert(clause);

        unsigned n_neg = 0;
        for (Lit lit : clause) {
            bytes += lit.sign() + ceil(log10((float)lit.var())) + 1;
            // resize vectors if necessary
            if (static_cast<unsigned>(lit.var()) > n_vars) {
                n_vars = lit.var();
                variable_horn.resize(n_vars + 1);
                variable_inv_horn.resize(n_vars + 1);
                literal_occurrences.resize(2 * n_vars + 2);
            }
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
    }
    // balance of positive and negative literals per variable
    for (unsigned v = 0; v < n_vars; v++) {
        double pos = (double)literal_occurrences[Lit(v, false)];
        double neg = (double)literal_occurrences[Lit(v, true)];
        if (std::max(pos, neg) > 0) {
            balance_variable.push_back(std::min(pos, neg) / std::max(pos, neg));
        }
    }
    ccs = uf.count_components();

    load_feature_record();
}

void CNF::BaseFeatures1::load_feature_record() {
    setFeature("clauses", (double)n_clauses);
    setFeature("variables", (double)n_vars);
    setFeature("bytes", (double)bytes);
    setFeature("ccs", (double)ccs);
    std::vector<double> clause_sizes_double(clause_sizes.begin(), clause_sizes.end());
    setFeatures({ "cls1", "cls2", "cls3", "cls4", "cls5", "cls6", "cls7", "cls8", "cls9", "cls10p" }, clause_sizes_double.begin()+1, clause_sizes_double.end());
    setFeature("horn", (double)horn);
    setFeature("invhorn", (double)inv_horn);
    setFeature("positive", (double)positive);
    setFeature("negative", (double)negative);
    std::vector<double> stats = getDistributionStats(variable_horn);
    setFeatures({ "hornvars_mean", "hornvars_variance", "hornvars_min", "hornvars_max", "hornvars_entropy" }, stats.begin(), stats.end());
    stats = getDistributionStats(variable_inv_horn);
    setFeatures({ "invhornvars_mean", "invhornvars_variance", "invhornvars_min", "invhornvars_max", "invhornvars_entropy" }, stats.begin(), stats.end());
    stats = getDistributionStats(balance_clause);
    setFeatures({ "balancecls_mean", "balancecls_variance", "balancecls_min", "balancecls_max", "balancecls_entropy" }, stats.begin(), stats.end());
    stats = getDistributionStats(balance_variable);
    setFeatures({ "balancevars_mean", "balancevars_variance", "balancevars_min", "balancevars_max", "balancevars_entropy" }, stats.begin(), stats.end());
}

CNF::BaseFeatures2::BaseFeatures2(const char* filename) : filename_(filename) { 
    initFeatures({ "vcg_vdegree_mean", "vcg_vdegree_variance", "vcg_vdegree_min", "vcg_vdegree_max", "vcg_vdegree_entropy" });
    initFeatures({ "vcg_cdegree_mean", "vcg_cdegree_variance", "vcg_cdegree_min", "vcg_cdegree_max", "vcg_cdegree_entropy" });
    initFeatures({ "vg_degree_mean", "vg_degree_variance", "vg_degree_min", "vg_degree_max", "vg_degree_entropy" });
    initFeatures({ "cg_degree_mean", "cg_degree_variance", "cg_degree_min", "cg_degree_max", "cg_degree_entropy" });
}

CNF::BaseFeatures2::~BaseFeatures2() { }

void CNF::BaseFeatures2::run() {
    StreamBuffer in(filename_);

    Cl clause;
    while (in.readClause(clause)) {
        vcg_cdegree.push_back(clause.size());

        for (Lit lit : clause) {
            // resize vectors if necessary
            if (static_cast<unsigned>(lit.var()) > n_vars) {
                n_vars = lit.var();
                vcg_vdegree.resize(n_vars + 1);
                vg_degree.resize(n_vars + 1);
            }
            // count variable occurrences
            ++vcg_vdegree[lit.var()];
            vg_degree[lit.var()] += clause.size();
        }
    }
    // clause graph features
    StreamBuffer in2(filename_);
    while (in2.readClause(clause)) {
        unsigned degree = 0;
        for (Lit lit : clause) {
            degree += vcg_vdegree[lit.var()];
        }
        clause_degree.push_back(degree);
    }

    load_feature_records();
}

void CNF::BaseFeatures2::load_feature_records() {
    std::vector<double> stats = getDistributionStats(vcg_cdegree);
    setFeatures({ "vcg_cdegree_mean", "vcg_cdegree_variance", "vcg_cdegree_min", "vcg_cdegree_max", "vcg_cdegree_entropy" }, stats.begin(), stats.end());
    stats = getDistributionStats(vcg_vdegree);
    setFeatures({ "vcg_vdegree_mean", "vcg_vdegree_variance", "vcg_vdegree_min", "vcg_vdegree_max", "vcg_vdegree_entropy" }, stats.begin(), stats.end());
    stats = getDistributionStats(vg_degree);
    setFeatures({ "vg_degree_mean", "vg_degree_variance", "vg_degree_min", "vg_degree_max", "vg_degree_entropy" }, stats.begin(), stats.end());
    stats = getDistributionStats(clause_degree);
    setFeatures({ "cg_degree_mean", "cg_degree_variance", "cg_degree_min", "cg_degree_max", "cg_degree_entropy" }, stats.begin(), stats.end());
}

CNF::BaseFeatures::BaseFeatures(const char* filename) : filename_(filename) { 
    BaseFeatures1 baseFeatures1(filename_);
    auto names1 = baseFeatures1.getNames();
    initFeatures(names1.begin(), names1.end());
    BaseFeatures2 baseFeatures2(filename_);
    auto names2 = baseFeatures2.getNames();
    initFeatures(names2.begin(), names2.end());
}

CNF::BaseFeatures::~BaseFeatures() { }

void CNF::BaseFeatures::run() {
    extractBaseFeatures1();
    extractBaseFeatures2();
}

void CNF::BaseFeatures::extractBaseFeatures1() {
    BaseFeatures1 baseFeatures1(filename_);
    baseFeatures1.run();
    for (auto name : baseFeatures1.getNames()) {
        setFeature(name, baseFeatures1.getFeature(name));
    }
}

void CNF::BaseFeatures::extractBaseFeatures2() {
    BaseFeatures2 baseFeatures2(filename_);
    baseFeatures2.run();
    for (auto name : baseFeatures2.getNames()) {
        setFeature(name, baseFeatures2.getFeature(name));
    }
}