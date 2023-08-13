/**
 * MIT License
 * 
 * Copyrigth (c) 2023 Markus Iser 
 */

#ifndef BASE_FEATURES_H_
#define BASE_FEATURES_H_

#include "IExtractor.h"
#include "src/util/StreamBuffer.h"
#include "src/extract/Util.h"

class CNFBaseFeatures1 : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;

    unsigned n_vars = 0, n_clauses = 0;

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

  public:
    CNFBaseFeatures1(const char* filename) : filename_(filename), features(), names() { 
        clause_sizes.fill(0);
        names.insert(names.end(), { "clauses", "variables" });
        names.insert(names.end(), { "cls1", "cls2", "cls3", "cls4", "cls5", "cls6", "cls7", "cls8", "cls9", "cls10", "cls10p" });
        names.insert(names.end(), { "horn", "invhorn", "positive", "negative" });
        names.insert(names.end(), { "hornvars_mean", "hornvars_variance", "hornvars_min", "hornvars_max", "hornvars_entropy" });
        names.insert(names.end(), { "invhornvars_mean", "invhornvars_variance", "invhornvars_min", "invhornvars_max", "invhornvars_entropy" });
        names.insert(names.end(), { "balancecls_mean", "balancecls_variance", "balancecls_min", "balancecls_max", "balancecls_entropy" });
        names.insert(names.end(), { "balancevars_mean", "balancevars_variance", "balancevars_min", "balancevars_max", "balancevars_entropy" });
    }

    virtual ~CNFBaseFeatures1() { }

    virtual void extract() {
        StreamBuffer in(filename_);

        Cl clause;
        while (in.readClause(clause)) {
            ++n_clauses;

            if (clause.size() < 10) {
                ++clause_sizes[clause.size()];
            } else {
                ++clause_sizes[10];
            }

            unsigned n_neg = 0;
            for (Lit lit : clause) {
                // resize vectors if necessary
                if (lit.var() > n_vars) {
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

        load_feature_record();
    }

    void load_feature_record() {
        features.insert(features.end(), { (double)n_clauses, (double)n_vars });
        for (unsigned i = 1; i < 11; ++i) {
            features.push_back(clause_sizes[i]);
        }
        features.insert(features.end(), { (double)horn, (double)inv_horn, (double)positive, (double)negative });
        push_distribution(features, variable_horn);
        push_distribution(features, variable_inv_horn);
        push_distribution(features, balance_clause);
        push_distribution(features, balance_variable);
    }

    virtual std::vector<double> getFeatures() const {
        return features;
    }

    virtual std::vector<std::string> getNames() const {
        return names;
    }
};

class CNFBaseFeatures2 : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;

    unsigned n_vars = 0, n_clauses = 0;

    // VCG Degree Distribution
    std::vector<unsigned> vcg_cdegree; // clause sizes
    std::vector<unsigned> vcg_vdegree; // occurence counts

    // VIG Degree Distribution
    std::vector<unsigned> vg_degree;

    // CG Degree Distribution
    std::vector<unsigned> clause_degree;

  public:
    CNFBaseFeatures2(const char* filename) : filename_(filename), features(), names() { 
        names.insert(names.end(), { "vcg_vdegree_mean", "vcg_vdegree_variance", "vcg_vdegree_min", "vcg_vdegree_max", "vcg_vdegree_entropy" });
        names.insert(names.end(), { "vcg_cdegree_mean", "vcg_cdegree_variance", "vcg_cdegree_min", "vcg_cdegree_max", "vcg_cdegree_entropy" });
        names.insert(names.end(), { "vg_degree_mean", "vg_degree_variance", "vg_degree_min", "vg_degree_max", "vg_degree_entropy" });
        names.insert(names.end(), { "cg_degree_mean", "cg_degree_variance", "cg_degree_min", "cg_degree_max", "cg_degree_entropy" });
    }
    virtual ~CNFBaseFeatures2() { }

    virtual void extract() {
        StreamBuffer in(filename_);

        Cl clause;
        while (in.readClause(clause)) {
            vcg_cdegree.push_back(clause.size());

            for (Lit lit : clause) {
                // resize vectors if necessary
                if (lit.var() > n_vars) {
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

    void load_feature_records() {
        push_distribution(features, vcg_vdegree);
        push_distribution(features, vcg_cdegree);
        push_distribution(features, vg_degree);
        push_distribution(features, clause_degree);
    }

    virtual std::vector<double> getFeatures() const {
        return features;
    }

    virtual std::vector<std::string> getNames() const {
        return names;
    }
};

class CNFBaseFeatures : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;

  public:
    CNFBaseFeatures(const char* filename) : filename_(filename), features(), names() { 
        CNFBaseFeatures1 baseFeatures1(filename_);
        std::vector<std::string> names1 = baseFeatures1.getNames();
        names.insert(names.end(), names1.begin(), names1.end());
        CNFBaseFeatures1 baseFeatures2(filename_);
        std::vector<std::string> names2 = baseFeatures2.getNames();
        names.insert(names.end(), names2.begin(), names2.end());
    }

    virtual ~CNFBaseFeatures() { }

    virtual void extract() {
        extractCNFBaseFeatures1();
        extractCNFBaseFeatures2();
    }

    void extractCNFBaseFeatures1() {
        CNFBaseFeatures1 baseFeatures1(filename_);
        baseFeatures1.extract();
        std::vector<double> feat = baseFeatures1.getFeatures();
        features.insert(features.end(), feat.begin(), feat.end());
    }

    void extractCNFBaseFeatures2() {
        CNFBaseFeatures2 baseFeatures2(filename_);
        baseFeatures2.extract();
        std::vector<double> feat = baseFeatures2.getFeatures();
        features.insert(features.end(), feat.begin(), feat.end());
    }

    virtual std::vector<double> getFeatures() const {
        return features;
    }

    virtual std::vector<std::string> getNames() const {
        return names;
    }
};



#endif // BASE_FEATURES_H_