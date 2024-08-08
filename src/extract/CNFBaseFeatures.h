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
#include <array>

namespace CNF {

template <template <typename> typename Alloc>
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
    std::vector<unsigned, Alloc<unsigned>> variable_horn, variable_inv_horn;

    // pos-neg literal balance (per clause)
    std::vector<double, Alloc<double>> balance_clause;

    // pos-neg literal balance (per variable)
    std::vector<double, Alloc<double>> balance_variable;

    // Literal Occurrences
    std::vector<unsigned, Alloc<unsigned>> literal_occurrences;

  public:
    BaseFeatures1(const char* filename) : filename_(filename), features(), names() { 
        clause_sizes.fill(0);
        names.insert(names.end(), { "clauses", "variables", "bytes", "ccs" });
        names.insert(names.end(), { "cls1", "cls2", "cls3", "cls4", "cls5", "cls6", "cls7", "cls8", "cls9", "cls10p" });
        names.insert(names.end(), { "horn", "invhorn", "positive", "negative" });
        names.insert(names.end(), { "hornvars_mean", "hornvars_variance", "hornvars_min", "hornvars_max", "hornvars_entropy" });
        names.insert(names.end(), { "invhornvars_mean", "invhornvars_variance", "invhornvars_min", "invhornvars_max", "invhornvars_entropy" });
        names.insert(names.end(), { "balancecls_mean", "balancecls_variance", "balancecls_min", "balancecls_max", "balancecls_entropy" });
        names.insert(names.end(), { "balancevars_mean", "balancevars_variance", "balancevars_min", "balancevars_max", "balancevars_entropy" });
    }

    virtual ~BaseFeatures1() { }

    virtual void extract() {
        StreamBuffer in(filename_);
        UnionFind uf;
        Cl clause;
        while (in.readClause(clause)) {
            ++n_clauses;            
            ++clause_sizes[std::min(clause.size(), 10UL)];
            // +1 for 0 at EOL and +1 for linebreak
            bytes += 2;

            uf.insert(clause);

            unsigned n_neg = 0;
            for (Lit lit : clause) {
                // +1 for whitespace after variables
                bytes += lit.sign() + numDigits(lit.var()) + 1;
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
        // subtract last linebreak
        bytes -= 1;

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

    void load_feature_record() {
        features.insert(features.end(), { (double)n_clauses, (double)n_vars, (double)bytes, (double)ccs });
        for (unsigned i = 1; i < 11; ++i) {
            features.push_back((double)clause_sizes[i]);
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

template <template <typename> typename Alloc>
class BaseFeatures2 : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;

    unsigned n_vars = 0, n_clauses = 0;

    // VCG Degree Distribution
    std::vector<unsigned, Alloc<unsigned>> vcg_cdegree; // clause sizes
    std::vector<unsigned, Alloc<unsigned>> vcg_vdegree; // occurence counts

    // VIG Degree Distribution
    std::vector<unsigned, Alloc<unsigned>> vg_degree;

    // CG Degree Distribution
    std::vector<unsigned, Alloc<unsigned>> clause_degree;

  public:
    BaseFeatures2(const char* filename) : filename_(filename), features(), names() { 
        names.insert(names.end(), { "vcg_vdegree_mean", "vcg_vdegree_variance", "vcg_vdegree_min", "vcg_vdegree_max", "vcg_vdegree_entropy" });
        names.insert(names.end(), { "vcg_cdegree_mean", "vcg_cdegree_variance", "vcg_cdegree_min", "vcg_cdegree_max", "vcg_cdegree_entropy" });
        names.insert(names.end(), { "vg_degree_mean", "vg_degree_variance", "vg_degree_min", "vg_degree_max", "vg_degree_entropy" });
        names.insert(names.end(), { "cg_degree_mean", "cg_degree_variance", "cg_degree_min", "cg_degree_max", "cg_degree_entropy" });
    }
    virtual ~BaseFeatures2() { }

    virtual void extract() {
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

template <template <typename> typename Alloc = std::allocator>
class BaseFeatures : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;

  public:
    BaseFeatures(const char* filename) : filename_(filename), features(), names() { 
        BaseFeatures1<Alloc> baseFeatures1(filename_);
        auto names1 = baseFeatures1.getNames();
        names.insert(names.end(), names1.begin(), names1.end());
        BaseFeatures2<Alloc> baseFeatures2(filename_);
        auto names2 = baseFeatures2.getNames();
        names.insert(names.end(), names2.begin(), names2.end());
    }

    virtual ~BaseFeatures() { }

    virtual void extract() {
        extractBaseFeatures1();
        extractBaseFeatures2();
    }

    void extractBaseFeatures1() {
        BaseFeatures1<Alloc> baseFeatures1(filename_);
        baseFeatures1.extract();
        auto feat = baseFeatures1.getFeatures();
        features.insert(features.end(), feat.begin(), feat.end());
    }

    void extractBaseFeatures2() {
        BaseFeatures2<Alloc> baseFeatures2(filename_);
        baseFeatures2.extract();
        auto feat = baseFeatures2.getFeatures();
        features.insert(features.end(), feat.begin(), feat.end());
    }

    virtual std::vector<double> getFeatures() const {
        return features;
    }

    virtual std::vector<std::string> getNames() const {
        return names;
    }
};

}; // namespace CNF


#endif // BASE_FEATURES_H_