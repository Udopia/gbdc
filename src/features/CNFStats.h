/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Markus Iser, KIT - Karlsruhe Institute of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/
#ifndef SRC_FEATURES_CNFSTATS_H_
#define SRC_FEATURES_CNFSTATS_H_

#include <math.h>

#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>

#include "src/util/SolverTypes.h"
#include "src/util/CNFFormula.h"
#include "src/util/ResourceLimits.h"

#include "src/features/Util.h"

// Calculate Subset of Satzilla Features + Other CNF Stats
// CF. 2004, Nudelmann et al., Understanding Random SAT - Beyond the Clause-to-Variable Ratio
class CNFStats {
    const CNFFormula& formula_;
    const ResourceLimits& limits_;
    std::vector<float> record;

 public:
    unsigned n_vars, n_clauses;

    explicit CNFStats(const CNFFormula& formula, const ResourceLimits& limits) :
     formula_(formula), limits_(limits), record(), n_vars(formula.nVars()), n_clauses(formula.nClauses()) {
    }

    void analyze_occurrences() {
        // Size
        std::array<unsigned, 10> clause_sizes;  // one entry per clause-size

        // Horn
        unsigned horn = 0, inv_horn = 0;  // number of (inv.) horn clauses
        unsigned positive = 0, negative = 0;  // number of positive / negative clauses
        std::vector<unsigned> variable_horn, variable_inv_horn;  // occurrences in (inv.) horn clauses (per variable)
        variable_horn.resize(n_vars + 1);
        variable_inv_horn.resize(n_vars + 1);

        // VCG Node Distribution:
        std::vector<unsigned> clause_occurrences;  // one entry per clause (its size)
        std::vector<unsigned> variable_occurrences;  // one entry per variable (its num. of occurrences)
        variable_occurrences.resize(n_vars + 1);

        for (Cl* clause : formula_) {
            if (clause->size() < 10) {
                ++clause_sizes[clause->size()];
            }

            clause_occurrences.push_back(clause->size());

            float neg = 0;
            for (Lit lit : *clause) {
                ++variable_occurrences[lit.var()];
                if (lit.sign()) ++neg;
            }

            // horn
            if (neg <= 1) {
                if (neg == 0) ++positive;
                ++horn;
                for (Lit lit : *clause) {
                    ++variable_horn[lit.var()];
                }
            }
            if (clause->size() - neg <= 1) {
                if (clause->size() - neg == 0) ++negative;
                ++inv_horn;
                for (Lit lit : *clause) {
                    ++variable_inv_horn[lit.var()];
                }
            }
        }

        // ## Problem Size Features ##
        record.push_back(n_clauses);
        record.push_back(n_vars);
        // DEL Ratio: C/V (including squared and cubic variants)
        // DEL Reciprocal Ratio: V/C (including squared and cubic variants)
        // DEL Linearized Ratio: |4.26 - C/V| (including squared and cubic variants)

        // changed to absolute numbers, not fraction of unary, binary and ternary clauses, also added more
        for (unsigned i = 1; i < 10; ++i) {
            record.push_back(clause_sizes[i]);
        }

        // ## Horn, Pos, Neg Proximity ##
        record.push_back(horn);
        record.push_back(inv_horn);
        record.push_back(positive);
        record.push_back(negative);
        push_distribution(&record, variable_horn);
        push_distribution(&record, variable_inv_horn);

        // ## Variable - Clause Graph Features ##
        // Variable Node Degree Statistics:
        push_distribution(&record, variable_occurrences);
        // Clause Node Degree Statistics:
        push_distribution(&record, clause_occurrences);
    }

    void analyze_degrees() {
        std::vector<unsigned> literal_occurrences;
        literal_occurrences.resize(2 * formula_.nVars() + 2);

        // VG Node Distribution
        std::vector<unsigned> variable_degree;  // one entry per variable (its degree in the VIG)
        variable_degree.resize(n_vars + 1, 0);

        // min(pos, neg) / max(pos, neg) literal occurence per clause
        std::vector<float> pos_neg_per_clause;  // one entry per clause

        for (Cl* clause : formula_) {
            float neg = 0;
            for (Lit lit : *clause) {
                ++literal_occurrences[lit];
                variable_degree[lit.var()] += 1.0 / pow(2, clause->size());
                if (lit.sign()) ++neg;
            }
            float pos = clause->size() - neg;
            pos_neg_per_clause.push_back(std::max(pos, neg) > 0 ? std::min(pos, neg) / std::max(pos, neg) : 0);
        }
        // ## Variable Graph Features ##
        push_distribution(&record, variable_degree);
        variable_degree.clear();
        std::vector<unsigned>().swap(variable_degree);

        push_distribution(&record, pos_neg_per_clause);  // min over max
        pos_neg_per_clause.clear();
        std::vector<float>().swap(pos_neg_per_clause);

        limits_.within_limits_or_throw();

        std::cout << "Pos/Neg per Variable" << std::endl;
        std::vector<float> pos_neg_per_variable;  // one entry per variable
        for (unsigned v = 0; v < n_vars; v++) {
            // divide min by max (not pos by neg as in satzilla)
            float pos = static_cast<float>(literal_occurrences[Lit(v, false)]);
            float neg = static_cast<float>(literal_occurrences[Lit(v, true)]);
            pos_neg_per_variable.push_back(std::max(pos, neg) > 0 ? std::min(pos, neg) / std::max(pos, neg) : 0);
        }
        push_distribution(&record, pos_neg_per_variable);  // min over max
        pos_neg_per_variable.clear();
        std::vector<float>().swap(pos_neg_per_variable);

        limits_.within_limits_or_throw();

        // ## Clause Graph Features ##
        std::cout << "Clause Graph Features" << std::endl;
        std::vector<unsigned> clause_degree;  // one entry per clause (number of neighbour clauses)
        clause_degree.resize(n_clauses, 0);
        unsigned cid = 0;
        for (Cl* clause : formula_) {
            clause_degree[cid] = 0;
            for (Lit lit : *clause) {
                clause_degree[cid] += literal_occurrences[~lit];
            }
            ++cid;
        }
        push_distribution(&record, clause_degree);
    }

    void analyze() {
        limits_.within_limits_or_throw();
        std::cout << "Analyzing Occurrences" << std::endl;
        analyze_occurrences();
        limits_.within_limits_or_throw();
        std::cout << "Analyzing Degrees" << std::endl;
        analyze_degrees();
        std::cout << "Done" << std::endl;

        // DEL Clustering Coefficient Statistics (FW: community structure features)

        // ## Missing: LP-Based Features, DPLL Search Space, Local Search Probes

        record.push_back(static_cast<float>(limits_.get_runtime()));
    }

    std::vector<float> BaseFeatures() {
        return record;
    }

    static std::vector<std::string> BaseFeatureNames() {
        return std::vector<std::string> { "clauses", "variables",
            "clause_size_1", "clause_size_2", "clause_size_3", "clause_size_4", "clause_size_5", "clause_size_6", "clause_size_7", "clause_size_8", "clause_size_9",
            "horn_clauses", "inv_horn_clauses", "positive_clauses", "negative_clauses",
            "horn_vars_mean", "horn_vars_variance", "horn_vars_min", "horn_vars_max", "horn_vars_entropy",
            "inv_horn_vars_mean", "inv_horn_vars_variance", "inv_horn_vars_min", "inv_horn_vars_max", "inv_horn_vars_entropy",
            "vg_degrees_mean", "vg_degrees_variance", "vg_degrees_min", "vg_degrees_max", "vg_degrees_entropy",
            "balance_clause_mean", "balance_clause_variance", "balance_clause_min", "balance_clause_max", "balance_clause_entropy",
            "balance_vars_mean", "balance_vars_variance", "balance_vars_min", "balance_vars_max", "balance_vars_entropy",
            "vcg_vdegrees_mean", "vcg_vdegrees_variance", "vcg_vdegrees_min", "vcg_vdegrees_max", "vcg_vdegrees_entropy",
            "vcg_cdegrees_mean", "vcg_cdegrees_variance", "vcg_cdegrees_min", "vcg_cdegrees_max", "vcg_cdegrees_entropy",
            "cg_degrees_mean", "cg_degrees_variance", "cg_degrees_min", "cg_degrees_max", "cg_degrees_entropy", "base_features_runtime"
        };
    }
};

#endif  // SRC_FEATURES_CNFSTATS_H_
