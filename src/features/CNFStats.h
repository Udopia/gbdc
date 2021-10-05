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
#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>

#include "src/util/SolverTypes.h"
#include "src/util/CNFFormula.h"
#include "src/util/Runtime.h"

#include "src/features/Util.h"

class CNFStats {
    CNFFormula formula_;
    Runtime runtime;

 public:
    unsigned n_vars, n_clauses, n_literals;
    std::vector<unsigned> clause_sizes;  // one entry per clause-size

    // VCG Node Distribution:
    std::vector<unsigned> clause_occurrences;  // one entry per clause (its size)
    std::vector<unsigned> variable_occurrences;  // one entry per variable (its num. of occurrences)
    // VG Node Distribution
    std::vector<unsigned> variable_degree;  // one entry per variable (its degree in the VIG)
    // CG Node Distribution
    std::vector<unsigned> clause_degree;  // one entry per clause (number of neighbour clauses)

    // Pos-Neg Literals Ratio
    std::vector<float> pos_neg_per_clause;  // one entry per clause
    std::vector<float> pos_neg_per_variable;  // one entry per variable

    // Horn
    unsigned horn;  // number of horn clauses
    unsigned inv_horn;  // number of inv. horn clauses
    unsigned positive, negative;  // number of positive / negative clauses
    std::vector<unsigned> variable_horn;  // occurrences in horn clauses (per variable)
    std::vector<unsigned> variable_inv_horn;  // occurrences in inv. horn clauses (per variable)

    explicit CNFStats(const CNFFormula& formula) :
     formula_(formula), n_vars(formula.nVars()), n_clauses(formula.nClauses()), n_literals(0),
     clause_sizes(), clause_occurrences(), variable_occurrences(),
     variable_degree(), clause_degree(),
     pos_neg_per_clause(), pos_neg_per_variable(),
     horn(0), inv_horn(0), positive(0), negative(0),
     variable_horn(), variable_inv_horn() {
        clause_sizes.resize(formula.nVars(), 0);
        variable_occurrences.resize(formula.nVars() + 1);
        variable_degree.resize(formula.nVars() + 1);
        clause_degree.resize(formula.nClauses(), 0);
        variable_horn.resize(formula.nVars() + 1);
        variable_inv_horn.resize(formula.nVars() + 1);
    }

    void analyze() {
        runtime.start();
        std::vector<unsigned> literal_occurrences;
        literal_occurrences.resize(2 * formula_.nVars() + 2);
        for (Cl* clause : formula_) {
            n_literals += clause->size();
            ++clause_sizes[clause->size()];
            clause_occurrences.push_back(clause->size());
            float neg = 0;
            for (Lit lit : *clause) {
                ++variable_occurrences[lit.var()];
                ++literal_occurrences[lit];
                variable_degree[lit.var()] += clause->size() - 1;
                // variable_jw_degree[lit.var()] += clause->size() / pow(2, clause->size());
                if (lit.sign()) ++neg;
            }
            // divide min by max (not pos by neg as in satzilla)
            float pos = clause->size() - neg;
            pos_neg_per_clause.push_back(std::max(pos, neg) > 0 ? std::min(pos, neg) / std::max(pos, neg) : 0);
            // horn
            if (neg <= 1) {
                if (neg == 0) ++positive;
                ++horn;
                for (Lit lit : *clause) {
                    ++variable_horn[lit.var()];
                }
            }
            if (pos <= 1) {
                if (pos == 0) ++negative;
                ++inv_horn;
                for (Lit lit : *clause) {
                    ++variable_inv_horn[lit.var()];
                }
            }
        }
        for (unsigned v = 0; v < n_vars; v++) {
            // divide min by max (not pos by neg as in satzilla)
            float pos = static_cast<float>(literal_occurrences[Lit(v, false)]);
            float neg = static_cast<float>(literal_occurrences[Lit(v, true)]);
            pos_neg_per_variable.push_back(std::max(pos, neg) > 0 ? std::min(pos, neg) / std::max(pos, neg) : 0);
        }
        unsigned cid = 0;
        for (Cl* clause : formula_) {
            clause_degree[cid] = std::accumulate(clause->begin(), clause->end(), 0,
                [this] (unsigned a, Lit lit) { return a + variable_occurrences[lit.var()]; } );
            clause_degree[cid] -= clause->size();
            ++cid;
        }
        runtime.stop();
    }

    // Subset of Satzilla Features + Other CNF Stats
    // CF. 2004, Nudelmann et al., Understanding Random SAT - Beyond the Clause-to-Variable Ratio
    std::vector<float> BaseFeatures() {
        std::vector<float> record;

        // ## Problem Size Features ##
        record.push_back(n_clauses);
        record.push_back(n_vars);
        // DEL Ratio: C/V (including squared and cubic variants)
        // DEL Reciprocal Ratio: V/C (including squared and cubic variants)
        // DEL Linearized Ratio: |4.26 - C/V| (including squared and cubic variants)

        // changed to absolute numbers, not fraction of unary, binary and ternary clauses, also added more
        record.push_back(clause_sizes[1]);
        record.push_back(clause_sizes[2]);
        record.push_back(clause_sizes[3]);
        record.push_back(clause_sizes[4]);
        record.push_back(clause_sizes[5]);
        record.push_back(clause_sizes[6]);
        record.push_back(clause_sizes[7]);
        record.push_back(clause_sizes[8]);
        record.push_back(clause_sizes[9]);

        // ## Horn, Pos, Neg Proximity ##
        record.push_back(horn);
        record.push_back(inv_horn);
        record.push_back(positive);
        record.push_back(negative);
        push_distribution(&record, variable_horn);
        push_distribution(&record, variable_inv_horn);

        // ## Balance Features ##
        push_distribution(&record, pos_neg_per_clause);  // min over max
        push_distribution(&record, pos_neg_per_variable);  // min over max

        // ## Variable - Clause Graph Features ##
        // Variable Node Degree Statistics:
        push_distribution(&record, variable_occurrences);
        // Clause Node Degree Statistics:
        push_distribution(&record, clause_occurrences);

        // ## Variable Graph Features ##
        push_distribution(&record, variable_degree);

        // ## Clause Graph Features ##
        push_distribution(&record, clause_degree);
        // DEL Clustering Coefficient Statistics (FW: community structure features)

        // ## Missing: LP-Based Features, DPLL Search Space, Local Search Probes

        record.push_back(static_cast<float>(runtime.get()));

        return record;
    }

    static std::vector<std::string> BaseFeatureNames() {
        return std::vector<std::string> { "clauses", "variables",
            "clause_size_1", "clause_size_2", "clause_size_3", "clause_size_4", "clause_size_5", "clause_size_6", "clause_size_7", "clause_size_8", "clause_size_9",
            "horn_clauses", "inv_horn_clauses", "positive_clauses", "negative_clauses",
            "horn_vars_mean", "horn_vars_variance", "horn_vars_min", "horn_vars_max", "horn_vars_entropy",
            "inv_horn_vars_mean", "inv_horn_vars_variance", "inv_horn_vars_min", "inv_horn_vars_max", "inv_horn_vars_entropy",
            "balance_clause_mean", "balance_clause_variance", "balance_clause_min", "balance_clause_max", "balance_clause_entropy",
            "balance_vars_mean", "balance_vars_variance", "balance_vars_min", "balance_vars_max", "balance_vars_entropy",
            "vcg_vdegrees_mean", "vcg_vdegrees_variance", "vcg_vdegrees_min", "vcg_vdegrees_max", "vcg_vdegrees_entropy",
            "vcg_cdegrees_mean", "vcg_cdegrees_variance", "vcg_cdegrees_min", "vcg_cdegrees_max", "vcg_cdegrees_entropy",
            "vg_degrees_mean", "vg_degrees_variance", "vg_degrees_min", "vg_degrees_max", "vg_degrees_entropy",
            "cg_degrees_mean", "cg_degrees_variance", "cg_degrees_min", "cg_degrees_max", "cg_degrees_entropy", "base_features_runtime"
        };
    }
};

#endif  // SRC_FEATURES_CNFSTATS_H_
