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
#ifndef SRC_FEATURES_CNFSTATSRED_H_
#define SRC_FEATURES_CNFSTATSRED_H_

#include <math.h>

#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>

#include "src/util/SolverTypes.h"
#include "src/util/CNFFormula.h"

#include "src/features/Util.h"

// Calculate Subset of Satzilla Features + Other CNF Stats
// CF. 2004, Nudelmann et al., Understanding Random SAT - Beyond the Clause-to-Variable Ratio
class CNFStatsRed {
    const CNFFormula& formula_;
    std::vector<float> record;

 public:
    unsigned n_vars, n_clauses;

    explicit CNFStatsRed(const CNFFormula& formula) :
     formula_(formula), record(), n_vars(formula.nVars()), n_clauses(formula.nClauses()) {
    }

    void analyze() {
        std::array<unsigned, 10> clause_sizes;  // one entry per clause-size
        clause_sizes.fill(0);

        for (Cl* clause : formula_) {
            if (clause->size() < 10) {
                ++clause_sizes[clause->size()];
            }
        }

        record.push_back(n_clauses);
        record.push_back(n_vars);
        
        for (unsigned i = 1; i < 10; ++i) {
            record.push_back(clause_sizes[i]);
        }
    }

    std::vector<float> BaseFeatures() {
        return record;
    }

    static std::vector<std::string> BaseFeatureNames() {
        return std::vector<std::string> { "clauses", "variables",
            "clause_size_1", "clause_size_2", "clause_size_3", "clause_size_4", "clause_size_5", "clause_size_6", "clause_size_7", "clause_size_8", "clause_size_9"
        };
    }
};

#endif  // SRC_FEATURES_CNFSTATSRED_H_
