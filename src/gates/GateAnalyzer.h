/*************************************************************************************************
CNFTools -- Copyright (c) 2015, Markus Iser, KIT - Karlsruhe Institute of Technology

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

#ifndef SRC_GATES_GATEANALYZER_H_
#define SRC_GATES_GATEANALYZER_H_

#include <cstdlib>
#include <algorithm>
#include <memory>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <climits>

#include "lib/ipasir.h"

#include "src/util/CNFFormula.h"
#include "src/util/ResourceLimits.h"

#include "src/gates/GateFormula.h"
#include "src/gates/BlockList.h"
#include "src/gates/OccurrenceList.h"


// T = BlockList has better root-selection heuristic but is slower in general
// BlockList can not be used atm as I changed the interface in favor of more hard wired code,
// i.e., specifically the "swap from index" in gate_formula.addGate
template<class T = OccurrenceList>
class GateAnalyzer {
    void* S;  // solver

    const CNFFormula& formula_;
    const ResourceLimits& limits_;

    GateFormula gate_formula;

    T index;  // occurence-list

    // analyzer configuration:
    bool patterns = false;
    bool semantic = false;
    unsigned max_ = 1;
    unsigned verbose_ = 0;

 public:
    GateAnalyzer(const CNFFormula& formula, const ResourceLimits& limits, bool patterns_, bool semantic_, unsigned max, unsigned verbose = 0) :
     formula_(formula), limits_(limits), gate_formula(formula.nVars(), verbose), index(formula),
     patterns(patterns_), semantic(semantic_), max_(max), verbose_(verbose) {
        if (semantic) S = ipasir_init();
    }

    ~GateAnalyzer() {
        if (semantic) ipasir_release(S);
    }

    GateFormula getGateFormula() const {
        return gate_formula;
    }

    /**
     * @brief Starting-point gate analysis: iterative root selection
     */
    void analyze() {
        std::vector<Cl*> root_clauses = index.estimateRoots();

        for (unsigned count = 0; count < max_ && !root_clauses.empty(); count++) {
            std::vector<Lit> candidates;
            for (Cl* clause : root_clauses) {
                gate_formula.addRoot(clause);
                candidates.insert(candidates.end(), clause->begin(), clause->end());
            }

            gate_recognition(candidates);

            root_clauses = index.estimateRoots();
        }

        std::unordered_set<Cl*> remainder;
        for (size_t lit = 0; lit < index.size(); lit++) {
            remainder.insert(index[lit].begin(), index[lit].end());
        }
        gate_formula.remainder.insert(gate_formula.remainder.end(), remainder.begin(), remainder.end());
    }

 private:
    /**
     * @brief Start hierarchical gate recognition with given root literals
     * 
     * @param roots 
     */
    void gate_recognition(std::vector<Lit> roots) {
        // std::cerr << "c Starting gate-recognition with roots: " << roots << std::endl;
        std::vector<Lit> candidates { roots.begin(), roots.end() };
        std::unordered_set<Lit> frontier;
        while (!candidates.empty()) {  // breadth_ first search is important here
            // std::cout << "Number of Candidates: " << candidates.size() << std::endl;
            limits_.within_limits_or_throw();
            for (Lit candidate : candidates) {
                if (checkAddGate(candidate)) {
                    Gate& gate = gate_formula.getGate(candidate);
                    index.remove(gate.fwd);
                    index.remove(gate.bwd);
                    frontier.insert(gate.inp.begin(), gate.inp.end());
                }
            }
            // std::cout << "frontier size " << frontier.size() << std::endl;
            candidates.clear();
            candidates.insert(candidates.end(), frontier.begin(), frontier.end());
            frontier.clear();
        }
    }

    std::vector<Lit> getInputLiterals(Lit output, const For& clauses) {
        std::vector<Lit> inp;
        for (Cl* clause : clauses) {
            unsigned pos = 0;  // reset insert position for each clause
            for (auto it = clause->begin(); it != clause->end(); ++it) {
                if (*it != output) {
                    while (pos < inp.size() && inp[pos] < *it) {  // clauses are sorted ;)
                        ++pos;
                    }
                    if (pos == inp.size()) {
                        // append all except for ~out and break
                        for (; *it < output; ++it) {
                            inp.insert(inp.end(), *it);
                        }
                        inp.insert(inp.end(), ++it, clause->end());
                        break;
                    } else if (inp[pos] > *it) {
                        inp.insert(inp.begin() + pos, *it);
                    }  // else: do not insert duplicate
                    ++pos;
                }
            }
        }
        return inp;
    }

    unsigned constrainSameInputVariables(Lit o, const For& fwd, const For& bwd) {
        // check if fwd and bwd constrain exactly the same inputs, return 0 on failure, otherwise return number of input variables
        std::unordered_set<Var> fwd_vars;
        std::unordered_set<Var> bwd_vars;
        for (Cl* c : fwd) for (Lit l : *c) if (l != ~o) fwd_vars.insert(l.var());
        for (Cl* c : bwd) for (Lit l : *c) if (l != o) {
            bool inserted = std::get<1>(bwd_vars.insert(l.var()));
            if (inserted && !fwd_vars.count(l.var())) {  // ensure: bwd_vars \subseteq fwd_vars
                return 0;
            }
        }
        if (fwd_vars.size() > bwd_vars.size()) {  // ensure: fwd_vars \subseteq bwd_vars
            return 0;
        }
        return fwd_vars.size();
    }

    /**
     * @brief checks if index contains a gate definition for the given candidate output and adds gate if positive
     * @return true if clauses encode gate, false otherwise
     */
    bool checkAddGate(Lit out) {
        // std::cout << "check add gate " << out << std::endl;
        if (index[~out].size() > 0 && index.isBlockedSet(out)) {
            GateType type = NONE;

            if (gate_formula.isNestedMonotonic(out)) {
                type = MONO;
            } else if (patterns) {
                unsigned input_size = constrainSameInputVariables(out, index[~out], index[out]);
                if (input_size > 0) {
                    type = fPattern(out, index[~out], index[out], input_size);
                }
            }

            if (type == NONE && semantic) {
                if (index[~out].size() > 1 && index[out].size() > 1) {  // case excluded by patterns
                    type = fSemantic(out, index[~out], index[out]);
                }
            }

            if (type != NONE) {
                gate_formula.addGate(type, out, index[~out], index[out], getInputLiterals(~out, index[~out]));
                return true;
            }
        }
        return false;
    }

    // clause patterns of full encoding
    // precondition: fwd blocks bwd on output literal o
    // fwd and bwd constrain same input variables
    GateType fPattern(Lit o, const For& fwd, const For& bwd, unsigned input_size) {
        // detect or gates
        if (fwd.size() == 1 && fixedClauseSize(bwd, 2)) {
            if (input_size == 1) return TRIV;
            return OR;
        }
        // detect and gates
        if (bwd.size() == 1 && fixedClauseSize(fwd, 2)) {
            return AND;
        }
        // under the preconditions (blocked set, same inputs, absence of redundancy) the follwing holds:
        // 2^n blocked clauses of size n+1 represent all input combinations with an output literal
        // Note that CNFFormula's input sanitizer takes care of inner-clause redundency but does not check for duplicate clauses!
        if (fwd.size() + bwd.size() == std::uint64_t(1) << input_size) {
            if (fixedClauseSize(fwd, input_size+1) && fixedClauseSize(bwd, input_size+1)) {
                if (input_size == 2 && fwd.size() == bwd.size()) return EQIV;
                return FULL;
            }
        }
        return NONE;
    }

    GateType fSemantic(Lit o, const For& fwd, const For& bwd) {
        // std::cout << "Semantic check for " << fwd.size() + bwd.size() << " clauses" << std::endl;
        // std::cout << fwd << std::endl;
        // std::cout << bwd << std::endl;
        for (const For& f : { fwd, bwd }) {
            for (Cl* cl : f) {
                for (Lit lit : *cl) {
                    if (lit.var() != o.var()) {
                        ipasir_add(S, lit.toDimacs());
                    } else {
                        ipasir_add(S, lit.positive().toDimacs());
                    }
                }
                ipasir_add(S, 0);
            }
        }
        ipasir_assume(S, o.negative().toDimacs());
        int result = ipasir_solve(S);
        ipasir_add(S, o.positive().toDimacs());
        return result == 20 ? GENERIC : NONE;
    }

    bool fixedClauseSize(const For& f, unsigned int n) {
        for (Cl* c : f) if (c->size() != n) return false;
        return true;
    }
};

#endif  // SRC_GATES_GATEANALYZER_H_
