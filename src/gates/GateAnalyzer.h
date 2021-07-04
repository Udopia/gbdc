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

#ifndef GATE_ANALYZER
#define GATE_ANALYZER

#include <cstdlib>
#include <algorithm>
#include <memory>
#include <cmath>
#include <vector>
#include <set>
#include <climits>

#include "util/CNFFormula.h"
#include "util/Runtime.h"

#include "gates/GateFormula.h"
#include "gates/BlockList.h"
#include "gates/OccurrenceList.h"

#include "ipasir.h"

template<class T = OccurrenceList>
class GateAnalyzer {
    
    void* S; // solver

    const CNFFormula& problem;
    GateFormula gate_formula;

    T index; // occurence-list

    // analyzer configuration:
    bool patterns = false;
    bool semantic = false;
    unsigned int max = 1;

public:
    GateAnalyzer(const CNFFormula& problem_, bool patterns_, bool semantic_, int tries_) :
        problem(problem_), gate_formula(problem_.nVars()), index(problem_), 
        patterns (patterns_), semantic (semantic_), max (tries_)
    {
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

        for (unsigned count = 0; count < max && !root_clauses.empty(); count++) {
            std::vector<Lit> candidates;
            for (Cl* clause : root_clauses) {            
                gate_formula.roots.push_back(clause);
                candidates.insert(candidates.end(), clause->begin(), clause->end());
                for (Lit l : *clause) gate_formula.setUsedAsInput(l);
            }

            gate_recognition(candidates);

            root_clauses = index.estimateRoots();
        }

        std::set<Cl*> remainder;
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
        std::cerr << "c Starting gate-recognition with roots: " << roots << std::endl;
        std::vector<Lit> candidates;
        std::vector<Lit> frontier { roots.begin(), roots.end() };
        while (!frontier.empty()) { // _breadth_ first search is important here
            candidates.swap(frontier);
            // visit each candidate output only once per pass:
            candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
            for (Lit candidate : candidates) {
                if (isGate(candidate, patterns, semantic)) { 
                    unsigned middle = frontier.size();
                    frontier.insert(frontier.end(), gate_formula.getGate(candidate).inp.begin(), gate_formula.getGate(candidate).inp.end());
                    // assert that gate.inp is sorted
                    std::inplace_merge(frontier.begin(), frontier.begin() + middle, frontier.end()); 
                }
            }
            candidates.clear();
        }
    }

    /**
     * @brief Test if index contains a gate definition for candidate output literal 'out'
     * 
     * @param out candidate output literal
     * @param pat use gate patterns
     * @param sem use semantic recognition
     * @return true 
     * @return false 
     */
    bool isGate(Lit out, bool pat, bool sem) {
        if (index[~out].size() > 0 && index.isBlockedSet(out)) {
            const For& fwd = index[~out];
            const For& bwd = index[out];
            bool monotonic = gate_formula.isNestedMonotonic(out);

            if (monotonic || pat && fPattern(out, fwd, bwd) || sem && fSemantic(out, fwd, bwd)) {
                gate_formula.addGate(out, fwd, bwd); 
                index.remove(out.var());
                return true;
            }
        }
        return false;
    }

    // clause patterns of full encoding
    // precondition: fwd blocks bwd on output literal o
    bool fPattern(Lit o, const For& fwd, const For& bwd) {
        // check if fwd and bwd constrain exactly the same inputs (in opposite polarity)
        std::set<Var> fwd_inp, bwd_inp;
        for (Cl* c : fwd) for (Lit l : *c) if (l != ~o) fwd_inp.insert(l.var());
        for (Cl* c : bwd) for (Lit l : *c) if (l != o) bwd_inp.insert(l.var());
        if (fwd_inp != bwd_inp) {
            return false;
        }
        // detect equivalence gates
        if (fwd.size() == 1 && bwd.size() == 1 && fwd.front()->size() == 2 && bwd.front()->size() == 2) {
            return true;
        }
        // detect or gates
        if (fwd.size() == 1 && fixedClauseSize(bwd, 2)) {
            return true;
        }
        // detect and gates
        if (bwd.size() == 1 && fixedClauseSize(fwd, 2)) {
            return true;
        }
        // 2^n blocked clauses of size n+1 represent all input combinations 
        // each combined with one output literal
        if (fwd.size() == bwd.size() && 2*fwd.size() == pow(2, fwd_inp.size()/2)) {
            std::set<Var> vars;
            for (Lit l : fwd_inp) vars.insert(l.var());
            return 2*vars.size() == fwd_inp.size();
        }
        return false;
    }

    bool fSemantic(Lit o, const For& fwd, const For& bwd) {
        CNFFormula constraint;
        Cl clause;
        for (const For& f : { fwd, bwd }) {
            for (Cl* cl : f) {
                for (Lit l : *cl) {
                    if (l.var() != o.var()) clause.push_back(l);
                    else clause.push_back(Lit(o.var(), false));
                }
                constraint.readClause(clause.begin(), clause.end());
                clause.clear();
            }
        }    
        //constraint.normalizeVariableNames();
        for (Cl* clause : constraint) {
            for (Lit lit : *clause) {
                ipasir_add(S, lit.toDimacs());
            }
            ipasir_add(S, 0);
        }
        ipasir_assume(S, Lit(o.var(), true).toDimacs());
        int result = ipasir_solve(S);
        ipasir_add(S, Lit(o.var(), false).toDimacs());
        return result == 20;
    }

    bool fixedClauseSize(const For& f, unsigned int n) {
        for (Cl* c : f) if (c->size() != n) return false;
        return true;
    }

};

#endif
