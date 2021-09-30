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
#ifndef SRC_FEATURES_GATESTATS_H_
#define SRC_FEATURES_GATESTATS_H_

#include <math.h>

#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>

#include "src/util/SolverTypes.h"
#include "src/util/Runtime.h"
#include "src/gates/GateFormula.h"
#include "src/gates/GateAnalyzer.h"

class GateStats {
    CNFFormula formula_;
    Runtime runtime;

 public:
    unsigned n_vars, n_gates, n_roots;
    unsigned n_none, n_generic, n_mono, n_and, n_or, n_triv, n_equiv, n_full;

    std::vector<unsigned> levels, levels_none, levels_generic, levels_mono, levels_and, levels_or, levels_triv, levels_equiv, levels_full;

    explicit GateStats(const CNFFormula& formula) :
     formula_(formula), n_vars(formula.nVars()), n_gates(), n_roots(),
     n_none(0), n_generic(0), n_mono(0), n_and(0), n_or(0), n_triv(0), n_equiv(0), n_full(0),
     levels(), levels_none(), levels_generic(), levels_mono(), levels_and(), levels_or(), levels_triv(), levels_equiv(), levels_full() {
    }

    void analyze(unsigned repeat) {
        runtime.start();
        GateAnalyzer<> analyzer(formula_, true, true, repeat);
        analyzer.analyze();
        GateFormula gates = analyzer.getGateFormula();
        n_gates = gates.nGates();
        n_roots = gates.nRoots();
        gates.normalizeRoots();
        ++n_vars;
        levels.resize(n_vars, 0);
        // BFS for level determination
        unsigned level = 0;
        std::vector<Lit> current({ gates.getRoot() });
        std::vector<Lit> next;
        while (!current.empty()) {
            ++level;
            for (Lit lit : current) {
                const Gate& gate = gates.getGate(lit);
                if (gate.isDefined() && levels[lit.var()] == 0) {
                    levels[lit.var()] = level;
                    next.insert(next.end(), gate.inp.begin(), gate.inp.end());
                }
            }
            current.clear();
            current.swap(next);
        }
        // Gate Type Counts and Levels
        for (unsigned i = 1; i <= n_vars; i++) {
            const Gate& gate = gates.getGate(Lit(Var(i)));
            switch (gate.type) {
                case NONE:  // input variable
                    ++n_none;
                    levels_none.push_back(levels[i]);
                    break;
                case GENERIC:  // generically recognized gate
                    ++n_generic;
                    levels_generic.push_back(levels[i]);
                    break;
                case MONO:  // monotonically nested gate
                    ++n_mono;
                    levels_mono.push_back(levels[i]);
                    break;
                case AND:  // non-monotonically nested and-gate
                    ++n_and;
                    levels_and.push_back(levels[i]);
                    break;
                case OR:  // non-monotonically nested or-gate
                    ++n_or;
                    levels_or.push_back(levels[i]);
                    break;
                case TRIV:  // non-monotonically nested trivial equivalence gate
                    ++n_triv;
                    levels_triv.push_back(levels[i]);
                    break;
                case EQIV:  // non-monotonically nested equiv- or xor-gate
                    ++n_equiv;
                    levels_equiv.push_back(levels[i]);
                    break;
                case FULL:  // non-monotonically nested full gate (=maxterm encoding) with more than two inputs
                    ++n_full;
                    levels_full.push_back(levels[i]);
                    break;
            }
        }
        runtime.stop();
    }

    template <typename T>
    float Mean(std::vector<T> counts) {
        float sum = static_cast<float>(std::accumulate(counts.begin(), counts.end(), 0));
        return sum / counts.size();
    }

    template <typename T>
    float Variance(std::vector<T> counts, float mean) {
        float sum = static_cast<float>(std::accumulate(counts.begin(), counts.end(), 0.0,
            [mean] (float a, unsigned b) { return a + pow(static_cast<float>(b - mean), 2); } ));
        return sum / counts.size();
    }

    template <typename T>
    float Entropy(std::vector<T> counts) {
        float entropy = 0;
        for (unsigned count : counts) {
            float p_x = static_cast<float>(count) / counts.size();
            if (p_x > 0) entropy -= p_x * log(p_x) / log(2);
        }
        return entropy;
    }

    template <typename T>
    void push_distribution(std::vector<float>* record, std::vector<T> distribution) {
        float mean = 0, variance = 0, min = 0, max = 0, entropy = 0;
        if (distribution.size() > 0) {
            mean = Mean(distribution);
            variance = Variance(distribution, mean);
            min = static_cast<float>(*std::min_element(distribution.begin(), distribution.end()));
            max = static_cast<float>(*std::max_element(distribution.begin(), distribution.end()));
            entropy = Entropy(distribution);
        }
        record->push_back(mean);
        record->push_back(variance);
        record->push_back(min);
        record->push_back(max);
        record->push_back(entropy);
    }

    // Gate Structural Features
    std::vector<float> GateFeatures() {
        std::vector<float> record;
        record.push_back(n_vars);
        record.push_back(n_gates);
        record.push_back(n_roots);
        // gate or variable types
        record.push_back(n_none);
        record.push_back(n_generic);
        record.push_back(n_mono);
        record.push_back(n_and);
        record.push_back(n_or);
        record.push_back(n_triv);
        record.push_back(n_equiv);
        record.push_back(n_full);
        push_distribution(&record, levels);
        push_distribution(&record, levels_none);
        push_distribution(&record, levels_generic);
        push_distribution(&record, levels_mono);
        push_distribution(&record, levels_and);
        push_distribution(&record, levels_or);
        push_distribution(&record, levels_triv);
        push_distribution(&record, levels_equiv);
        push_distribution(&record, levels_full);
        record.push_back(static_cast<float>(runtime.get()));
        return record;
    }

    static std::vector<std::string> GateFeatureNames() {
        return std::vector<std::string> { "n_vars", "n_gates", "n_roots",
            "n_none", "n_generic", "n_mono", "n_and", "n_or", "n_triv", "n_equiv", "n_full",
            "levels_mean", "levels_variance", "levels_min", "levels_max", "levels_entropy",
            "levels_none_mean", "levels_none_variance", "levels_none_min", "levels_none_max", "levels_none_entropy",
            "levels_generic_mean", "levels_generic_variance", "levels_generic_min", "levels_generic_max", "levels_generic_entropy",
            "levels_mono_mean", "levels_mono_variance", "levels_mono_min", "levels_mono_max", "levels_mono_entropy",
            "levels_and_mean", "levels_and_variance", "levels_and_min", "levels_and_max", "levels_and_entropy",
            "levels_or_mean", "levels_or_variance", "levels_or_min", "levels_or_max", "levels_or_entropy",
            "levels_triv_mean", "levels_triv_variance", "levels_triv_min", "levels_triv_max", "levels_triv_entropy",
            "levels_equiv_mean", "levels_equiv_variance", "levels_equiv_min", "levels_equiv_max", "levels_equiv_entropy",
            "levels_full_mean", "levels_full_variance", "levels_full_min", "levels_full_max", "levels_full_entropy", "gate_extraction_time"
        };
    }
};

#endif  // SRC_FEATURES_GATESTATS_H_
