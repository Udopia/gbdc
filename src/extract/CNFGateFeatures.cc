/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#include "src/extract/CNFGateFeatures.h"

#include "src/util/SolverTypes.h"
#include "src/extract/gates/GateFormula.h"
#include "src/extract/gates/GateAnalyzer.h"
#include "src/util/CaptureDistribution.h"

CNF::GateFeatures::GateFeatures(const char* filename) : filename_(filename) { 
    initFeatures({ "n_vars", "n_gates", "n_roots" });
    initFeatures({ "n_none", "n_generic", "n_mono" });
    initFeatures({ "n_and", "n_or", "n_triv", "n_equiv", "n_full" });
    initFeatures({ "levels_mean", "levels_variance", "levels_min", "levels_max", "levels_entropy" });
    initFeatures({ "levels_none_mean", "levels_none_variance", "levels_none_min", "levels_none_max", "levels_none_entropy" });
    initFeatures({ "levels_generic_mean", "levels_generic_variance", "levels_generic_min", "levels_generic_max", "levels_generic_entropy" });
    initFeatures({ "levels_mono_mean", "levels_mono_variance", "levels_mono_min", "levels_mono_max", "levels_mono_entropy" });
    initFeatures({ "levels_and_mean", "levels_and_variance", "levels_and_min", "levels_and_max", "levels_and_entropy" });
    initFeatures({ "levels_or_mean", "levels_or_variance", "levels_or_min", "levels_or_max", "levels_or_entropy" });
    initFeatures({ "levels_triv_mean", "levels_triv_variance", "levels_triv_min", "levels_triv_max", "levels_triv_entropy" });
    initFeatures({ "levels_equiv_mean", "levels_equiv_variance", "levels_equiv_min", "levels_equiv_max", "levels_equiv_entropy" });
    initFeatures({ "levels_full_mean", "levels_full_variance", "levels_full_min", "levels_full_max", "levels_full_entropy" });
}

CNF::GateFeatures::~GateFeatures() { }

void CNF::GateFeatures::run() {
    CNFFormula formula(filename_);
    GateAnalyzer analyzer(formula, true, true, formula.nVars() / 3, false);
    analyzer.analyze();
    GateFormula gates = analyzer.getGateFormula();
    n_vars = formula.nVars();
    n_gates = gates.nGates();
    n_roots = gates.nRoots();
    levels.resize(n_vars + 1, 0);
    // BFS for level determination
    unsigned level = 0;
    std::vector<Lit> current = gates.getRoots();
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
    load_feature_records();
}

void CNF::GateFeatures::load_feature_records() {
    setFeature("n_vars", (double)n_vars);
    setFeature("n_gates", (double)n_gates);
    setFeature("n_roots", (double)n_roots);
    setFeature("n_none", (double)n_none);
    setFeature("n_generic", (double)n_generic);
    setFeature("n_mono", (double)n_mono);
    setFeature("n_and", (double)n_and);
    setFeature("n_or", (double)n_or);
    setFeature("n_triv", (double)n_triv);
    setFeature("n_equiv", (double)n_equiv);
    setFeature("n_full", (double)n_full);
    std::vector<double> stats = getDistributionStats(levels);
    setFeatures({ "levels_mean", "levels_variance", "levels_min", "levels_max", "levels_entropy" }, stats.begin(), stats.end());
    std::vector<double> stats_none = getDistributionStats(levels_none);
    setFeatures({ "levels_none_mean", "levels_none_variance", "levels_none_min", "levels_none_max", "levels_none_entropy" }, stats_none.begin(), stats_none.end());
    std::vector<double> stats_generic = getDistributionStats(levels_generic);
    setFeatures({ "levels_generic_mean", "levels_generic_variance", "levels_generic_min", "levels_generic_max", "levels_generic_entropy" }, stats_generic.begin(), stats_generic.end());
    std::vector<double> stats_mono = getDistributionStats(levels_mono);
    setFeatures({ "levels_mono_mean", "levels_mono_variance", "levels_mono_min", "levels_mono_max", "levels_mono_entropy" }, stats_mono.begin(), stats_mono.end());
    std::vector<double> stats_and = getDistributionStats(levels_and);
    setFeatures({ "levels_and_mean", "levels_and_variance", "levels_and_min", "levels_and_max", "levels_and_entropy" }, stats_and.begin(), stats_and.end());
    std::vector<double> stats_or = getDistributionStats(levels_or);
    setFeatures({ "levels_or_mean", "levels_or_variance", "levels_or_min", "levels_or_max", "levels_or_entropy" }, stats_or.begin(), stats_or.end());
    std::vector<double> stats_triv = getDistributionStats(levels_triv);
    setFeatures({ "levels_triv_mean", "levels_triv_variance", "levels_triv_min", "levels_triv_max", "levels_triv_entropy" }, stats_triv.begin(), stats_triv.end());
    std::vector<double> stats_equiv = getDistributionStats(levels_equiv);
    setFeatures({ "levels_equiv_mean", "levels_equiv_variance", "levels_equiv_min", "levels_equiv_max", "levels_equiv_entropy" }, stats_equiv.begin(), stats_equiv.end());
    std::vector<double> stats_full = getDistributionStats(levels_full);
    setFeatures({ "levels_full_mean", "levels_full_variance", "levels_full_min", "levels_full_max", "levels_full_entropy" }, stats_full.begin(), stats_full.end());
}