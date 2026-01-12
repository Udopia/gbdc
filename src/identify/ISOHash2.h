/*************************************************************************************************
CNFTools -- Copyright (c) 2026, Ashlin Iser, Frederick Gehm, KIT - Karlsruhe Institute of Technology

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

#ifndef ISOHASH2_H_
#define ISOHASH2_H_

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

// #ifndef XXH_INLINE_ALL
#define XXH_INLINE_ALL
// #endif
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"

namespace CNF {

struct Settings {
    uint64_t max_iterations = 6;
    bool print_stats = false;
};

class IsoHash2 {
public:
    using Hash = std::uint64_t;
    using Clause = Cl; 
    using Literal = Lit;

    struct Stats {
        Hash hash = 0;
        uint64_t round = 0;
        bool stabilized = false;
    };

    struct Fingerprint {
        Hash sum_hash;
        Hash xor_hash;
        
        bool operator==(const Fingerprint& other) const {
            return sum_hash == other.sum_hash && xor_hash == other.xor_hash;
        }
    };

private:
    struct LitColors {
        Hash p = 1;
        Hash n = 1;
    };

    struct ColorFunction {
        std::vector<LitColors> colors_by_var;

        explicit ColorFunction(size_t n_vars) : colors_by_var(n_vars + 1) {}

        inline Hash& operator()(Literal lit) {
            return lit.sign() ? colors_by_var[lit.var()].n : colors_by_var[lit.var()].p;
        }
        inline const Hash& operator()(Literal lit) const {
            return lit.sign() ? colors_by_var[lit.var()].n : colors_by_var[lit.var()].p;
        }
    };

    const Settings& settings;
    const CNFFormula& cnf;

    ColorFunction color_functions[2];

    inline ColorFunction& current_color(std::uint64_t round) { return color_functions[round % 2]; }
    inline const ColorFunction& current_color(std::uint64_t round) const { return color_functions[round % 2]; }

    inline ColorFunction& next_color(std::uint64_t round) { return color_functions[(round + 1) % 2]; }
    inline const ColorFunction& next_color(std::uint64_t round) const { return color_functions[(round + 1) % 2]; }

    static inline Hash hash_bytes(const void* data, std::size_t len) {
        return XXH3_64bits(data, len);
    }

    Hash variable_hash(const LitColors& lc) const {
        Hash data[2];
        if (lc.p > lc.n) { data[0] = lc.n; data[1] = lc.p; }
        else { data[0] = lc.p; data[1] = lc.n; }
        return hash_bytes(data, sizeof(data));
    }

    Hash clause_hash(const Clause& clause, std::uint64_t round) const {
        Hash combined = 0;
        for (const Literal lit : clause) {
            combined+=current_color(round)(lit);
        }
        return hash_bytes(&combined, sizeof(combined));
    }

    Fingerprint finalize(std::uint64_t round) {
        auto& next_vec = next_color(round).colors_by_var;
        const auto& cur_vec = current_color(round).colors_by_var;

        Hash acc_sum = 0;
        Hash acc_xor = 0;

        for (size_t var_id = 1; var_id <= cnf.nVars(); ++var_id) {
            auto& agg_lc = next_vec[var_id];
            const auto& cur_lc = cur_vec[var_id];

            const Hash features_p[3] = {cur_lc.p, agg_lc.p, cur_lc.n};
            agg_lc.p = hash_bytes(features_p, sizeof(features_p));

            const Hash features_n[3] = {cur_lc.n, agg_lc.n, cur_lc.p};
            agg_lc.n = hash_bytes(features_n, sizeof(features_n));

            const Hash vh = variable_hash(agg_lc);
            acc_sum += vh;
            acc_xor ^= vh;
        }

        return {acc_sum, acc_xor};
    }


    Fingerprint iteration_step(std::uint64_t round) {
        auto& next_vec = next_color(round).colors_by_var;
        std::fill(next_vec.begin(), next_vec.end(), LitColors{0, 0});

        for (const auto& clause_ptr : cnf) {
            const Clause& clause = *clause_ptr;
            const Hash ch = clause_hash(clause, round);
            for (const Literal lit : clause) {
                next_color(round)(lit) += ch;
            }
        }
        return finalize(round);
    }

public:
    IsoHash2(const CNFFormula& formula, const Settings& s) :
        settings(s),
        cnf(formula),
        color_functions{ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())}
    {}

    Stats run() {
        Stats stats{};
        Fingerprint prev_fingerprint = {0,0};

        while (stats.round < settings.max_iterations) {
            Fingerprint cur_fingerprint = iteration_step(stats.round);
            stats.round++;

            if (stats.round > 1 && cur_fingerprint == prev_fingerprint) {
                stats.stabilized = true;
                stats.hash = hash_bytes(&cur_fingerprint, sizeof(cur_fingerprint));
                if (settings.print_stats) {
                    std::cerr << "c Stabilized after " << stats.round << " rounds." << std::endl;
                }
                return stats;
            }

            if (settings.print_stats) {
                std::cerr << "c Round " << stats.round << " Hash: " << cur_fingerprint.sum_hash << std::endl;
            }

            prev_fingerprint = cur_fingerprint;
        }

        stats.hash = hash_bytes(&prev_fingerprint, sizeof(prev_fingerprint));

        if (!stats.stabilized && settings.print_stats) {
            std::cerr << "c Reached max iterations (" << settings.max_iterations << ")." << std::endl;
        }

        return stats;
    }
};

inline IsoHash2::Stats isohash2_stats(const char* filename, const Settings& s = {}) {
    CNFFormula cnf(filename);
    IsoHash2 hasher(cnf, s);
    return hasher.run();
}

inline std::string isohash2(const char* filename) {
    const Settings& s = {};
    const auto stats = isohash2_stats(filename, s);

    std::ostringstream oss;
    oss << std::hex << std::setw(16) << std::setfill('0') << std::nouppercase << stats.hash;

    return oss.str();
}


} // namespace CNF

#endif // ISOHASH2_H_