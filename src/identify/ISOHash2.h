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
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <utility>

#ifndef XXH_INLINE_ALL
#define XXH_INLINE_ALL
#endif
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"

namespace CNF {

struct IsoHash2Settings {
    std::uint64_t max_iterations = 6;
    bool print_stats = false;
};

class IsoHash2 {
public:
    using Hash = std::uint64_t;
    using Clause = Cl; 
    using Literal = Lit;

    struct Stats {
        Hash hash = 0;
        std::uint64_t round = 0;
        bool stabilized = false;
    };

private:
    struct LitColors {
        Hash p = 1;
        Hash n = 1;
    };

    struct Fingerprint {
        Hash sum_hash = 0;
        Hash xor_hash = 0;
        
        bool operator==(const Fingerprint& other) const {
            return sum_hash == other.sum_hash && xor_hash == other.xor_hash;
        }
    };

    struct Signature {
        Hash p;
        Hash n;
        std::uint32_t original_index;

        bool operator<(const Signature& other) const {
            if (p != other.p) return p < other.p;
            return n < other.n;
        }
        bool operator==(const Signature& other) const {
            return p == other.p && n == other.n;
        }
    };

    IsoHash2Settings settings;
    const CNFFormula& cnf;

    std::vector<LitColors> colors_[2];
    std::vector<Signature> sort_buffer;

    inline Hash& color_of(std::vector<LitColors>& vec, Literal lit) {
        return lit.sign() ? vec[lit.var()].n : vec[lit.var()].p;
    }

    inline Hash color_of(const std::vector<LitColors>& vec, Literal lit) const {
        return lit.sign() ? vec[lit.var()].n : vec[lit.var()].p;
    }

    static inline Hash hash_bytes(const void* data, std::size_t len) {
        return XXH3_64bits(data, len);
    }

    // mix64variant13 [Steele et al. 2014]
    static inline Hash fast_mix(Hash k) {
        k ^= k >> 30; k *= 0xbf58476d1ce4e5b9ULL;
        k ^= k >> 27; k *= 0x94d049bb133111ebULL;
        k ^= k >> 31;
        return k;
    }

    Hash clause_hash(const Clause& clause, const std::vector<LitColors>& cur_vec) const {
        Hash combined = 0;
        for (const Literal lit : clause) {
            combined += color_of(cur_vec, lit);
        }
        return fast_mix(combined);
    }

    Fingerprint finalize(const std::vector<LitColors>& cur_vec, std::vector<LitColors>& next_vec) {
        const std::size_t num_vars = cnf.nVars();

        for (std::size_t var_id = 1; var_id <= num_vars; ++var_id) {
            auto& agg_lc = next_vec[var_id];
            const auto& cur_lc = cur_vec[var_id];

            const Hash features_p[3] = {cur_lc.p, agg_lc.p, cur_lc.n};
            Hash h_p = hash_bytes(features_p, sizeof(features_p));

            const Hash features_n[3] = {cur_lc.n, agg_lc.n, cur_lc.p};
            Hash h_n = hash_bytes(features_n, sizeof(features_n));

            sort_buffer[var_id-1] = {h_p, h_n, static_cast<std::uint32_t>(var_id)};
        }

        std::sort(sort_buffer.begin(), sort_buffer.end());

        Hash acc_sum = 0;
        Hash acc_xor = 0;
        Hash current_rank = 0;

        for (std::size_t i = 0; i < num_vars; ++i) {
            if (i > 0 && !(sort_buffer[i] == sort_buffer[i-1])) {
                current_rank++;
            }

            Hash stable_color = fast_mix(current_rank);
            std::uint32_t original_idx = sort_buffer[i].original_index;

            next_vec[original_idx].p = stable_color;
            next_vec[original_idx].n = stable_color;

            acc_sum += stable_color;
            acc_xor ^= stable_color;
        }

        return {acc_sum, acc_xor};
    }

    Fingerprint iteration_step(int cur_idx, int next_idx) {
        auto& cur_vec = colors_[cur_idx];
        auto& next_vec = colors_[next_idx];

        std::fill(next_vec.begin(), next_vec.end(), LitColors{0, 0});

        for (const auto& clause_ptr : cnf) {
            const Clause& clause = *clause_ptr;
            const Hash ch = clause_hash(clause, cur_vec);
            for (const Literal lit : clause) {
                color_of(next_vec, lit) += ch;
            }
        }
        return finalize(cur_vec, next_vec);
    }

public:
    IsoHash2(const CNFFormula& formula, const IsoHash2Settings& s) :
        settings(s),
        cnf(formula),
        colors_{std::vector<LitColors>(cnf.nVars() + 1), std::vector<LitColors>(cnf.nVars() + 1)},
        sort_buffer(cnf.nVars())
    {}

    Stats run() {
        Stats stats{};
        Fingerprint prev_fingerprint = {0,0};

        int cur_idx = 0;
        int next_idx = 1;

        while (stats.round < settings.max_iterations) {
            Fingerprint cur_fingerprint = iteration_step(cur_idx, next_idx);
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
            std::swap(cur_idx, next_idx);
        }

        stats.hash = hash_bytes(&prev_fingerprint, sizeof(prev_fingerprint));

        if (!stats.stabilized && settings.print_stats) {
            std::cerr << "c Reached max iterations (" << settings.max_iterations << ")." << std::endl;
        }

        return stats;
    }
};

inline IsoHash2::Stats isohash2_stats(const char* filename, const IsoHash2Settings& s = {}) {
    CNFFormula cnf(filename);
    IsoHash2 hasher(cnf, s);
    return hasher.run();
}

inline std::string isohash2(const char* filename, const IsoHash2Settings& s = {}) {
    const auto stats = isohash2_stats(filename, s);

    std::ostringstream oss;
    oss << std::hex << std::setw(16) << std::setfill('0') << std::nouppercase << stats.hash;

    return oss.str();
}

} // namespace CNF

#endif // ISOHASH2_H_
