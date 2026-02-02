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

#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <sstream>

#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"

namespace CNF {

struct IsoHash2Settings {
    int max_iterations = 20; // 0 = until stabilized
    bool print_stats = false;
};

class IsoHash2 {
public:
    using Hash = uint64_t;
    using Clause = Cl;
    using Literal = Lit;

    struct Stats {
        Hash hash = 0;
        int round = 0;
        bool stabilized = false;
    };

private:
    struct LitColors {
        Hash val[2] = {1, 1};
    };

    struct ColorFunction {
        std::vector<LitColors> colors_by_var;
        explicit ColorFunction(size_t n_vars) : colors_by_var(n_vars + 1) {}

        inline Hash& operator()(Lit lit) {
            return colors_by_var[lit.var()].val[(int)lit.sign()];
        }
        inline const Hash& operator()(Lit lit) const {
            return colors_by_var[lit.var()].val[(int)lit.sign()];
        }
    };

    // mix64variant13 [Steele et al. 2014]
    static inline Hash fast_mix(Hash k) {
        k ^= k >> 30; k *= 0xbf58476d1ce4e5b9ULL;
        k ^= k >> 27; k *= 0x94d049bb133111ebULL;
        k ^= k >> 31;
        return k;
    }

    static inline Hash rotl64(Hash x, int r) {
        return (x << r) | (x >> (64 - r));
    }

    const IsoHash2Settings& settings;
    const CNFFormula& cnf;
    ColorFunction color_functions[2];
    Stats stats;

    std::vector<Hash> partition_buffer;
    size_t prev_partition_count = 0;

    inline ColorFunction& old_color() { return color_functions[stats.round % 2]; }
    inline const ColorFunction& old_color() const { return color_functions[stats.round % 2]; }

    inline ColorFunction& new_color() { return color_functions[(stats.round + 1) % 2]; }
    inline const ColorFunction& new_color() const { return color_functions[(stats.round + 1) % 2]; }

    inline Hash state_hash_oriented(const LitColors& lc) const {
        Hash p = lc.val[0];
        Hash n = lc.val[1];
        Hash x = p ^ rotl64(n, 1);
        return fast_mix(x + 0x9e3779b97f4a7c15ULL);
    }

    inline Hash state_hash_canonical(const LitColors& lc) const {
        Hash p = lc.val[0];
        Hash n = lc.val[1];
        if (p > n) std::swap(p, n);
        Hash x = p ^ rotl64(n, 1);
        return fast_mix(x + 0x9e3779b97f4a7c15ULL);
    }

    Hash clause_hash(const Clause& clause) const {
        Hash a = 0, b = 0;
        const auto& c_func = old_color();

        for (const Literal lit : clause) {
            Hash y = fast_mix(c_func(lit) + 0x9e3779b97f4a7c15ULL);
            a += y;
            b ^= rotl64(y, 23);
        }

        Hash combined = a ^ fast_mix(b + 0xbf58476d1ce4e5b9ULL) ^ (Hash)clause.size();
        return fast_mix(combined);
    }

    void finalize_literal_colors() {
        const size_t num_vars = cnf.nVars();
        auto* agg_vec = &new_color().colors_by_var;
        const auto* old_vec = &old_color().colors_by_var;

        for (size_t i = 1; i <= num_vars; ++i) {
            Hash old_p = (*old_vec)[i].val[0];
            Hash old_n = (*old_vec)[i].val[1];
            Hash agg_p = (*agg_vec)[i].val[0];
            Hash agg_n = (*agg_vec)[i].val[1];

            Hash x_p = old_p + fast_mix(agg_p) + rotl64(old_n, 1);
            Hash x_n = old_n + fast_mix(agg_n) + rotl64(old_p, 1);

            (*agg_vec)[i].val[0] = fast_mix(x_p);
            (*agg_vec)[i].val[1] = fast_mix(x_n);
        }
    }

    void iteration_step() {
        auto& nc_vec = new_color().colors_by_var;
        std::memset(nc_vec.data(), 0, nc_vec.size() * sizeof(LitColors));

        for (const auto& clause_ptr : cnf) {
            Hash ch = clause_hash(*clause_ptr);
            for (const Literal lit : *clause_ptr) {
                new_color()(lit) += ch;
            }
        }
        finalize_literal_colors();
    }

    bool check_stabilization() {
        const size_t n = cnf.nVars();
        if (partition_buffer.size() != n) partition_buffer.resize(n);

        const auto& current_colors = old_color().colors_by_var;
        for (size_t i = 1; i <= n; ++i) {
            partition_buffer[i - 1] = state_hash_oriented(current_colors[i]);
        }

        std::sort(partition_buffer.begin(), partition_buffer.end());

        size_t current_partition_count = 0;
        if (n > 0) {
            current_partition_count = 1;
            for (size_t i = 1; i < n; ++i) {
                current_partition_count += (partition_buffer[i] != partition_buffer[i - 1]);
            }
        }

        bool stable = (current_partition_count == prev_partition_count);
        prev_partition_count = current_partition_count;
        return stable;
    }

public:
    IsoHash2(const CNFFormula& formula, const IsoHash2Settings& s) :
        settings(s),
        cnf(formula),
        color_functions{ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())},
        partition_buffer(cnf.nVars())
    {}

    Stats run() {
        stats = Stats{};
        prev_partition_count = 0;

        while (stats.round < settings.max_iterations || settings.max_iterations == 0) {
            iteration_step();
            stats.round++;

            bool stable = check_stabilization();

            if (settings.print_stats) std::cerr << "c Round " << stats.round << " partitions: " << prev_partition_count << "\n";

            if (stable) {
                stats.stabilized = true;
                if (settings.print_stats) std::cerr << "c Stabilized after " << stats.round << " rounds.\n";
                break;
            }
        }

        if (!stats.stabilized && settings.print_stats) std::cerr << "c Reached max iterations (" << settings.max_iterations << ").\n";

        // FINAL HASH
        const size_t n = cnf.nVars();
        if (partition_buffer.size() != n) partition_buffer.resize(n);

        const auto& final_colors = old_color().colors_by_var;
        for (size_t i = 1; i <= n; ++i) {
            partition_buffer[i - 1] = state_hash_canonical(final_colors[i]);
        }
        std::sort(partition_buffer.begin(), partition_buffer.end());
        stats.hash = XXH3_64bits(partition_buffer.data(), partition_buffer.size() * sizeof(Hash));
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
