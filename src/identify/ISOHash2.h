/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Markus Iser, KIT - Karlsruhe Institute of Technology
ISOHash2 -- Copyright (c) 2025, Timon Passlick, KIT - Karlsruhe Institute of Technology (https://github.com/TimonPasslick/gbdc/blob/master/src/identify/ISOHash2.h)
ISOHash2 -- Copyright (c) 2025, Frederick Gehm, KIT - Karlsruhe Institute of Technology

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

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <optional>
#include <type_traits>
#include <unordered_set>
#include <vector>
#include <memory>
#include <fstream>
#include "src/external/md5/md5.h"
#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"

// debugging
#ifdef ISOHASH2_DEBUG
  #include <iostream>
  #define ISOH_LOG(MSG)  do { std::cerr << MSG << '\n'; } while (0)
#else
  #define ISOH_LOG(MSG)  do {} while (0)
#endif

namespace CNF {

struct WLHRuntimeConfig {
    unsigned depth;
    bool cross_reference_literals;
    bool rehash_clauses;
    bool optimize_first_iteration;
    unsigned progress_check_iteration;
    // bool shrink_to_fit;
    bool return_measurements;
    bool sort_for_clause_hash;
    bool use_xxh3;
    std::optional<unsigned> prime_ring_modulus;
};

class WeisfeilerLemanHasher {
public:
    using Clock = std::chrono::high_resolution_clock;
    using Clause = ::Cl;
    using ClausePtr = const ::Cl*;
    using Lit = ::Lit;
    using Hash = std::uint64_t;

    struct LitColors {
        Hash p = 0;
        Hash n = 0;
        void flip() { std::swap(n, p); }
        Hash variable_hash(const WeisfeilerLemanHasher* parent) const {
            LitColors copy = *this;
            if (n > p) copy.flip();
            return parent->hash(copy);
        }
        void cross_reference(WeisfeilerLemanHasher* parent) {
            const Hash pcr = parent->hash(*this);
            flip();
            const Hash ncr = parent->hash(*this);
            p = pcr;
            n = ncr;
        }
    };

    struct ColorFunction {
        std::vector<LitColors> colors;
        explicit ColorFunction(std::size_t n) : colors(n+1, {1, 1}) {}
        Hash operator()(const ::Lit& lit) const { return const_cast<ColorFunction*>(this)->operator()(lit); }
        Hash& operator()(const ::Lit& lit) { 
            unsigned var_id = lit.var().id;
            if (var_id >= colors.size()) { // Check if var_id is too high
                throw std::out_of_range("Variable ID " + std::to_string(var_id) + " is out of range for color access.");
            }
            if (lit.sign()) { // true if negative
                return colors[var_id].n;
            } else { // false if positive
                return colors[var_id].p;
            }
            // return reinterpret_cast<Hash*>(&colors[0])[lit]; 
        }
    };

    WeisfeilerLemanHasher(const char* filename, const WLHRuntimeConfig& cfg)
        : cfg(cfg),
          parsing_start_mem(get_mem_usage()),
          parsing_start_time(Clock::now()),
          cnf(filename),
          start_mem(get_mem_usage()),
          start_time(Clock::now()),
          color_functions{ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())}
    {
        if (cfg.prime_ring_modulus && *cfg.prime_ring_modulus < 2) {
            throw std::invalid_argument("prime_ring_modulus must be >= 2 if specified.");
    }
    ISOH_LOG("Loaded CNF  vars=" << cnf.nVars()
                << " clauses=" << cnf.nClauses()
                << " lits=" << cnf.nLits()
                << " maxClauseLen=" << cnf.maxClauseLength());
    }

    std::string operator()() {
        std::string result = std::to_string(run());
        if (cfg.return_measurements) {
            const auto calculation_time = std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start_time).count();
            const auto parsing_time = std::chrono::duration_cast<std::chrono::nanoseconds>(start_time - parsing_start_time).count();
            long parsing_mem_usage = start_mem;
            long mem_usage = get_mem_usage();
            if (parsing_start_mem == -1 || parsing_mem_usage == -1 || mem_usage == -1) {
                mem_usage = -1;
                parsing_mem_usage = -1;
            } else {
                mem_usage -= start_mem;
                parsing_mem_usage -= parsing_start_mem;
            }
            const double iteration_count = std::min<double>(iteration, cfg.depth / 2.);
            result +=
                "," + std::to_string(parsing_time) +
                "," + std::to_string(calculation_time) +
                "," + std::to_string(parsing_mem_usage) +
                "," + std::to_string(mem_usage) +
                "," + std::to_string(iteration_count) +
                "," + std::to_string(cnf.nVars()) +
                "," + std::to_string(cnf.nClauses()) +
                "," + std::to_string(cnf.nLits()) +
                "," + std::to_string(cnf.maxClauseLength());
        }
        return result;
    }

private:
    const WLHRuntimeConfig cfg;
    long parsing_start_mem;
    Clock::time_point parsing_start_time;
    CNFFormula cnf;
    long start_mem;
    Clock::time_point start_time;
    ColorFunction color_functions[2];
    unsigned iteration = 0;
    std::unordered_set<Hash> unique_hashes;
    unsigned previous_unique_hashes = 1;

    ColorFunction& old_color() { return color_functions[iteration % 2]; }
    ColorFunction& new_color() { return color_functions[(iteration + 1) % 2]; }

    template <typename T>
    Hash hash(const T& t) const {
        if (!cfg.prime_ring_modulus.has_value()) {
            if (cfg.use_xxh3)
                return XXH3_64bits(&t, sizeof(t));
            MD5 md5;
            md5.consume_binary(t);
            return md5.finish();
        }
        constexpr std::uint64_t max = std::numeric_limits<std::uint64_t>::max();
        const std::uint64_t mod = cfg.prime_ring_modulus.value();
        std::uint64_t hash = max;
        const std::uint64_t first_problem = max - (max % mod);
        for (std::uint16_t seed = 0; hash >= first_problem; ++seed) {
            if (cfg.use_xxh3)
                hash = XXH3_64bits_withSecret(&t, sizeof(t), &seed, sizeof(seed));
            else {
                MD5 md5;
                md5.consume_binary(seed);
                md5.consume_binary(t);
                hash = md5.finish();
            }
        }
        return hash % mod;
    }

    static void combine(Hash* acc, Hash in, std::optional<unsigned> mod) {
        if (mod.has_value()) {
            *acc= (*acc + (in % *mod)) % *mod;
        } else {
            *acc += in;
        }
    }

    template <typename T, typename C>
    Hash hash_sum(const C& c, const std::function<Hash(const T&)>& f) const {
        Hash h = 0;
        for (const T& t : c)
            combine(&h, f(t), cfg.prime_ring_modulus);
        return h;
    }

    bool in_optimized_iteration() const {
        return iteration == 0 && cfg.optimize_first_iteration;
    }

    Hash variable_hash() {
        if (cfg.cross_reference_literals) {
            return hash_sum<LitColors>(old_color().colors, 
                [this](const LitColors& lc) { return lc.variable_hash(this); });
        }
        Hash h = 0;
        for (unsigned var = 1; var <= cnf.nVars(); ++var) {
            combine(&h, old_color()(Lit(var, false)), cfg.prime_ring_modulus); // p
            combine(&h, old_color()(Lit(var, true )), cfg.prime_ring_modulus); // n
        }
        return h;
    }

    // Main iteration step
    void iteration_step() {
        ISOH_LOG("Iteration " << iteration);

        std::fill(new_color().colors.begin(), new_color().colors.end(), LitColors{0,0});

        cross_reference();

        for (ClausePtr clp : cnf.clauses()) {
            const Clause& cl = *clp;
            const Hash clh = (!in_optimized_iteration()) 
                ? clause_hash(cl) 
                : cfg.rehash_clauses ? hash(cl.size()) : cl.size();

            for (const auto& lit : cl)
                combine(&new_color()(lit), clh, cfg.prime_ring_modulus);
        }
        ++iteration;
    }

    void cross_reference() {
        if (!cfg.cross_reference_literals || in_optimized_iteration())
            return;
        for (LitColors& lc : old_color().colors)
            lc.cross_reference(this);
    }

    Hash clause_hash(const Clause& cl) {
        if (!cfg.sort_for_clause_hash) {
            Hash h = hash_sum<typename Clause::value_type>(cl, [this](const auto& lit) { return old_color()(lit); });
            if (cfg.rehash_clauses) h = hash(h);
            return h;
        } else {
            std::vector<Hash> sorted;
            sorted.reserve(cl.size());
            for (const auto& lit : cl)
                sorted.push_back(old_color()(lit));
            std::sort(sorted.begin(), sorted.end());
            return XXH3_64bits(sorted.data(), sorted.size() * sizeof(Hash));
        }
    }

    std::optional<Hash> check_progress() {
        if ((iteration != cfg.progress_check_iteration && iteration != cfg.progress_check_iteration + 1 && iteration < 6) || iteration == 0)
            return std::nullopt;

        unique_hashes.reserve(previous_unique_hashes);
        const Hash vh = hash_sum<LitColors>(old_color().colors, [this](const LitColors& lc) {
            const Hash vh = lc.variable_hash(this);
            unique_hashes.insert(vh);
            return vh;
        });
        if (unique_hashes.size() <= previous_unique_hashes) {
            ISOH_LOG("Stabilised at iter=" << iteration
                    << " uniqueVars=" << unique_hashes.size()
                    << "  early hash=" << vh);
            return vh;
        }
        previous_unique_hashes = unique_hashes.size();
        unique_hashes.clear();
        return std::nullopt;
    }

    Hash cnf_hash() {
        cross_reference();
        Hash h = 0;
        for (ClausePtr clp : cnf.clauses()) {
            combine(&h, clause_hash(*clp), cfg.prime_ring_modulus);
        }
        return h;
    }

    Hash run() {
        ISOH_LOG("Start hashing, maxDepth=" << cfg.depth);
        while (iteration < cfg.depth / 2) {
            if (const auto result = check_progress())
                return *result;
            iteration_step();
        }
        Hash final = cfg.depth % 2 == 0 ? variable_hash() : cnf_hash();
        ISOH_LOG("Finished at iter=" << iteration << "  finalHash=" << final);
        return final;
    }
};

inline std::string weisfeiler_leman_hash(
    const char* filename,
    unsigned depth = 100,
    bool cross_reference_literals = true,
    bool rehash_clauses = true,
    bool optimize_first_iteration = true,
    unsigned progress_check_iteration = 6,
    // bool shrink_to_fit = false,
    bool return_measurements = false,
    bool sort_for_clause_hash = false,
    bool use_xxh3 = true,
    std::optional<unsigned> prime_ring_modulus = std::nullopt
    // bool use_half_word_hash = true,
    // bool use_prime_ring = false 
    )
    {
    WLHRuntimeConfig cfg{
        depth,
        cross_reference_literals,
        rehash_clauses,
        optimize_first_iteration,
        progress_check_iteration,
        // shrink_to_fit,
        return_measurements,
        sort_for_clause_hash,
        use_xxh3,
        prime_ring_modulus // instead of use_half_word_hash and use_prime_ring
        // use_half_word_hash,
        // use_prime_ring
    };
    WeisfeilerLemanHasher hasher(filename, cfg);
    return hasher();
}

} // namespace CNF

#endif  // ISOHASH2_H_