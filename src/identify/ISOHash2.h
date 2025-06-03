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

namespace CNF {

struct WLHRuntimeConfig {
    unsigned depth;
    bool cross_reference_literals;
    bool rehash_clauses;
    bool optimize_first_iteration;
    unsigned progress_check_iteration;
    // bool shrink_to_fit;
    bool no_measurements;
    bool sort_for_clause_hash;
    bool use_xxh3;
    std::optional<unsigned> prime_ring_modulus;
};

struct WLHResult {
    std::string filename;
    uint64_t hash;
    std::string hexhash;
    unsigned iterations;
    std::string status;
    long total_runtime;
    std::vector<long> iterations_time;
    std::vector<unsigned> unique_variables;

    unsigned nVars, nClauses, nLits, maxClauseLen;

    WLHRuntimeConfig cfg;
};


class WeisfeilerLemanHasher {
public:
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
            if (var_id >= colors.size()) {
                throw std::out_of_range("Variable ID " + std::to_string(var_id) + " is out of range for color access.");
            }
            if (lit.sign()) { // true if n
                return colors[var_id].n;
            } else { // false if p
                return colors[var_id].p;
            }
        }
    };

    WeisfeilerLemanHasher(const char* filename, const WLHRuntimeConfig& cfg)
        : cfg(cfg),
          cnf(filename),
          color_functions{ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())}
    {
        if (cfg.prime_ring_modulus && *cfg.prime_ring_modulus < 2) {
            throw std::invalid_argument("prime_ring_modulus must be >= 2 if specified.");
    }
    }

    WLHResult collect_measurements(const std::string& original_filename) {
        iter_us.clear();
        uniq_cnt.clear();
        unique_hashes.clear();
        iteration = 0;
        previous_unique_hashes = 1;
        status = "running";
        uint64_t final_hash = run();
        std::ostringstream hexhash;
        hexhash << std::hex << std::setw(16) << std::setfill('0') << final_hash;

        WLHResult r{};

        r.filename = original_filename;
        r.hash = final_hash;
        r.hexhash = hexhash.str();
        r.iterations = iteration;
        r.status = status;
        r.iterations_time = iter_us;
        r.unique_variables = uniq_cnt;
        r.nVars = cnf.nVars();
        r.nClauses = cnf.nClauses();
        r.nLits = cnf.nLits();
        r.maxClauseLen = cnf.maxClauseLength();
        r.cfg = cfg;

        return r;
    }

    std::string operator()() {
        return std::to_string(run());
    }

private:
    const WLHRuntimeConfig cfg;
    CNFFormula cnf;
    ColorFunction color_functions[2];
    unsigned iteration = 0;
    std::unordered_set<Hash> unique_hashes;
    unsigned previous_unique_hashes = 1;

    std::vector<long> iter_us;
    std::vector<unsigned> uniq_cnt;
    std::string status = "running";

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
        // measurements
        auto t0 = !cfg.no_measurements
              ? std::chrono::high_resolution_clock::now()
              : std::chrono::high_resolution_clock::time_point{};
        // measurements

        // reset new color function
        std::fill(new_color().colors.begin(), new_color().colors.end(), LitColors{0,0});

        // cross reference positive and negative literals
        cross_reference();

        // clause processing
        for (ClausePtr clp : cnf.clauses()) {
            const Clause& cl = *clp;
            const Hash clh = (!in_optimized_iteration()) 
                ? clause_hash(cl) 
                : cfg.rehash_clauses ? hash(cl.size()) : cl.size();

            // accumulate new lit color hash
            for (const auto& lit : cl)
                combine(&new_color()(lit), clh, cfg.prime_ring_modulus);
        }

        // measurements
        if (!cfg.no_measurements) {
            auto t1 = std::chrono::high_resolution_clock::now();
            long dt = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
            iter_us.push_back(dt);
        }
        // measurements

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
            status = "stabilized";
            return vh;
        }
        previous_unique_hashes = unique_hashes.size();
        if (!cfg.no_measurements)
            uniq_cnt.push_back(previous_unique_hashes -1);
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
        while (iteration < cfg.depth / 2) {
            if (const auto result = check_progress())
                return *result;
            iteration_step();
        }
        if (iteration >= cfg.depth / 2) {
            status = "depth_reached";
        }
        Hash final = cfg.depth % 2 == 0 ? variable_hash() : cnf_hash();


        return final;
    }
};

// could be deleted => keeping it for now as backup
inline std::string weisfeiler_leman_hash(
    const char* filename,
    unsigned depth = 100,
    bool cross_reference_literals = true,
    bool rehash_clauses = true,
    bool optimize_first_iteration = true,
    unsigned progress_check_iteration = 1,
    // bool shrink_to_fit = false,
    bool no_measurements = false,
    bool sort_for_clause_hash = false,
    bool use_xxh3 = true,
    std::optional<unsigned> prime_ring_modulus = std::nullopt //std::nullopt
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
        no_measurements,
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