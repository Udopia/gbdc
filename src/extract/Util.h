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

#ifndef SRC_FEATURES_UTIL_H_
#define SRC_FEATURES_UTIL_H_

#include <math.h>

#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include "src/util/threadpool/TrackingAllocator.h"

template <typename Container>
double Mean(Container&& distribution) {
    double mean = 0.0;
    for (size_t i = 0; i < distribution.size(); i++) {
        mean += (distribution[i] - mean) / (i + 1);
    }
    return mean;
}

template <typename Container>
double Variance(Container&& distribution, double mean) {
    double vari = 0.0;
    for (size_t i = 0; i < distribution.size(); i++) {
        double diff = distribution[i] - mean;
        vari += (diff*diff - vari) / (i + 1);
    }
    return vari;
}

double ScaledEntropyFromOccurenceCounts(std::unordered_map<int64_t, int64_t> occurence, size_t total) {
    // collect and sort summands
    std::vector<long double> summands;
    for (auto& pair : occurence) {
        long double p_x = (long double)pair.second / (long double)total;
        long double summand = p_x * log2(p_x);
        // long double summand = (pair.second * log2(pair.second) - pair.second * log2(total)) / total;
        summands.push_back(summand);
    }
    std::sort(summands.begin(), summands.end(), [] (long double a, long double b) { return abs(a) < abs(b); });
    // calculate entropy
    long double entropy = 0;
    for (long double summand : summands) {
        entropy -= summand;
    }
    // scale by log of number of categories    
    return log2(summands.size()) == 0 ? 0 : (double)entropy / log2(summands.size());
}

double ScaledEntropy(std::vector<double> distribution) {
    std::unordered_map<int64_t, int64_t> occurence;
    for (double value : distribution) {
        // snap to 3 digits after decimal point
        int64_t snap = static_cast<int64_t>(std::round(1000 * value));
        if (occurence.count(value)) {
            occurence[value] = occurence[value] + 1;
        } else {
            occurence[value] = 1;
        }
    }
    return ScaledEntropyFromOccurenceCounts(occurence, distribution.size());
}

template <typename Container>
double ScaledEntropy(Container&& distribution) {
    std::unordered_map<int64_t, int64_t> occurence;
    for (unsigned value : distribution) {
        if (occurence.count(value)) {
            occurence[value] = occurence[value] + 1;
        } else {
            occurence[value] = 1;
        }
    }
    return ScaledEntropyFromOccurenceCounts(occurence, distribution.size());
}

template <typename V, typename W>
void push_distribution(V&& record, W&& distribution) {
    if (distribution.size() == 0) {
        record.insert(record.end(), { 0, 0, 0, 0, 0 });
        return;
    }
    std::sort(distribution.begin(), distribution.end());
    double mean = Mean(distribution);
    double variance = Variance(distribution, mean);
    double min = distribution.front();
    double max = distribution.back();
    double entropy = ScaledEntropy(distribution);
    record.insert(record.end(), { mean, variance, min, max, entropy });
}

inline size_t numDigits(unsigned x)
{
    return ceil(log10(x));
}

class UnionFind
{
private:
    template<typename T>
    struct vwrapper
    {
        std::vector<T, TrackingAllocator<T>> v;
        T &operator[](size_t idx)
        {
            if (idx >= v.size())
            {
                auto old_end = v.size();
                v.resize(idx + 1);
                std::generate(v.begin() + old_end, v.end(), [&]
                              { return T(old_end++); });
            }
            return v[idx];
        }

        size_t size()
        {
            return v.size();
        }
    };
    vwrapper<Var> ccs;
public:
    UnionFind() : ccs() {}
    
    /**
     * @brief Inserts all variables of the given clause cl into the data structure and
     * merges their components by selecting the representative with the lowest variable ID
     * as the only remaining representative.
     * @param cl Clause for which to insert all variables.
     * @return
     */
    inline void insert(const Cl &cl)
    {
        Var min_var = Var(cl.front().var()), par;
        for (const Lit &lit : cl) {
            par = find(lit.var());
            if (min_var > par){
                ccs[min_var] = par;
                min_var = par;
            } else {
                ccs[par] = min_var;
            }
        }
    }

    inline Var find(Var var)
    {
        return var == ccs[var] ? var : (ccs[var] = find(ccs[var]));
    }

    inline unsigned count_components()
    {
        unsigned num_components = 0;
        for (unsigned i = 1; i < ccs.size(); ++i)
        {
            num_components += i == static_cast<unsigned>(find(ccs[i]));
        }
        return num_components;
    }
};

#endif // SRC_FEATURES_UTIL_H_
