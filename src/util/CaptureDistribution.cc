/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#include "CaptureDistribution.h"

#include <algorithm>
#include <cmath>

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
            occurence[snap] = occurence[snap] + 1;
        } else {
            occurence[snap] = 1;
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

// Explicit template instantiations
template void push_distribution<std::vector<double, std::allocator<double> >&, std::vector<unsigned int, std::allocator<unsigned int> >&>(std::vector<double, std::allocator<double> >&, std::vector<unsigned int, std::allocator<unsigned int> >&);
template void push_distribution<std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&>(std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&);
template void push_distribution<std::vector<double, std::allocator<double> >&, std::vector<unsigned long, std::allocator<unsigned long> >&>(std::vector<double, std::allocator<double> >&, std::vector<unsigned long, std::allocator<unsigned long> >&);