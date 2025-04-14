/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cstddef>

template <typename T> 
double Mean(std::vector<T> distribution);

template <typename T> 
double Variance(std::vector<T> distribution, double mean);

double ScaledEntropyFromOccurenceCounts(std::unordered_map<int64_t, int64_t> occurence, size_t total);
double ScaledEntropy(std::vector<double> distribution);

template <typename T> 
double ScaledEntropy(std::vector<T> distribution);

template <typename T>
std::vector<double> getDistributionStats(std::vector<T> distribution);