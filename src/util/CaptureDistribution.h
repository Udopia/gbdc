/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cstddef>

template <typename Container> 
double Mean(Container&& distribution);

template <typename Container> 
double Variance(Container&& distribution, double mean);

double ScaledEntropyFromOccurenceCounts(std::unordered_map<int64_t, int64_t> occurence, size_t total);
double ScaledEntropy(std::vector<double> distribution);

template <typename Container> 
double ScaledEntropy(Container&& distribution);

template <typename V, typename W> 
void push_distribution(V&& record, W&& distribution);