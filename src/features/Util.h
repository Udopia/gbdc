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

template <typename T>
float Mean(std::vector<T> distribution) {
    float sum = static_cast<float>(std::accumulate(distribution.begin(), distribution.end(), 0));
    return sum / distribution.size();
}

template <typename T>
float Variance(std::vector<T> distribution, float mean) {
    float sum = static_cast<float>(std::accumulate(distribution.begin(), distribution.end(), 0.0,
        [mean] (float a, unsigned b) { return a + pow(static_cast<float>(b - mean), 2); } ));
    return sum / distribution.size();
}

float Entropy(std::vector<unsigned> distribution) {
    float entropy = 0;
    std::vector<unsigned> frequency;
    for (unsigned value : distribution) {
        if (value >= frequency.size()) {
            frequency.resize(value, 0);
        }
        ++frequency[value];
    }
    for (unsigned freq : frequency) if (freq > 0) {
        float p_x = static_cast<float>(freq) / distribution.size();
        if (p_x > 0) entropy -= p_x * log(p_x) / log(2);
    }
    return entropy;
}

float Entropy(std::vector<float> distribution) {
    float entropy = 0;
    std::vector<unsigned> frequency;
    for (float value : distribution) {
        if (value >= frequency.size()) {
            frequency.resize(value, 0);
        }
        ++frequency[static_cast<unsigned>(std::round(10 * value))];
    }
    for (unsigned freq : frequency) if (freq > 0) {
        float p_x = static_cast<float>(freq) / distribution.size();
        if (p_x > 0) entropy -= p_x * log(p_x) / log(2);
    }
    return entropy;
}


template <typename T>
void push_distribution(std::vector<float>* record, std::vector<T> distribution) {
    float mean = 0, variance = 0, min = 0, max = 0, entropy = 0;
    if (distribution.size() > 0) {
        mean = Mean(distribution);
        variance = Variance(distribution, mean);
        min = static_cast<float>(*std::min_element(distribution.begin(), distribution.end()));
        max = static_cast<float>(*std::max_element(distribution.begin(), distribution.end()));
        entropy = Entropy(distribution);
    }
    record->push_back(mean);
    record->push_back(variance);
    record->push_back(min);
    record->push_back(max);
    record->push_back(entropy);
}

#endif  // SRC_FEATURES_UTIL_H_
