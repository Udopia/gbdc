/*************************************************************************************************
GBDHash -- Copyright (c) 2020, Markus Iser, KIT - Karlsruhe Institute of Technology

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

#pragma once
#include <cstdio>
#include <cstdint>
#include <string>
#include <chrono>
#include <regex>
#include <fstream>
#include <memory>
#include <sstream>
#include <future>
#include <unordered_map>
#include <variant>
#include <filesystem>

#include "src/identify/GBDHash.h"
#include "src/identify/ISOHash.h"

#include "src/util/threadpool/ThreadPool.h"
#include "src/util/threadpool/TrackingAllocator.h"

#include "src/extract/CNFBaseFeatures.h"
#include "src/extract/CNFGateFeatures.h"
#include "src/extract/WCNFBaseFeatures.h"
#include "src/extract/OPBBaseFeatures.h"

#include "src/transform/IndependentSet.h"
#include "src/transform/Normalize.h"
#include "src/util/ResourceLimits.h"

namespace tp = threadpool;

static std::string version()
{
    std::ifstream file("setup.py");
    if (!file.is_open())
    {
        return "Error: Could not open setup.py";
    }

    std::string line;
    std::regex version_regex(R"(version\s*=\s*['"]([^'"]+)['"])");
    std::smatch match;

    while (std::getline(file, line))
    {
        if (std::regex_search(line, match, version_regex))
        {
            return match[1].str();
        }
    }

    return "Error: Version not found in setup.py";
}

template <int type_id>
constexpr auto get_runtime_desc()
{
    static_assert((type_id == 0) || (type_id == 1));
    switch (type_id)
    {
    case 0:
        return "base_features_runtime";
    case 1:
        return "gate_features_runtime";
    default:
        return "base_features_runtime";
    }
}

template <typename Extractor>
static auto tp_extract_features(const std::string filepath)
{
    Extractor stats(filepath.c_str());
    stats.extract();
    tp::extract_ret_t dict;
    const auto names = stats.getNames();
    const auto features = stats.getFeatures();
    for (size_t i = 0; i < features.size(); ++i)
    {
        dict[names[i]] = features[i];
    }
    return dict;
}

template <typename Extractor, int type_id = 0>
static auto feature_names()
{
    auto names = Extractor("").getNames();
    names.push_back(get_runtime_desc<type_id>());
    return names;
}

template <template <template <typename> typename> typename Extractor>
static auto tp_extract(const std::uint64_t _mem_max, const std::uint32_t _jobs_max, const std::vector<std::tuple<tp::extract_arg_t>> args)
{
    auto tp = std::make_shared<tp::ThreadPool<tp::extract_ret_t, tp::extract_arg_t>>(_mem_max, _jobs_max, &tp_extract_features<Extractor<TrackingAllocator>>, args);
    auto q = tp->get_result_queue();
    std::thread([tp]()
                { tp->start_threadpool(); 
                std::cerr << "Destroying Threadpool!" << std::endl; })
        .detach();
    return q;
}