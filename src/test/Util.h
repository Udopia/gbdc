#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <random>
#include <cstdlib>
#include <algorithm>
#include "../util/threadpool/ThreadPool.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

namespace tp = threadpool;

static double string_to_double(const std::string &s)
{
    std::istringstream os(s);
    double d;
    os >> d;
    return d;
}

static bool is_number(const std::string &s)
{
    return !s.empty() && std::find_if(s.begin(),
                                      s.end(), [](unsigned char c)
                                      { return !std::isdigit(c); }) == s.end();
}


static bool fequal(const double a, const double b)
{
    const double epsilon = fmax(fabs(a), fabs(b)) * 1e-5;
    return fabs(a - b) <= epsilon;
}

static void check_subset(tp::extract_ret_t &expected, tp::extract_ret_t &actual)
{
    for (const auto [key, val] : expected)
    {
        if (val.index() == 0)
            CHECK(fequal(std::get<double>(val), std::get<double>(actual[key])));
    }
}

template <typename Val>
static std::unordered_map<std::string, Val> record_to_map(std::string record_file_name)
{
    std::ifstream record_file(record_file_name);
    std::unordered_map<std::string, Val> map;

    std::string line;
    while (std::getline(record_file, line))
    {
        std::istringstream iss(line);
        std::string key, value;

        if (std::getline(iss, key, '=') && std::getline(iss, value))
        {
            map[key] = string_to_double(value);
        }
    }

    record_file.close();

    return map;
}
