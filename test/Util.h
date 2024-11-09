#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <random>
#include <cstdlib>
#include <algorithm>
#include <cstdint>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

static std::string tmp_filename(std::string dir, std::string ext, unsigned length = 32U)
{
    const char hex_chars[] = "0123456789abcdef";
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 15);

    std::string filename;
    filename.reserve(length);

    for (std::size_t i = 0; i < length; ++i)
    {
        filename += hex_chars[dis(gen)];
    }

    return dir + "/" + filename + ext;
}

template <typename Entry>
static bool has_extension(Entry entry, std::string ex)
{
    std::string filename = entry.path().filename().string();
    ex = "." + ex + ".";
    return (filename.find(ex) != std::string::npos);
}

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

template <typename Container>
static void check_subset(Container &&subset, Container &&superset)
{
    for (size_t i = 0; i < subset.size(); ++i)
    {
        CHECK((subset[i] == superset[i]));
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
