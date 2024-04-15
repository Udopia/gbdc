#pragma once
#include <string>
#include <algorithm>
#include <queue>
#include <thread>
#include <filesystem>
#include <iostream>
#include <fstream>

class csv_t
{
private:
    std::ofstream of;
    std::string s = "";

public:
    csv_t(const std::string &filename, std::initializer_list<std::string> columns) : of(filename)
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error opening file " << filename << std::endl;
            return;
        }
        write_to_file(columns);
    }

    void write_to_file(std::initializer_list<std::string> data)
    {
        for (const auto &d : data)
        {
            s.append(d + ",");
        }
        s.pop_back();
        of << s << "\n";
        s.clear();
    }

    template <typename Data>
    void write_to_file(Data &&data)
    {
        for (const auto &d : data)
        {
            s.append(std::to_string(d) + ",");
        }
        s.pop_back();
        of << s << "\n";
        s.clear();
    }

    void write_to_file(std::initializer_list<size_t> data)
    {
        for (const auto &d : data)
        {
            s.append(std::to_string(d) + ",");
        }
        s.pop_back();
        of << s << "\n";
        s.clear();
    }

    void close_file()
    {
        of.close();
    }
};

class TerminationRequest : public std::runtime_error
{
public:
    size_t memnbt;
    explicit TerminationRequest(size_t _memnbt) : std::runtime_error(""), memnbt(_memnbt) {}
};