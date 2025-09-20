/**
 * MIT License
 * Copyright (c) 2025 Ashlin Iser 
 */

#pragma once

#include <iostream>
#include <fstream>
#include <memory>
#include <string>

class OutputWrapper {
public:
    // Constructor takes optional path; nullptr or empty string => std::cout
    explicit OutputWrapper(const std::string* output_path = nullptr) {
        if (output_path && !output_path->empty()) {
            stream_ = std::shared_ptr<std::ostream>(
                new std::ofstream(*output_path, std::ofstream::out),
                [](std::ostream* p) { delete p; }
            );
        } else {
            stream_ = std::shared_ptr<std::ostream>(
                &std::cout,
                [](std::ostream*) { /* do nothing */ }
            );
        }
    }

    // Templated stream insertion operator
    template<typename T>
    OutputWrapper& operator<<(const T& value) {
        (*stream_) << value;
        return *this;
    }

    // Support manipulators like std::endl
    OutputWrapper& operator<<(std::ostream& (*manip)(std::ostream&)) {
        (*stream_) << manip;
        return *this;
    }

    // Optionally expose the underlying stream (read-only)
    std::ostream& stream() const {
        return *stream_;
    }

private:
    std::shared_ptr<std::ostream> stream_;
};
