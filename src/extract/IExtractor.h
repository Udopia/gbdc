/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>

class IExtractor {
protected:
    std::unordered_map<std::string, double> dict;
    std::vector<std::string> names;

    void initFeatures(std::initializer_list<std::string> featureNames) {
        for (const auto& name : featureNames) {
            setFeature(name, 0.0);    
        }
    }

    template<typename Iterator>
    void initFeatures(Iterator begin, Iterator end) {
        for (auto it = begin; it != end; ++it) {
            setFeature(*it, 0.0);
        }
    }

    template<typename T>
    void setFeature(const std::string& name, T value) {
        if (dict.find(name) == dict.end()) {
            names.push_back(name);
        }
        dict[name] = static_cast<double>(value);
    }

    template<typename Iterator>
    void setFeatures(std::initializer_list<std::string> featureNames, Iterator begin, Iterator end) {
        const std::string* nameiter = featureNames.begin();
        for (auto it = begin; it != end; ++it, ++nameiter) {
            setFeature(*nameiter, *it);
        }
    }

public:
    virtual ~IExtractor() { }
    virtual void run() = 0;

    std::vector<std::string> getNames() const {
        return names;
    }

    std::vector<double> getFeatures() const {
        std::vector<double> features;
        for (const auto& name : names) {
            features.push_back(dict.at(name));
        }
        return features;
    }

    double getFeature(const std::string& name) const {
        return dict.at(name);
    }
};