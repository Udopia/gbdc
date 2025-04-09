/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#ifndef EXTRACTOR_INTERFACE_H_
#define EXTRACTOR_INTERFACE_H_

#include <string>
#include <vector>

class IExtractor {
public:
    virtual ~IExtractor() { }
    virtual void run() = 0;
    virtual std::vector<double> getFeatures() const = 0;
    virtual std::vector<std::string> getNames() const = 0;
};

#endif // EXTRACTOR_INTERFACE_H_