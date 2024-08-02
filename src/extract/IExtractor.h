/**
 * MIT License
 * 
 * Copyrigth (c) 2023 Markus Iser 
 */

#ifndef EXTRACTOR_INTERFACE_H_
#define EXTRACTOR_INTERFACE_H_

#include <string>
#include <vector>

class IExtractor {
public:
    virtual ~IExtractor() { }
    virtual void extract() = 0;
    virtual std::vector<double> getFeatures() const = 0;
    virtual std::vector<std::string> getNames() const = 0;
};

#endif // EXTRACTOR_INTERFACE_H_