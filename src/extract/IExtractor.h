/**
 * MIT License
 *
 * Copyrigth (c) 2023 Markus Iser
 */

#ifndef EXTRACTOR_INTERFACE_H_
#define EXTRACTOR_INTERFACE_H_

#include <string>
#include <vector>
#include "src/util/threadpool/TrackingAllocator.h"

template <typename T>
using talloc = TrackingAllocator<T>;

// template <typename T>
// using tvector = std::vector<T, talloc<T>>;

template <typename T>
class tvector : public std::vector<T, talloc<T>>
{
public:
    using Base = std::vector<T, talloc<T>>; // Type alias for the base class

    // Inherit constructors from the base class
    using Base::Base;

    tvector(const std::vector<T> &v) : Base(v.begin(), v.end()) {}

    template <typename Container>
    tvector &operator=(const Container &v)
    {
        static_assert(std::is_same_v<T, typename Container::value_type>,
                      "Container value type must match tvector value type");
        this->clear();
        this->reserve(v.size());
        for (const auto u : v)
        {
            this->push_back(u);
        }
        return *this;
    }
};

class IExtractor
{
public:
    virtual ~IExtractor() {}
    virtual void extract() = 0;
    virtual std::vector<double> getFeatures() const = 0;
    virtual std::vector<std::string> getNames() const = 0;
};

#endif // EXTRACTOR_INTERFACE_H_