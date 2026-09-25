#pragma once

#include <vector>

#include "src/extract/IExtractor.h"
#include "src/util/CNFFormula.h"

namespace CNF {

// Literal Incidence Graph (LIG) - encodes polarities
class cnf2lig : public IExtractor {
 private:
    CNFFormula F;
    const char* filename_;
    const char* output_;

 public:
    cnf2lig(const char* filename, const char* output = nullptr);
    virtual ~cnf2lig();
    virtual void run();
};

}  // namespace CNF