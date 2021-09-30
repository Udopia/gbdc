/*************************************************************************************************
CNFTools -- Copyright (c) 2015, Markus Iser, KIT - Karlsruhe Institute of Technology

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

#ifndef SRC_UTIL_CNFFORMULA_H_
#define SRC_UTIL_CNFFORMULA_H_

#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "src/StreamBuffer.h"
#include "src/util/SolverTypes.h"

class CNFFormula {
    std::shared_ptr<For> formula;
    unsigned variables;

 public:
    CNFFormula() : formula(new For()), variables(0) { }

    explicit CNFFormula(const For& formula) : formula(), variables(0) {
        readClauses(formula);
    }

    explicit CNFFormula(Cl& clause) : formula(), variables(0) {
        readClause(clause.begin(), clause.end());
    }

    // copy constructor using shared-ptr to same formula
    CNFFormula(const CNFFormula& other) : formula(other.formula), variables(other.variables) { }

    ~CNFFormula() { }

    typedef For::const_iterator const_iterator;

    inline const_iterator begin() const {
        return formula->begin();
    }

    inline const_iterator end() const {
        return formula->end();
    }

    inline const Cl* operator[] (int i) const {
        return (*formula)[i];
    }

    inline size_t nVars() const {
        return variables;
    }

    inline size_t nClauses() const {
        return formula->size();
    }

    inline int newVar() {
        return ++variables;
    }

    inline void clear() {
        formula->clear();
    }

    // create gapless representation of variables
    void normalizeVariableNames() {
        std::vector<unsigned> name;
        name.resize(variables+1, 0);
        unsigned int max = 0;
        for (Cl* clause : *formula) {
            for (Lit& lit : *clause) {
                if (name[lit.var()] == 0) name[lit.var()] = max++;
                lit = Lit(name[lit.var()], lit.sign());
            }
        }
        variables = max;
    }

    void readDimacsFromFile(const char* filename) {
        StreamBuffer in(filename);
        Cl clause;
        while (!in.eof()) {
            in.skipWhitespace();
            if (in.eof()) {
                break;
            }
            if (*in == 'p' || *in == 'c') {
                in.skipLine();
            } else {
                for (int plit = in.readInteger(); plit != 0; plit = in.readInteger()) {
                    clause.push_back(Lit(abs(plit), plit < 0));
                }
                readClause(clause.begin(), clause.end());
                clause.clear();
            }
        }
    }

    void readClause(std::initializer_list<Lit> list) {
        readClause(list.begin(), list.end());
    }

    void readClauses(const For& formula) {
        for (Cl* clause : formula) {
            readClause(clause->begin(), clause->end());
        }
    }

    template <typename Iterator>
    void readClause(Iterator begin, Iterator end) {
        Cl* clause = new Cl { begin, end };
        if (clause->size() > 0) {
            // remove redundant literals
            std::sort(clause->begin(), clause->end());
            clause->erase(std::unique(clause->begin(), clause->end()), clause->end());
            // skip tautologies
            bool tautology = clause->end() != std::unique(clause->begin(), clause->end(), [](Lit l1, Lit l2) {
                return l1.var() == l2.var();
            });
            if (tautology) {
                delete clause;
                return;
            }
            // record maximal variable
            variables = std::max(variables, (unsigned int)clause->back().var());
        }
        formula->push_back(clause);
    }
};

#endif  // SRC_UTIL_CNFFORMULA_H_

