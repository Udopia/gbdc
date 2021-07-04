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

#ifndef CANDY_GATES_BLOCKLIST_H_
#define CANDY_GATES_BLOCKLIST_H_

#include <vector>
#include <set>
#include <limits>

#include "util/CNFFormula.h"

class BlockList {

private:
    const CNFFormula& problem;

    std::vector<For> index;
    std::vector<Cl*> unitc;
    std::vector<uint16_t> num_blocked;

    #define CLAUSES_ARE_SORTED 
#ifdef CLAUSES_ARE_SORTED
    bool isBlocked(Lit o, Cl& c1, Cl& c2) const { // assert o \in c1 and ~o \in c2
        for (unsigned i = 0, j = 0; i < c1.size() && j < c2.size(); c1[i] < c2[j] ? ++i : ++j) {
            if (c1[i] != o && c1[i] == ~c2[j]) return true;
        }
        return false;
    }
#else
    bool isBlocked(Lit o, Cl& c1, Cl& c2) const { // assert o \in c1 and ~o \in c2
        for (Lit l1 : c1) if (l1 != o) for (Lit l2 : c2) if (l1 == ~l2) return true;
        return false;
    }
#endif

    bool isBlocked(Lit o, Cl* clause) const { // assert o \in clause
        for (Cl* c2 : index[~o]) if (!isBlocked(o, *clause, *c2)) return false;
        return true;
    }

    void initBlockingCounter(Lit o) {
        int i = 0; 
        int j = index[o].size()-1;
        while (i <= j) {
            //std::cout << i << " <= " << j << " < " << index[o].size() << std::endl;
            if (isBlocked(o, index[o][i])) {
                ++i;
            }
            else {
                if (i < j) std::swap(index[o][i], index[o][j]);
                --j;
            }
        }
        num_blocked[o] = i;
        if (num_blocked[o] == index[o].size()) {
            num_blocked[~o] == index[~o].size();
        }
    }

public:
    BlockList(const CNFFormula& problem_) : problem(problem_), unitc() { 
        index.resize(2 + 2 * problem.nVars());
        num_blocked.resize(2 + 2 * problem.nVars(), 0);

        for (Cl* clause : problem_) {
            if (clause->size() == 1) {
                unitc.push_back(clause);
            }
            else {
                for (Lit lit : *clause) {
                    index[lit].push_back(clause);
                }
            }
        }
    }

    ~BlockList() { }

    void remove(Var o) {
        std::set<Lit> literals;
        for (Lit olit : { Lit(o, false), Lit(o, true) }) {
            for (Cl* clause : index[olit]) {
                for (Lit lit : *clause) {
                    if (lit != olit) {
                        unsigned pos = 0;
                        for (auto it = index[lit].begin(); it < index[lit].end(); it++, pos++) {
                            if (*it == clause) {
                                index[lit].erase(it);
                                break;
                            }
                        }
                        if (pos < num_blocked[lit]) { // removed clause was blocked by lit
                            --num_blocked[lit];
                            if (num_blocked[lit] == index[lit].size()) {
                                num_blocked[~lit] == index[~lit].size();
                            }
                        }
                        else if (num_blocked[lit] == index[lit].size()) {
                            num_blocked[~lit] == index[~lit].size();
                        }
                        else {
                            literals.insert(~lit);
                        }
                    }
                }
            }
            index[olit].clear();
            num_blocked[olit] = 0;
        }
        for (Lit lit : literals) {
            if (num_blocked[lit] != index[lit].size()) initBlockingCounter(lit);
        }
    }

    inline const For& operator [](size_t o) const {
        return index[o];
    }

    inline size_t size() const {
        return index.size();
    }

    inline bool isBlockedSet(Lit o) {
        return index[o].size() == num_blocked[o];
    }

    For estimateRoots() {
        For result {};

        if (unitc.size() > 0) {
            std::swap(result, unitc);
        }
        else {
            Lit lit = getMinimallyUnblockedLiteral();
            if (lit != lit_Undef) {
                result = stripUnblockedClauses(lit);
            }
        }

        for (Cl* c : result) for (Lit l : *c) if (num_blocked[l] == 0) initBlockingCounter(l);

        return result;
    }

    Lit getMinimallyUnblockedLiteral() {
        Lit result = lit_Undef;
        uint16_t min = std::numeric_limits<uint16_t>::max();
        for (int v = problem.nVars()-1; v >= 0 && min > 1; v--) {
            for (Lit lit : { Lit(v, true), Lit(v, false) }) {
                size_t total = index[lit].size();
                if (num_blocked[lit] == 0) {
                    initBlockingCounter(lit);
                }
                size_t diff = total - num_blocked[lit];
                if (diff > 0 && diff < min) {
                    min = (uint16_t)diff;
                    result = lit;
                }
            }
        }
        return result;
    }

    For stripUnblockedClauses(Lit o) {
        For result;
        
        for (Cl* clause : index[o]) {
            if (!isBlocked(o, clause)) {
                result.push_back(clause);
            }
        }

        for (Cl* clause : result) {
            for (Lit lit : *clause) {
                For& h = index[lit];
                h.erase(std::remove(h.begin(), h.end(), clause), h.end());
                if (lit != o) {
                    num_blocked[lit] = 0;
                    num_blocked[~lit] = 0;
                }
                else {
                    num_blocked[~lit] = index[~lit].size();
                }
            }
        }

        return result;
    }

};

#endif