#ifndef SRC_UTIL_CSRCNFFORMULA_H_
#define SRC_UTIL_CSRCNFFORMULA_H_

#include <vector>
#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

class CSRCNFFormula {
    std::vector<Lit> lits_;
    std::vector<uint32_t> start_;
    uint32_t nVars_ = 0;
    uint32_t nLits_ = 0;

 public:
    explicit CSRCNFFormula(const char* filename) {
        readDimacsFromFile(filename);
        normalizeVariableNames();
    }

    struct Clause {
        const Lit* begin_;
        const Lit* end_;
        const Lit* begin() const {
            return begin_;
        };
        const Lit* end() const {
            return end_;
        };
        uint32_t size() const {
            return end_ - begin_;
        };
    };

    struct ClauseIterator {
        const CSRCNFFormula* cnf;
        uint32_t idx;

        bool operator!=(ClauseIterator o) const {
            return idx != o.idx;
        }
        ClauseIterator& operator++() { 
            ++idx; return *this; 
        }
        Clause operator*() const {
            return { 
                cnf->lits_.data() + cnf->start_[idx], 
                cnf->lits_.data() + cnf->start_[idx+1] 
            };
        }
    };

    struct ClausesRange {
        ClauseIterator b,e;
        ClauseIterator begin() const { 
            return b; 
        }
        ClauseIterator end() const { 
            return e; 
        }
    };
    ClausesRange clauses() const {
        return {
            {this, 0}, 
            {this, static_cast<uint32_t>(start_.size()-1)}
        };
    }

    // get useful statistics
    uint32_t nVars() const {
        return nVars_;
    }
    uint32_t nClauses() const {
        return start_.size()-1;
    }
    uint32_t nLiterals() const {
        return nLits_;
    }

    // create gapless representation of variables
    void normalizeVariableNames() {
        std::vector<unsigned> map(nVars_+1, UINT32_MAX);
        unsigned next = 0;
        for (Lit& lit : lits_) {
            unsigned v = lit.var().id;
            if (map[v] == UINT32_MAX) map[v] = next++;
            lit = Lit(map[v], lit.sign());
        }
        nVars_ = next;
    }

    void readDimacsFromFile(const char* filename)
    {
        StreamBuffer in(filename);
        std::vector<Lit> clause;

        start_.push_back(0);

        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
                continue;
            }

            int plit;
            while (in.readInteger(&plit)) {
                if (plit == 0) break;
                uint32_t absID = std::abs(plit);
                Var v(absID);
                clause.emplace_back(v, plit < 0);

                nVars_ = std::max(nVars_, absID);
            }

            // store clause
            lits_.insert(lits_.end(), clause.begin(), clause.end());
            nLits_ += clause.size();
            start_.push_back(lits_.size());   // begin index of the next clause

            clause.clear();
        }
    }
};

#endif  // SRC_UTIL_CSRCNFFORMULA_H_
