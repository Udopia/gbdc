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
    uint32_t maxClauseLen_ = 0;

 public:
    // constructor
    explicit CSRCNFFormula(const char* filename, bool shrink_to_fit = false) {
        readDimacsFromFile(filename);
        normalizeVariableNames();
        canonicalise();
        if (shrink_to_fit) {
            lits_.shrink_to_fit();
            start_.shrink_to_fit();
        }
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
    uint32_t maxClauseLength() const { 
        return maxClauseLen_; 
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

            // store clause and modify stats
            lits_.insert(lits_.end(), clause.begin(), clause.end());
            nLits_ += clause.size();
            maxClauseLen_  = std::max<uint32_t>(maxClauseLen_, clause.size());
            start_.push_back(lits_.size());

            clause.clear();
        }
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

    void canonicalise()
    {
        // setup
        std::vector<Lit>&  L = lits_;
        std::vector<uint32_t>& S = start_;

        // prepare
        uint32_t out = 0;
        std::vector<uint32_t> newStart; // fresh CSR index array
        newStart.reserve(S.size());
        newStart.push_back(0);

        // temporary container for one clause at a time
        std::vector<Lit> buf;
        // find largest clause size to reserve capacity
        uint32_t maxClauseLen = 0;
        for (uint32_t i = 0; i + 1 < S.size(); ++i)
            maxClauseLen = std::max(maxClauseLen, S[i+1] - S[i]);
        buf.reserve(maxClauseLen);

        // custom comparator
        auto litLess = [](const Lit& a, const Lit& b) {
            return (a.var().id < b.var().id) || 
                (a.var().id == b.var().id && a.sign() < b.sign());
        };

        // main loop
        for (uint32_t c = 0; c + 1 < S.size(); ++c) {
            const uint32_t begin = S[c], end = S[c+1];

            // sort
            buf.assign(L.begin() + begin, L.begin() + end);
            std::sort(buf.begin(), buf.end(), litLess);

            // tautology
            bool taut = false;
            uint32_t w = 0;
            for (uint32_t r = 0; r < buf.size(); ++r) {
                if (w && buf[r] == buf[w-1])
                    continue;
                if (w && buf[r].var() == buf[w-1].var()) {
                    taut = true;
                    break;
                }
                buf[w++] = buf[r];
            }
            buf.resize(w);

            // drop the whole clause
            if (taut || buf.empty())
                continue;

            // copy cleaned clause back into L
            std::copy(buf.begin(), buf.end(), L.begin() + out);
            out += static_cast<uint32_t>(buf.size());
            newStart.push_back(out);
        }

        // finalize
        L.resize(out);
        S.swap(newStart);
        nLits_ = out;
    }
};

#endif  // SRC_UTIL_CSRCNFFORMULA_H_
