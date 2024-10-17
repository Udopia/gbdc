#ifndef UNIONFIND_H
#define UNIONFIND_H

#include <vector>
#include <algorithm>

#include "src/util/SolverTypes.h"

class UnionFind
{
private:
    struct vwrapper
    {
        std::vector<Var> v;
        Var &operator[](size_t idx)
        {
            if (idx >= v.size())
            {
                auto old_end = v.size();
                v.resize(idx + 1);
                std::generate(v.begin() + old_end, v.end(), [&]
                              { return Var(old_end++); });
            }
            return v[idx];
        }

        size_t size()
        {
            return v.size();
        }
    };

    vwrapper ccs;

public:
    UnionFind();
    void insert(const Cl &cl);
    Var find(Var var);
    unsigned count_components();
};

#endif // UNIONFIND_H