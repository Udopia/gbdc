#include "UnionFind.h"
#include <algorithm>
#include <vector>


UnionFind::UnionFind() : ccs() {}

void UnionFind::insert(const Cl &cl)
{
    Var min_var = Var(cl.front().var()), par;
    for (const Lit &lit : cl) {
        par = find(lit.var());
        if (min_var > par){
            ccs[min_var] = par;
            min_var = par;
        } else {
            ccs[par] = min_var;
        }
    }
}

Var UnionFind::find(Var var)
{
    return var == ccs[var] ? var : (ccs[var] = find(ccs[var]));
}

unsigned UnionFind::count_components()
{
    unsigned num_components = 0;
    for (unsigned i = 1; i < ccs.size(); ++i)
    {
        num_components += i == static_cast<unsigned>(find(ccs[i]));
    }
    return num_components;
}