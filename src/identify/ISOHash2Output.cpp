#include "ISOHash2.h"
#include <ostream>
#include <sstream>
#include <fstream>

namespace CNF {

void write_csv(const WLHResult& r, std::ostream& out) {
    static bool wrote_header = false;
    if (!wrote_header) {
        out << "filename,hash,iterations,total_runtime,status,"
               "iterations_time,unique_variables,"
               "nVars,nClauses,nLits,maxClauseLen,"
               "depth,cross_ref,rehash,opt_first,"
               "prog_iter,sort,use_xxh3,prime_ring\n";
        wrote_header = true;
    }

    auto join = [&](auto const& v) {
        std::ostringstream s;
        for (size_t i = 0; i < v.size(); ++i) {
            if (i) s << '|';
            s << v[i];
        }
        return s.str();
    };

    out  << r.filename << ',' << r.hash << ',' << r.iterations << ',' << r.total_runtime << ',' << r.status << ','
         << '"' << join(r.iterations_time) << '"' << ','
         << '"' << join(r.unique_variables) << '"' << ','
         << r.nVars << ',' << r.nClauses << ',' << r.nLits << ',' << r.maxClauseLen << ','
         << r.cfg.depth << ',' << (r.cfg.cross_reference_literals ? "yes" : "no") << ',' << (r.cfg.rehash_clauses ? "yes" : "no") << ',' << (r.cfg.optimize_first_iteration ? "yes" : "no") << ','
         << r.cfg.progress_check_iteration << ',' << (r.cfg.sort_for_clause_hash ? "yes" : "no") << ',' << (r.cfg.use_xxh3 ? "yes" : "no") << ','
         << (r.cfg.prime_ring_modulus
                 ? std::to_string(*r.cfg.prime_ring_modulus)
                 : "0")
         << '\n';
}

} // namespace CNF
