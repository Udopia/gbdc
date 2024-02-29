/**
 * MIT License
 * 
 * Copyrigth (c) 2023 Markus Iser 
 */

#ifndef OPB_BASE_FEATURES_H_
#define OPB_BASE_FEATURES_H_

#include "IExtractor.h"
#include "src/util/StreamBuffer.h"
#include "src/extract/Util.h"

namespace OPB {
    
class Constr;
class BaseFeatures;
    
class TermSum {
    friend Constr;
    friend BaseFeatures;

    std::vector<int> coeffs{};
    int max = 0;
    int min = 0;
    int abs_min_coeff = std::numeric_limits<int>::max();
    Var max_var{0};

  public:
    TermSum(StreamBuffer &in) {
        for (in.skipWhitespace(); *in != ';' && *in != '>' && *in != '='; in.skipWhitespace()) {
            int coeff;
            in.readInteger(&coeff);
            in.skipWhitespace();
            if (*in == 'x') {
                in.skip();
            } else {
                assert(*in == '~');
                in.skip();
                in.skipWhitespace();
                in.skip();
            }
            if (coeff < 0) {
                min += coeff;
            } else {
                max += coeff;
            }
            abs_min_coeff = std::min(std::abs(coeff), abs_min_coeff);
            int var;
            in.readInteger(&var);
            if (Var(var + 1) > max_var) max_var = Var(var + 1);
            coeffs.push_back(coeff);
        }
    }
    
    inline size_t nTerms() {
        return coeffs.size();
    }
    
    inline int const minVal() {
        return min;
    }
    
    inline int const maxVal() {
        return max;
    }
    
    inline Var const maxVar() {
        return max_var;
    }
    
    inline int const minCoeff() {
        return abs_min_coeff;
    }
};

class Constr {
    friend BaseFeatures;
    
  public:
    enum Rel { GE, EQ };
    struct Analysis {
        bool tautology : 1;
        bool unsat : 1;
        bool assignment : 1;
        bool clause : 1;
        bool card : 1;
    };

  private:
    TermSum terms;
    Rel rel;
    int bound;

  public:
    Constr(StreamBuffer &in) : terms(in) {
        if (*in == '>') {
            rel = GE;
            in.skipString(">=");
        } else {
            assert(*in == '=');
            in.skip();
        }
        in.readInteger(&bound);
        in.skipWhitespace();
        if (*in == ';') in.skip();
    }
    
    inline Analysis analyse() {
        Analysis a {};
        if (terms.nTerms()) {
            int multiplier = abs(terms.coeffs.front());
            a.card = true;
            for (int coeff: terms.coeffs) {
                if (std::abs(coeff) != multiplier) {
                    a.card = false;
                    break;
                }
            }
        }
        switch (rel) {
            case GE:
                a.tautology = terms.minVal() >= bound;
                a.unsat = terms.maxVal() < bound;
                a.assignment = terms.maxVal() - terms.minCoeff() < bound && terms.maxVal() > bound;
                a.clause = bound > terms.minVal() && bound <= terms.minVal() + terms.minCoeff();
                break;
            case EQ:
            default:
                a.tautology = terms.minVal() == terms.maxVal() && terms.minVal() == bound;
                a.unsat = terms.minVal() > bound || terms.maxVal() < bound;
                a.assignment = bound == terms.maxVal() || bound == terms.minVal();
                a.clause = false;
        }

        return a;
    }
    
    inline Var maxVar() {
        return terms.maxVar();
    }
};

class BaseFeatures : public IExtractor {
    const char* filename_;
    std::vector<double> features;
    std::vector<std::string> names;

    unsigned n_vars = 0, n_constraints = 0;
    unsigned n_pbs_ge = 0, n_pbs_eq = 0;
    unsigned n_cards_ge = 0, n_cards_eq = 0;
    unsigned n_clauses = 0, n_assignments = 0;
    bool trivially_unsat = false;
    
    unsigned obj_terms = 0;
    int obj_max_val = 0, obj_min_val = 0;
    std::vector<int> obj_coeffs{};

  public:
    BaseFeatures(const char* filename) : filename_(filename), features(), names() { 
        names.insert(names.end(), { "constraints", "variables" });
        names.insert(names.end(), { "pbs_ge", "pbs_eq", "cards_ge", "cards_eq" });
        names.insert(names.end(), { "clauses", "assignments", "trivially_unsat" });
        names.insert(names.end(), { "obj_terms", "obj_max_val", "obj_min_val" });
        names.insert(names.end(), { "obj_coeffs_mean", "obj_coeffs_variance" });
        names.insert(names.end(), { "obj_coeffs_min", "obj_coeffs_max", "obj_coeffs_entropy" });
    }

    virtual ~BaseFeatures() { }

    virtual void extract() {
        StreamBuffer in(filename_);

        bool seen_obj = false;
        while (in.skipWhitespace()) {
            if (*in == '*') {
                in.skipLine();
            } else if (*in == 'm') {
                in.skipString("min:");
                // if multiple objective lines are encountered, the first will be used
                if (seen_obj) {
                    in.skipLine();
                    continue;
                }
                seen_obj = true;
                TermSum obj(in);
                obj_terms = obj.nTerms();
                obj_max_val = obj.maxVal();
                obj_min_val = obj.minVal();
                obj_coeffs = obj.coeffs;
                if (obj.maxVar() > n_vars) n_vars = obj.maxVar();
                in.skipWhitespace();
                if (*in == ';') in.skip();
            } else {
                n_constraints++;
                
                Constr constr(in);
                if (constr.maxVar() > n_vars) n_vars = constr.maxVar();
                Constr::Analysis a = constr.analyse();
                if (a.unsat) {
                    trivially_unsat = true;
                }
                if (a.assignment) {
                    n_assignments++;
                }
                if (a.clause) {
                    n_clauses++;
                } else if (a.card) {
                    switch (constr.rel) {
                        case Constr::GE:
                            n_cards_ge++;
                            break;                        
                        case Constr::EQ:
                            n_cards_eq++;
                    }
                } else {
                    switch (constr.rel) {
                        case Constr::GE:
                            n_pbs_ge++;
                            break;
                        case Constr::EQ:
                            n_pbs_eq++;
                    }
                }
            }
        }
        
        load_feature_record();
    }

    void load_feature_record() {
        features.insert(features.end(), { (double)n_constraints, (double)n_vars });
        features.insert(features.end(), { (double)n_pbs_ge, (double)n_pbs_eq });
        features.insert(features.end(), { (double)n_cards_ge, (double)n_cards_eq });
        features.insert(features.end(), { (double)n_clauses, (double)n_assignments });
        features.insert(features.end(), { (double)trivially_unsat });
        features.insert(features.end(), { (double)obj_terms, (double)obj_max_val, (double)obj_min_val });
        push_distribution(features, obj_coeffs);
    }

    virtual std::vector<double> getFeatures() const {
        return features;
    }

    virtual std::vector<std::string> getNames() const {
        return names;
    }
};

}; // namespace OPB


#endif // OPB_BASE_FEATURES_H_