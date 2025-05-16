/**
 * MIT License
 * Copyright (c) 2024 Markus Iser
 */

#include "CNFSaniCheck.h"

/**
 * @brief Determines if a CNF formula is sanitised, normalised, etc.
 */
void CNF::SaniCheck::run() {
    checkNormalised();
    if (sanicheck) {
        checkSanitised();
    }
}

void CNF::SaniCheck::checkNormalised() {
    StreamBuffer in(filename_);

    int head_vars = 0;
    int head_clauses = 0;
    int norm_vars = 0;
    int norm_clauses = 0;
    bool normalised = true;
    bool comment = false;
    unsigned count = 0;
    bool start = true;

    while (count = in.skipAndCountWhitespace() || start) {
        start = false;
        normalised &= (1 == count);
        if (*in == 'p') {
            in.skipString("p");
            normalised &= (*in == ' ' && 1 == in.skipAndCountWhitespace());
            in.skipString("cnf");
            normalised &= (*in == ' ' && 1 == in.skipAndCountWhitespace());
            in.readInteger(&head_vars);
            normalised &= (*in == ' ' && 1 == in.skipAndCountWhitespace());
            in.readInteger(&head_clauses);
            normalised &= (*in == '\n');
        } else if (*in == 'c') {
            comment = true;
            while (*in != '\n' && *in != '\r') {
                if (!in.skip()) break;
            }
        } else {
            normalised &= isdigit(*in) || *in == '-';
            int plit;
            int len = 0;
            while (in.readInteger(&plit)) {
                if (plit == 0) break;
                ++len;
                norm_vars = std::max(abs(plit), norm_vars);
                normalised &= (*in == ' ' && 1 == in.skipAndCountWhitespace());
            }
            if (len > 0) norm_clauses++;
            normalised &= (*in == '\n');
        }
    }

    setFeature("head_vars", head_vars);
    setFeature("head_clauses", head_clauses);
    setFeature("norm_vars", norm_vars);
    setFeature("norm_clauses", norm_clauses);
    setFeature("whitespace_normalised", normalised ? 1.0 : 0.0);
    setFeature("has_comment", comment ? 1.0 : 0.0);
}

void CNF::SaniCheck::checkSanitised() {
    int norm_vars = (unsigned)getFeature("norm_vars");
    int sani_vars = 0;
    int sani_clauses = 0;
    bool has_taut = false;
    bool has_dupl = false;
    bool has_empty = false;

    StreamBuffer in2(filename_);
    unsigned *mask = (unsigned *)calloc(norm_vars * 2 + 2, sizeof(unsigned));
    mask += norm_vars + 1;

    unsigned stamp = 0;
    while (in2.skipWhitespace()) {
        if (*in2 == 'c' || *in2 == 'p') {
            if (!in2.skipLine()) break;
        } else {
            bool tautological = false;
            int clausemax = 0;
            int plit;
            ++stamp;
            while (in2.readInteger(&plit)) {
                if (abs(plit) > norm_vars) {
                    throw ParserException(
                        std::string(filename_) + ": variable " +
                        std::to_string(abs(plit)) + " out of range");
                }
                if (plit == 0) break;
                if (mask[-plit] == stamp) {
                    tautological = true;
                    has_taut = true;
                    break;
                } else if (mask[plit] != stamp) {
                    mask[plit] = stamp;
                    clausemax = std::max(abs(plit), clausemax);
                } else {
                    has_dupl = true;
                }
            }
            if (!tautological) {
                ++sani_clauses;
                if (clausemax == 0) {
                    has_empty = true;
                } else {
                    sani_vars = std::max(clausemax, sani_vars);
                }
            } else {
                in2.skipLine();
            }
        }
    }

    setFeature("sani_vars", sani_vars);
    setFeature("sani_clauses", sani_clauses);
    setFeature("has_tautological_clause", has_taut ? 1.0 : 0.0);
    setFeature("has_duplicate_literals", has_dupl ? 1.0 : 0.0);
    setFeature("has_empty_clause", has_empty ? 1.0 : 0.0);
}