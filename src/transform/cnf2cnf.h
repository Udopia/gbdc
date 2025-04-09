/**
 * MIT License
 * Copyright (c) 2024 Markus Iser 
 */

#ifndef SRC_TRANSFORM_NORMALIZE_H_
#define SRC_TRANSFORM_NORMALIZE_H_

#include <vector>
#include <algorithm>

#include "src/util/StreamBuffer.h"


void determine_counts(const char* filename, int& nvars, int& nclauses) {
    StreamBuffer in(filename);
    nvars = 0; 
    nclauses = 0;
    while (in.skipWhitespace()) {
        if (*in == 'c' || *in == 'p') {
            if (!in.skipLine()) break;
        } else {
            int plit;
            while (in.readInteger(&plit)) {
                if (plit == 0) break;
                nvars = std::max(abs(plit), nvars);
            }
            nclauses++;
        }
    }
}

/**
 * @brief Normalizes a CNF formula by removing comments and 
 * generating a header based on the real number of clauses
 * and the maximum variable index.
 * 
 * @param filename 
 */
void normalize(const char* filename) {
    StreamBuffer in(filename);
    int vars, clauses;
    determine_counts(filename, vars, clauses);
    std::cout << "p cnf " << vars << " " << clauses << std::endl;
    while (in.skipWhitespace()) {
        if (*in == 'c' || *in == 'p') {
            if (!in.skipLine()) break;
        } else {
            int plit;
            while (in.readInteger(&plit)) {
                if (plit == 0) break;
                std::cout << plit << " ";
            }
            std::cout << "0" << std::endl;
        }
    }
}

/**
 * @brief Sanitizes a CNF formula by removing comments and generating a normalized header.
 * Removes duplicate literals from clauses and removes tautological clauses while preserving
 * the order of clauses and literals.
 * 
 * @param filename
 */
void sanitize(const char* filename) {
    StreamBuffer in(filename);

    int vars, clauses;
    determine_counts(filename, vars, clauses);
    std::cout << "p cnf " << vars << " " << clauses << std::endl;

    // set mask[lit] to clause number if lit is present in clause
    unsigned* mask = (unsigned*)calloc(2*vars + 2, sizeof(unsigned));
    mask += vars + 1;

    std::vector<int> clause;
    unsigned stamp = 0;
    while (in.skipWhitespace()) {
        if (*in == 'c' || *in == 'p') {
            if (!in.skipLine()) break;
        } else {
            ++stamp;
            bool tautological = false;
            int plit;
            while (in.readInteger(&plit)) {
                if (plit == 0) break;
                if (mask[-plit] == stamp) {
                    tautological = true;
                    break;
                }
                else if (mask[plit] != stamp) {
                    mask[plit] = stamp;
                    clause.push_back(plit);
                }
            }
            if (!tautological) {
                for (int plit : clause) {
                    std::cout << plit << " ";
                }
                std::cout << "0" << std::endl;
            }
            else {
                in.skipLine();
            }
            clause.clear();
        }
    }
}


/**
 * @brief Checks if a CNF formula is already sanitized.
 * 
 * @param filename 
 * @return true if neither duplicate literals nor tautological clauses are present
 * @return false otherwise
 */
bool check_sanitized(const char* filename) {
    StreamBuffer in(filename);

    int vars, clauses;
    determine_counts(filename, vars, clauses);
    
    // set mask[lit] to clause number if lit is present in clause
    int* mask = (int*)calloc(2*vars + 2, sizeof(int));
    mask += vars + 1;

    unsigned stamp = 0;
    while (in.skipWhitespace()) {
        if (*in == 'c' || *in == 'p') {
            if (!in.skipLine()) break;
        } else {
            ++stamp;
            int plit;
            while (in.readInteger(&plit)) {
                if (plit == 0) break;
                if (static_cast<unsigned>(mask[plit]) == stamp || static_cast<unsigned>(mask[-plit]) == stamp) {
                    return false;
                }
            }
        }
    }
    return true;
}



#endif  // SRC_TRANSFORM_NORMALIZE_H_
