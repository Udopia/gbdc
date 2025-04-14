/**
 * MIT License
 * Copyright (c) 2024 Markus Iser
 */

#include "cnf2cnf.h"
#include "src/extract/CNFSaniCheck.h"

/**
 * @brief Normalises a CNF Formula
 * - Removes comments and generates a normalised header
 * - Replaces sequences of whitespace with a single space
 * - Prints one clause per line
 */
void CNF::Normaliser::run() {
    StreamBuffer in(filename_);

    CNF::SaniCheck ana(filename_, false);
    ana.run();
    std::cout << "p cnf " << (unsigned)ana.getFeature("norm_vars") << " "
              << (unsigned)ana.getFeature("norm_clauses") << std::endl;

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
 * @brief Sanitises a CNF formula
 * - Removes duplicate literals from clauses while preserving literal order
 * - Removes tautological clauses
 * - Outputs the normalised formula (cf. CNF::Normaliser)
 */
void CNF::Sanitiser::run() {
    StreamBuffer in(filename_);

    CNF::SaniCheck ana(filename_, true);
    ana.run();
    std::cout << "p cnf " << (unsigned)ana.getFeature("sani_vars") << " "
              << (unsigned)ana.getFeature("sani_clauses") << std::endl;

    // set mask[lit] to clause number if lit is present in clause
    unsigned *mask = (unsigned *)calloc(
        2 * (unsigned)ana.getFeature("norm_vars") + 2, sizeof(unsigned));
    mask += (unsigned)ana.getFeature("norm_vars") + 1;

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
                } else if (mask[plit] != stamp) {
                    mask[plit] = stamp;
                    clause.push_back(plit);
                }
            }
            if (!tautological) {
                for (int plit : clause) {
                    std::cout << plit << " ";
                }
                std::cout << "0" << std::endl;
            } else {
                in.skipLine();
            }
            clause.clear();
        }
    }
}