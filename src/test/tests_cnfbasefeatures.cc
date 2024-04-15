#include <iostream>
#include <array>
#include <cstdio>
#include <unordered_map>
#include <filesystem>
#include <string>

#include "src/test/Util.h"
#include "src/extract/CNFBaseFeatures.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

TEST_CASE("CNFBaseFeatures")
{
    SUBCASE("Basefeature extraction")
    {
        const char *cnf_file = "src/test/resources/01bd0865ab694bc71d80b7d285d5777d-shuffling-2-s1480152728-of-bench-sat04-434.used-as.sat04-711.cnf.xz";
        const char *expected_record_file = "src/test/resources/expected_record.txt";
        auto expected_record = record_to_map<double>(expected_record_file);
        CNF::BaseFeatures stats(cnf_file);
        stats.extract();
        std::vector<double> record = stats.getFeatures();
        std::vector<std::string> names = stats.getNames();
        CHECK(record.size() == expected_record.size());
        for (unsigned i = 0; i < record.size(); i++)
        {
            CHECK_MESSAGE(fequal(expected_record[names[i]], record[i]), ("\nUnexpected record for feature '" + names[i] + "'\nExpected: " + std::to_string(expected_record[names[i]]) + "\nActual: " + std::to_string(record[i])));
        }
    }
}
