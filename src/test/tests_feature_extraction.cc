#include <iostream>
#include <array>
#include <cstdio>
#include <unordered_map>
#include <filesystem>
#include <string>

#include "src/test/Util.h"
#include "src/extract/CNFBaseFeatures.h"
#include "src/extract/OPBBaseFeatures.h"
#include "src/extract/WCNFBaseFeatures.h"
#include "src/extract/CNFGateFeatures.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

template <typename Extractor>
void extract(const char *test_file, const char *expected_record_file)
{
    auto expected_record = record_to_map<double>(expected_record_file);
    Extractor stats(test_file);
    stats.extract();
    auto record = stats.getFeatures();
    auto names = stats.getNames();
    CHECK(record.size() == expected_record.size());
    for (unsigned i = 0; i < record.size(); i++)
    {
        CHECK_MESSAGE(fequal(expected_record[names[i]], record[i]), ("\nUnexpected record for feature '" + names[i] + "'\nExpected: " + std::to_string(expected_record[names[i]]) + "\nActual: " + std::to_string(record[i])));
    }
}

const std::string test_dir = "src/test/resources/test_files/";
const std::string records_dir = "src/test/resources/expected_records/";

TEST_CASE("Feature extraction")
{
    SUBCASE("CNF base")
    {
        const auto test_file = test_dir + "cnf_test.cnf.xz";
        const auto expected_record_file = records_dir + "cnf_base.txt";
        extract<CNF::BaseFeatures<>>(test_file.c_str(), expected_record_file.c_str());
    }

    SUBCASE("CNF gates")
    {
        const auto test_file = test_dir + "cnf_test.cnf.xz";
        const auto expected_record_file = records_dir + "cnf_gates.txt";
        extract<CNFGateFeatures<>>(test_file.c_str(), expected_record_file.c_str());
    }

    SUBCASE("WCNF base")
    {
        const auto test_file = test_dir + "wcnf_test.wcnf.xz";
        const auto expected_record_file = records_dir + "wcnf_base.txt";
        extract<WCNF::BaseFeatures<>>(test_file.c_str(), expected_record_file.c_str());
    }

    SUBCASE("OPB base")
    {
        const auto test_file = test_dir + "opb_test.opb.xz";
        const auto expected_record_file = records_dir + "opb_base.txt";
        extract<OPB::BaseFeatures<>>(test_file.c_str(), expected_record_file.c_str());
    }
}
