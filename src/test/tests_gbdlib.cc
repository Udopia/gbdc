#include <iostream>
#include <array>
#include <cstdio>
#include <unordered_map>
#include <filesystem>
#include <string>

#include "src/test/Util.h"
#include "src/extract/CNFBaseFeatures.h"
#include "src/util/threadpool/TrackingAllocator.h"
#include "src/gbdlib.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

namespace tp = threadpool;
namespace fs = std::filesystem;

TEST_CASE("GBDLib")
{
    SUBCASE("BaseFeature_Names")
    {
        auto super = feature_names<CNF::BaseFeatures<>>();
        auto sub = CNF::BaseFeatures("").getNames();
        check_subset(sub, super);
        CHECK_EQ(sub.size(), super.size() - 1);
        CHECK_EQ(super.back(), get_runtime_desc<0>());
    }

    SUBCASE("BaseFeatures_ThreadPool_Extract_Correct")
    {
        std::string cnf_file = "src/test/resources/test_files/cnf_test.cnf.xz";
        const char *expected_record_file = "src/test/resources/expected_records/cnf_base.txt";
        auto map = record_to_map<std::variant<double, std::string>>(expected_record_file);

        auto dict = tp_extract_features<CNF::BaseFeatures<TrackingAllocator>>(cnf_file);
        check_subset(map, dict);
    }

    SUBCASE("BaseFeatures_ThreadPool_Async")
    {
        const std::string folderPath = "src/test/resources/test_files";
        std::vector<std::tuple<std::string>> paths;
        for (const auto &entry : fs::directory_iterator(folderPath))
        {
            if (has_extension(entry,"cnf"))
            {
                paths.push_back(entry.path());
            }
        }
        auto q = tp_extract<CNF::BaseFeatures>(1UL << 26UL, 4U, paths);
        std::cerr << q.use_count() << "\n";
        size_t job_counter = 0;
        while (!q->done())
        {
            if (!q->empty())
            {
                ++job_counter;
                auto f = q->pop();
                CHECK_EQ(std::get<1>(f), !std::get<0>(f).empty());
            }
            else
                std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
        CHECK_EQ(paths.size(), job_counter);
    }

    SUBCASE("Version")
    {
        fs::current_path("..");
        CHECK_EQ(version(), "0.2.45");
    }

    // SUBCASE("BaseFeatures_ProcessPool_Correct")
    // {
    //     std::string cnf_file = "src/test/resources/01bd0865ab694bc71d80b7d285d5777d-shuffling-2-s1480152728-of-bench-sat04-434.used-as.sat04-711.cnf.xz";
    //     const char *expected_record_file = "src/test/resources/expected_record.txt";
    //     auto map = record_to_map<std::variant<double, std::string>>(expected_record_file);
    //     auto dict = extract_features<CNF::BaseFeatures>(cnf_file, 1UL << 36UL, 1UL << 36UL);
    //     check_subset(map, dict);
    //     CHECK(dict.find("base_features_runtime") != dict.end());
    // }
}
