#include <iostream>
#include <array>
#include <cstdio>
#include <unordered_map>
#include <filesystem>
#include <string>

#include "src/test/Util.h"
#include "src/extract/CNFBaseFeatures.h"
#include "src/util/threadpool/TrackingAllocator.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

namespace tp = threadpool;
namespace fs = std::filesystem;

// TEST_CASE("GBDLib")
// {
//     SUBCASE("BaseFeature_Names")
//     {
//         auto ex = CNF::BaseFeatures("");
//         auto super = feature_names<CNF::BaseFeatures<>>();
//         auto sub = ex.getNames();
//         check_subset(sub, super);
//         CHECK_EQ(sub.size(), super.size() - 1);
//         CHECK_EQ(super.back(), ex.getRuntimeDesc());
//     }

//     SUBCASE("BaseFeatures_ThreadPool_Extract_Correct")
//     {
//         std::string cnf_file = "src/test/resources/test_files/cnf_test.cnf.xz";
//         const char *expected_record_file = "src/test/resources/expected_records/cnf_base.txt";
//         auto map = record_to_map<std::variant<double, std::string>>(expected_record_file);

//         auto dict = tp_extract_features<CNF::BaseFeatures<TrackingAllocator>>(cnf_file);
//         check_subset(map, dict);
//     }

//     SUBCASE("BaseFeatures_ThreadPool_Async")
//     {
//         const std::string folderPath = "src/test/resources/test_files";
//         std::vector<std::tuple<std::string>> paths;
//         for (const auto &entry : fs::directory_iterator(folderPath))
//         {
//             if (has_extension(entry,"cnf"))
//             {
//                 paths.push_back(entry.path());
//             }
//         }
//         auto q = tp_extract<CNF::BaseFeatures>(1UL << 26UL, 4U, paths);
//         std::cerr << q.use_count() << "\n";
//         size_t job_counter = 0;
//         while (!q->done())
//         {
//             if (!q->empty())
//             {
//                 ++job_counter;
//                 auto f = q->pop();
//                 // if successful, feature vector must not be empty
//                 CHECK_EQ(std::get<1>(f), !std::get<0>(f).empty());
//             }
//             else
//                 std::this_thread::sleep_for(std::chrono::milliseconds(500));
//         }
//         CHECK_EQ(paths.size(), job_counter);
//     }
// }
