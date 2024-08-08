#include <iostream>
#include <array>
#include <cstdio>
#include <unordered_map>
#include <filesystem>
#include <string>
#include <thread>

#include "src/test/Util.h"
#include "src/extract/CNFBaseFeatures.h"
#include "src/util/threadpool/ThreadPool.h"
#include "src/util/threadpool/TrackingAllocator.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

namespace tp = threadpool;

static std::vector<double> test_extract(std::string filepath)
{
    CNF::BaseFeatures<TrackingAllocator> stats(filepath.c_str());
    stats.extract();
    return stats.getFeatures();
};

TEST_CASE("Threadpool_Extract")
{
    SUBCASE("Basefeature extraction")
    {
        namespace fs = std::filesystem;
        const std::string folderPath = "src/test/resources/test_files";
        std::vector<std::tuple<std::string>> paths;
        for (const auto &entry : fs::directory_iterator(folderPath))
        {
            if (has_extension(entry, "cnf"))
            {
                paths.push_back(entry.path());
            }
        }
        tp::ThreadPool<std::vector<double>, std::string> tp(1UL << 26UL, 4U, test_extract, paths);
        auto q = tp.get_result_queue();
        std::thread tp_thread(&tp::ThreadPool<std::vector<double>, std::string>::start_threadpool, &tp);
        auto names = CNF::BaseFeatures("").getNames();
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
        tp_thread.join();
        CHECK_EQ(paths.size(), job_counter);
    }
}
