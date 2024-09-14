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

using result_t = threadpool::result_t<threadpool::extract_ret_t>;
std::string desc = CNF::BaseFeatures<>("").getRuntimeDesc();

auto test_extract(std::string filepath)
{
    tp::extract_ret_t m;
    CNF::BaseFeatures<TrackingAllocator> stats(filepath.c_str());
    try
    {
        stats.extract();
        auto names = stats.getNames();
        auto features = stats.getFeatures();
        for (size_t i = 0; i < names.size(); ++i)
        {
            m[names[i]] = features[i];
        }
    }
    catch (MemoryLimitExceeded e)
    {
        m[desc] = "memout";
    }
    catch (TimeLimitExceeded e)
    {
        m[desc] = "timeout";
    }
    m["local"] = filepath;
    return m;
};

auto create_threadpool(std::uint64_t _mem_max, std::uint32_t _jobs_max, std::string folderPath, std::uint64_t _runtime_max = 100ULL)
{
    namespace fs = std::filesystem;
    std::vector<std::tuple<std::string>> paths;
    for (const auto &entry : fs::directory_iterator(folderPath))
    {
        if (has_extension(entry, "cnf"))
        {
            paths.push_back(entry.path());
        }
    }
    auto tp = std::make_shared<tp::ThreadPool<tp::extract_ret_t, std::string>>(_mem_max, _jobs_max, _runtime_max);
    for (const auto &path : paths)
    {
        tp->push_job(test_extract, path);
    }
    return std::make_pair(tp, paths.size());
}

template <typename TpPointer, typename Func>
void run_threadpool(TpPointer &&tp, Func &&f, size_t num_jobs, size_t timeout = UINT64_MAX)
{
    auto begin = std::chrono::steady_clock::now().time_since_epoch().count();
    while (num_jobs != 0)
    {
        if (tp->result_ready())
        {
            auto r = tp->pop_result();
            f(r);
            --num_jobs;
        }
        else
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        auto end = std::chrono::steady_clock::now().time_since_epoch().count();
        if ((end - begin) > timeout)
        {
            throw;
        }
    }
}

TEST_CASE("Insufficient max memory")
{
    size_t timeout = 15UL * 1e9; // 15 seconds timeout
    const std::string folderPath = "src/test/resources/test_files";
    auto [tp, num_jobs] = create_threadpool(1, 4, folderPath);
    auto func = [&](result_t &r)
    {
        CHECK_EQ(std::get<std::string>(std::get<tp::RETURN_VALUE>(r)[desc]), "memout");
    };
    run_threadpool(tp, func, num_jobs, timeout);
}

std::unordered_map<std::string, int> status;
std::unordered_map<std::string, tp::extract_ret_t> first_run;

TEST_CASE("Threadpool_Extract")
{
    SUBCASE("Basefeature extraction")
    {
        const std::string folderPath = "src/test/resources/test_files";
        auto [tp, num_jobs] = create_threadpool(45, 3, folderPath);
        auto func = [&](result_t &r)
        {
            std::string path = std::get<std::string>(std::get<tp::RETURN_VALUE>(r)["local"]);
            status[path] = std::get<tp::STATUS>(r);
            CHECK_EQ(std::get<tp::STATUS>(r) == tp::SUCCESS, !std::get<tp::RETURN_VALUE>(r).empty());
            first_run[path] = std::get<tp::RETURN_VALUE>(r);
        };
        run_threadpool(tp, func, num_jobs);
    }

    SUBCASE("Assert equal behavior of runs with the same set of instances")
    {
        const std::string folderPath = "src/test/resources/test_files";
        auto [tp, num_jobs] = create_threadpool(45, 5, folderPath);
        auto func = [&](result_t &r)
        {
            std::string path = std::get<std::string>(std::get<tp::RETURN_VALUE>(r)["local"]);
            // check if extraction succeeded/failed in both runs
            CHECK_EQ(std::get<tp::STATUS>(r), status[path]);
            // check if feature record is the same
            check_eqset(std::get<tp::RETURN_VALUE>(r), first_run[path]);
        };
        run_threadpool(tp, func, num_jobs);
    }

    SUBCASE("Assert correct TIMEOUT behavior")
    {
        const std::string folderPath = "src/test/resources/test_files";
        auto [tp, num_jobs] = create_threadpool(45, 5, folderPath, 0);
        auto func = [&](result_t &r)
        {
            CHECK_EQ(std::get<std::string>(std::get<tp::RETURN_VALUE>(r)[desc]), "timeout");
        };
        run_threadpool(tp, func, num_jobs);
    }
}
