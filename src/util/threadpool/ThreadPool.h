#pragma once

#include <string>
#include <tuple>
#include <vector>
#include <cstdint>
#include <mutex>
#include <atomic>
#include <barrier>
#include "MPSCQueue.h"
#include <unordered_map>
#include "Util.h"
#include <tuple>
#include <variant>
#include "MemoryManagement.h"

static void debug_msg(const std::string &msg);

namespace threadpool
{
    template <typename... Args>
    struct job_t
    {
        std::tuple<Args...> args;
        size_t termination_count = 0;
        size_t emn = 0; // estimated memory needed
        size_t memnbt;  // memory needed before termination
        size_t idx;
        job_t() {}
        job_t(std::tuple<Args...> _args, size_t buffer_per_job = 0, size_t _idx = 0) : args(_args), memnbt(buffer_per_job), idx(_idx) {}
        void terminate_job(size_t _memnbt)
        {
            memnbt = _memnbt;
            ++termination_count;
        }
    };

    using compute_ret_t = std::vector<std::tuple<std::string, std::string, std::variant<double, long, std::string>>>;
    using compute_arg_t = std::tuple<std::string, std::string, std::unordered_map<std::string, long>>;
    using compute_t = std::tuple<compute_ret_t, compute_arg_t>;
    using extract_ret_t = std::unordered_map<std::string, std::variant<double, std::string>>;
    using extract_arg_t = std::string;

    // result consists of result of computation, bool indicating whether computation was successful and string to file
    template <typename T>
    using result_t = std::tuple<T, bool, std::string>;

    template <typename Ret, typename... Args>
    class ThreadPool
    {
    private:
        MPSCQueue<job_t<Args...>> jobs;
        std::vector<std::thread> threads;
        std::unordered_map<size_t, job_t<Args...>> in_process;
        std::atomic<size_t> termination_counter = 0;
        std::atomic<size_t> ready = 0;
        std::mutex jobs_m;
        std::shared_ptr<MPSCQueue<result_t<Ret>>> results;
        std::function<Ret(Args...)> func;

        /**
         * @brief Work routine for worker threads.
         */
        void work()
        {
            bool terminated = false;
            init_thread_datum();
            ready.fetch_add(1);
            while (ready.load() != (threads.size() + 1U))
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            while (next_job())
            {
                if (terminated)
                {
                    cleanup_termination();
                    terminated = false;
                }
                try
                {
                    wait_for_starting_permission();
                    const auto result = std::apply(func, in_process[tl_id].args);
                    output_result(result, true);
                    finish_job();
                }
                catch (const TerminationRequest &tr)
                {
                    if (tr.what())
                        debug_msg(tr.what());
                    requeue_job(tr.memnbt);
                    terminated = true;
                }
            }
            if (terminated)
            {
                cleanup_termination();
            }
            finish_work();
        }

        bool next_job()
        {
            std::unique_lock<std::mutex> l(jobs_m);
            if (!jobs.empty())
            {
                in_process[tl_id] = jobs.pop();
                return true;
            }
            return false;
        }

        /**
         * @brief Resets the state of the thread such that a new job can be started.
         */
        void termination_penalty()
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(50) * termination_counter.load(relaxed));
        }

        /**
         * @brief Called after processing all jobs. Makes threads joinable.
         */
        void finish_work()
        {
            results->quit();
        }
        /**
         * @brief Subroutine for wait_for_starting_permission()
         *
         * @param size Size of required memory.
         * @return True if memory has been reserved and the job can be started.
         * False if thread has to wait until memory is freed up.
         */
        bool can_start(size_t size)
        {
            if (can_alloc(size))
            {
                tl_data->set_reserved(size);
                return true;
            }
            return false;
        }
        /**
         * @brief Determines whether the current job of a thread can be started,
         * based on currently available memory and the (estimated) memory to be required by the job.
         *
         */
        void wait_for_starting_permission()
        {
            debug_msg("Waiting...\n");
            auto mem_needed = in_process[tl_id].memnbt;
            if (mem_needed > mem_max){
                terminate(mem_needed);
            }
            while (!can_start(mem_needed))
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(0));
            }
            threads_working.fetch_add(1, relaxed);
            debug_msg("Start job with index " + std::to_string(in_process[tl_id].idx) + " ...");
        }
        /**
         * @brief Called after forceful termination of a job. It stores the peak memory allocated
         * during the execution of the job. This information can later be used to determine
         * more accurately, whether a job can be fit into memory or not. Finally, it locks the job queue
         * and places job at the end of the queue.
         *
         * @param memnbt Memory needed by the job before its termination.
         */
        void requeue_job(size_t memnbt)
        {
            in_process[tl_id].terminate_job(memnbt);
            if (can_fit(in_process[tl_id].memnbt))
            {
                debug_msg("Cannot requeue job!\nMemory needed before termination: " + std::to_string(in_process[tl_id].memnbt) + "\nMaximum amount of memory available: " + std::to_string(mem_max - (size_t)1e6));
                output_result({}, false);
                return;
            }
            debug_msg("Requeuing job!\nMemory needed before termination: " + std::to_string(in_process[tl_id].memnbt) + "\nMaximum amount of memory available: " + std::to_string(mem_max - (size_t)1e6));
            std::unique_lock<std::mutex> lock(jobs_m);
            jobs.push(in_process[tl_id]);
        }
        /**
         * @brief General cleanup function. Resets the state of the thread such that a new job can be started.
         */
        void cleanup_termination()
        {
            finish_job();
            termination_counter.fetch_add(1, relaxed);
            termination_ongoing.unlock();
            termination_penalty();
        }

        void finish_job()
        {
            threads_working.fetch_sub(1, relaxed);
            size_t to_be_unreserved = tl_data->unreserve_memory();
            unreserve_memory(to_be_unreserved);
            tl_data->reset();
        }

        void output_result(const Ret &result, const bool success)
        {
            debug_msg("Extraction" + std::string(success ? " " : " not ") + "succesful!\n");
            results->push({result, success, std::get<std::string>(in_process[tl_id].args)});
        }

        void init_jobs(const std::vector<std::tuple<Args...>> &args)
        {
            size_t idx = 0;
            for (const auto &arg : args)
            {
                jobs.push(job_t<Args...>(arg, buffer_per_job, idx++));
            }
        }

        /**
         * @brief Initializes worker threads.
         *
         * @param num_threads Number of worker threads
         */
        void init_threads(std::uint32_t num_threads)
        {
            for (std::uint32_t i = 0; i < num_threads; ++i)
            {
                threads.emplace_back(std::thread(&ThreadPool::work, this));
            }
        }

    public:
        /**
         * @brief Called by master thread at the beginning of the program. Frequently samples execution data.
         */
        void start_threadpool()
        {
            ready.fetch_add(1);
            while (ready.load() != (threads.size() + 1U))
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            csv_t csv("data.csv", {"time", "allocated", "reserved", "jobs"});
            size_t tp, total_allocated;
            size_t begin = std::chrono::steady_clock::now().time_since_epoch().count();
            while (!results->done())
            {
                total_allocated = 0;
                std::for_each(thread_data.begin(), thread_data.end(), [&](const std::shared_ptr<thread_data_t> &td)
                              { total_allocated += td->mem_allocated; });
                tp = std::chrono::steady_clock::now().time_since_epoch().count();
                csv.write_to_file({tp - begin, total_allocated, reserved.load(relaxed), threads_working.load(relaxed)});
                std::this_thread::sleep_for(std::chrono::milliseconds(20));
            }
            for (uint i = 0; i < threads.size(); ++i)
                threads[i].join();
            debug_msg("Number of references to queue left: " + std::to_string(results.use_count()));
        }

        auto get_result_queue()
        {
            return results;
        }

        /**
         * @brief Construct a new ThreadPool object. Initializes the job queue and its worker threads.
         *
         * @param paths List of paths to instances, for which feature extraction shall be done.
         * @param _mem_max Maximum amount of memory available to the ThreadPool.
         * @param _jobs_max Maximum amount of jobs that can be executed in parallel, i.e. number of threads.
         */
        ThreadPool(std::uint64_t _mem_max, std::uint32_t _jobs_max, std::function<Ret(Args...)> _func, std::vector<std::tuple<Args...>> args) : jobs(_jobs_max), in_process(_jobs_max), func(_func)
        {
            results = std::make_shared<MPSCQueue<result_t<Ret>>>(_jobs_max);
            init_memory_management(_mem_max);
            init_jobs(args);
            init_threads(_jobs_max);
        }
    };
};

static void debug_msg(const std::string &msg)
{
#ifndef NDEBUG
    std::cerr << "[Thread " + std::to_string(tl_id) + "] " + msg << std::endl;
#endif
};