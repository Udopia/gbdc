#pragma once

#include <string>
#include <tuple>
#include <vector>
#include <cstdint>
#include <mutex>
#include <atomic>
#include "LockedQueue.h"
#include <unordered_map>
#include <tuple>
#include <variant>
#include <functional>
#include "MemoryManagement.h"

static void debug_msg(const std::string &msg);

namespace threadpool
{
    template <typename Ret, typename... Args>
    struct job_t
    {
        std::function<Ret(Args...)> func;
        std::tuple<Args...> args;
        size_t termination_count = 0;
        size_t emn = 0; // estimated memory needed
        size_t memnbt;  // memory needed before termination

        job_t() {}
        job_t(std::function<Ret(Args...)> _func, std::tuple<Args...> _args, size_t buffer_per_job) : func(_func), args(_args), memnbt(buffer_per_job) {}
        job_t(const job_t &other)
        {
            args = other.args;
            func = other.func;
            termination_count = other.termination_count;
            memnbt = other.memnbt;
        }

        job_t &operator=(const job_t &other)
        {
            args = other.args;
            func = other.func;
            termination_count = other.termination_count;
            memnbt = other.memnbt;
            return *this;
        }

        Ret do_job()
        {
            return std::apply(func, args);
        }

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
    using extract_arg_t = std::tuple<std::string, std::string>;

    enum R_STATUS
    {
        SUCCESS,
        MEMOUT,
        TIMEOUT,
        FILEOUT
    };

    enum R_INDICES
    {
        RETURN_VALUE,
        STATUS
    };

    template <typename T>
    using result_t = std::tuple<T, int>;

    template <typename Ret, typename... Args>
    class ThreadPool
    {
    private:
        LockedQueue<job_t<Ret, Args...>> jobs;
        LockedQueue<result_t<Ret>> results;
        std::vector<std::thread> threads;
        std::vector<job_t<Ret, Args...>> in_process;
        size_t wait_ms = 5;
        bool done = false;

        /**
         * @brief Work routine for worker threads.
         */
        void work()
        {
            init_thread_datum();
            while (!done)
            {
                if (!try_get_job())
                {
                    std::this_thread::sleep_for(std::chrono::milliseconds(wait_ms));
                    continue;
                }
                try
                {
                    wait_for_starting_permission();
                    start_job();
                    const auto result = in_process[tl_id].do_job();
                    push_result(result, SUCCESS);
                    finish_job();
                }
                catch (const TerminationRequest &tr)
                {
                    in_process[tl_id].terminate_job(tr.memnbt);
                    if (!try_requeue_job())
                        push_result({}, MEMOUT);
                    finish_job();
                    termination_ongoing.unlock();
                }
                catch (const MemoryLimitExceeded &tl)
                {
                    push_result({}, MEMOUT);
                    finish_job();
                }
                catch (const TimeLimitExceeded &tl)
                {
                    push_result({}, TIMEOUT);
                    finish_job();
                }
            }
        }

        bool try_get_job()
        {
            return jobs.try_pop(in_process[tl_id]);
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
         */
        void wait_for_starting_permission()
        {
            debug_msg("Waiting...\n");
            auto mem_needed = in_process[tl_id].memnbt;
            if (!can_fit(mem_needed))
            {
                tl_data->flag_for_memout();
                return;
            }
            while (!can_start(mem_needed))
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(wait_ms));
            }
        }

        /**
         * @brief Requeues job if memory requirements allow it.
         *
         * @return True if the job has been requeued, false if its memory requirements are too high.
         */
        bool try_requeue_job()
        {
            size_t mem_available = mem_max > buffer_per_job ? mem_max - buffer_per_job : 0;
            if (!can_fit(in_process[tl_id].memnbt))
            {
                debug_msg("Cannot requeue job!\nMemory needed before termination: " + std::to_string(in_process[tl_id].memnbt) + "\nMaximum amount of memory available: " + std::to_string(mem_available));
                return false;
            }
            debug_msg("Requeuing job!\nMemory needed before termination: " + std::to_string(in_process[tl_id].memnbt) + "\nMaximum amount of memory available: " + std::to_string(mem_available));
            jobs.push(in_process[tl_id]);
            return true;
        }

        void start_job()
        {
            debug_msg("Start job ...");
            tl_data->begin_time();
        }

        void finish_job()
        {
            size_t to_be_unreserved = tl_data->unreserve_memory();
            unreserve_memory(to_be_unreserved);
            tl_data->reset();
        }

        void push_result(const Ret &result, const int status)
        {
            debug_msg("Extraction" + std::string(status == SUCCESS ? " " : " not ") + "succesful!\n");
            results.push({result, status});
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
         * @brief Construct a new ThreadPool object. Initializes the job queue and its worker threads.
         *
         * @param paths List of paths to instances, for which feature extraction shall be done.
         * @param _mem_max Maximum amount of memory available to the ThreadPool (in bytes).
         * @param _jobs_max Maximum amount of jobs that can be executed in parallel, i.e. number of threads.
         * @param _runtime_max Maximum amount of time a job may run (in seconds).
         */
        ThreadPool(std::uint64_t _mem_max, std::uint32_t _jobs_max, std::uint64_t _runtime_max) : jobs(), results(), in_process(_jobs_max)
        {
            init_memory_management(_mem_max, _runtime_max);
            init_threads(_jobs_max);
        }

        void stop()
        {
            done = true;
            for (uint i = 0; i < threads.size(); ++i)
                threads[i].join();
            reset_memory_management();
        }

        ~ThreadPool()
        {
            if (!done)
                stop();
        }

        void push_job(std::function<Ret(Args...)> func, std::tuple<Args...> args)
        {
            jobs.push(job_t<Ret, Args...>(func, args, buffer_per_job));
        }

        auto pop_result()
        {
            return results.pop();
        }

        bool result_ready()
        {
            return !results.empty();
        }
    };
};

static void debug_msg(const std::string &msg)
{
#ifndef NDEBUG
    std::cerr << "[Thread " + std::to_string(tl_id) + "] " + msg << std::endl;
#endif
};