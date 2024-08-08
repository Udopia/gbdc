#pragma once
#include <cstdint>

#include <cstdint>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <locale.h>
#include <dlfcn.h>
#include <thread>
#include <cassert>
#include <functional>
#include <atomic>
#include <bit>
#include <mutex>
#include "../ResourceLimits.h"

class TerminationRequest : public std::runtime_error
{
public:
    size_t memnbt;
    explicit TerminationRequest(size_t _memnbt) : std::runtime_error(""), memnbt(_memnbt) {}
};

struct thread_data_t
{
    /* currently allocated memory */
    size_t mem_allocated = 0UL;
    /* currently reserved memory */
    size_t mem_reserved = 0UL;
    /* peak allocated memory for current job */
    size_t peak_mem_allocated = 0UL;
    /* number of calls to malloc for current job */
    size_t num_allocs = 0UL;
    /* runtime */
    size_t runtime_begin = 0UL;
    bool memout = false;

    thread_data_t() {}

    void begin_time()
    {
        runtime_begin = std::chrono::steady_clock::now().time_since_epoch().count();
    }

    void flag_for_memout(){
        memout = true;
    }

    void inc_allocated(const size_t size)
    {
        ++num_allocs;
        mem_allocated += size;
        peak_mem_allocated = std::max(mem_allocated, peak_mem_allocated);
    }

    void dec_allocated(const size_t size)
    {
        mem_allocated -= size;
    }

    void set_reserved(const size_t size)
    {
        mem_reserved = size;
    }

    size_t unreserve_memory()
    {
        size_t tmp = mem_reserved;
        set_reserved(0UL);
        return tmp < mem_allocated ? 0UL : tmp - mem_allocated;
    }

    void reset()
    {
        memout = false;
        peak_mem_allocated = 0UL;
        num_allocs = 0UL;
    }

    size_t rmem_needed(const size_t size) const
    {
        // if reserved memory is not sufficient to harbor allocation size,
        // additional reserved memory is needed
        if (mem_reserved > mem_allocated)
        {
            return (mem_reserved - mem_allocated) > size ? 0UL : size - (mem_reserved - mem_allocated);
        }
        return size;
    }

    size_t rmem_not_needed(const size_t size) const
    {
        // can return allocated memory down to initially reserved amount
        return mem_allocated > mem_reserved ? std::min(mem_allocated - mem_reserved, size) : 0UL;
    }
};

constexpr std::uint64_t buffer_per_job = static_cast<std::uint64_t>(2e7); // 20MB buffer to accommodate StreamBuffers
constexpr std::uint64_t UNTRACKED = UINT8_MAX;
constexpr std::memory_order relaxed = std::memory_order_relaxed;

thread_local std::uint64_t tl_id = UNTRACKED;
thread_local std::shared_ptr<thread_data_t> tl_data;

std::mutex termination_ongoing;
std::mutex m_thread_data;

size_t mem_max = UINT64_MAX;
size_t runtime_max = UINT64_MAX;

std::atomic<std::uint32_t> thread_id_counter = 0;
std::atomic<size_t> peak = 0, reserved = 0;

std::vector<std::shared_ptr<thread_data_t>> thread_data;

void reset_memory_management()
{
    peak.store(0);
    reserved.store(0);
    thread_id_counter.store(0);
    thread_data.clear();
    thread_data.shrink_to_fit();
}

void init_memory_management(size_t _mem_max, size_t _runtime_max)
{
    mem_max = _mem_max * static_cast<size_t>(1e6);         // bytes to megabytes
    runtime_max = _runtime_max * static_cast<size_t>(1e9); // seconds to nanoseconds
}

void init_thread_datum()
{
    tl_data = std::make_shared<thread_data_t>();
    std::unique_lock<std::mutex> l(m_thread_data);
    tl_id = thread_id_counter.fetch_add(1);
    thread_data.push_back(tl_data);
}

bool has_time()
{
    size_t current_time_ns = std::chrono::steady_clock::now().time_since_epoch().count();
    return (current_time_ns - tl_data->runtime_begin) < runtime_max;
}

bool can_fit(size_t size)
{
    return size <= mem_max;
}

/**
 * @brief Gives back reserved memory.
 *
 * @param size Amount of memory to be given back.
 */
void unreserve_memory(size_t size)
{
    assert((reserved.load(relaxed) < mem_max));
    reserved.fetch_sub(size, relaxed);
    assert((reserved.load(relaxed) < mem_max));
}

void inc_allocated(size_t size)
{
    tl_data->inc_allocated(size);
}

void dec_allocated(size_t size)
{
    unreserve_memory(tl_data->rmem_not_needed(size));
    tl_data->dec_allocated(size);
}

/**
 * @brief Determines whether the requested amount of memory can be allocated or not,
 * based on current memory consumption.
 *
 * @param size Amount of memory requested
 * @return True if the requested amount of memory has been reserved and can subsequently be allocated.
 * False if the requested amount of memory is too large to fit into the currently available.
 */
bool can_alloc(size_t size)
{
    if (size == 0)
        return true;
    bool exchanged = false;
    auto _reserved = reserved.load(relaxed);
    while ((_reserved + size) <= mem_max &&
           !(exchanged = reserved.compare_exchange_weak(_reserved, _reserved + size, relaxed, relaxed)))
        ;
    return exchanged;
}

void terminate(size_t lrq)
{
    throw TerminationRequest(std::max(tl_data->mem_allocated + lrq, std::max(tl_data->peak_mem_allocated, tl_data->mem_reserved)));
}
