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
#include "Util.h"

constexpr std::uint64_t buffer_per_job = 2e7; // 20MB buffer to accommodate StreamBuffers
constexpr std::uint64_t UNTRACKED = UINT8_MAX;
constexpr std::memory_order relaxed = std::memory_order_relaxed;
size_t mem_max = 1ULL << 30;

thread_local std::uint64_t tl_id = UNTRACKED;
std::atomic<std::uint32_t> thread_id_counter = 0;

std::mutex termination_ongoing;
std::atomic<uint16_t> threads_working = 0;
std::atomic<size_t> peak = 0, reserved = 0;

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

    thread_data_t() {}

    void inc_allocated(const size_t size)
    {
        ++num_allocs;
        mem_allocated += size;
        peak_mem_allocated = std::max(mem_allocated, peak_mem_allocated);
    }

    void dec_allocated(const size_t size)
    {
        mem_allocated -= size;
        if (mem_allocated > (1 << 30))
        {
            throw;
        }
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

std::mutex m_thread_data;
std::vector<std::shared_ptr<thread_data_t>> thread_data;
thread_local std::shared_ptr<thread_data_t> tl_data;

void init_memory_management(size_t _mem_max)
{
    mem_max = _mem_max;
}

void init_thread_datum()
{
    tl_data = std::make_shared<thread_data_t>();
    tl_id = thread_id_counter.fetch_add(1,relaxed);
    std::unique_lock<std::mutex> l(m_thread_data);
    thread_data.push_back(tl_data);
}

bool can_fit(size_t size)
{
    return size > mem_max - buffer_per_job;
}

/**
 * @brief Gives back reserved memory.
 *
 * @param size Amount of memory to be given back.
 */
void unreserve_memory(size_t size)
{
    assert((reserved.load(std::memory_order_acquire) >= size));
    assert((reserved.load(std::memory_order_acquire) <= mem_max));
    reserved.fetch_sub(size, std::memory_order_relaxed);
    assert((reserved.load(std::memory_order_acquire) <= mem_max));
}

void inc_allocated(size_t tid, size_t size)
{
    if (tid != tl_id)
    {
        return;
    }
    tl_data->inc_allocated(size);
}

void dec_allocated(size_t tid, size_t size)
{
    if (tid != tl_id)
    {
        return;
    }
    unreserve_memory(tl_data->rmem_not_needed(size));
    tl_data->dec_allocated(size);
}

size_t threshold()
{
    return 0UL;
    //    return buffer_per_job * threads_working.load(relaxed);
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
    if ((reserved.load(relaxed) > mem_max))
        throw;
    if (size == 0)
        return true;
    bool exchanged = false;
    auto _reserved = reserved.load(relaxed);
    while ((_reserved + size + threshold()) <= mem_max &&
           !(exchanged = reserved.compare_exchange_weak(_reserved, _reserved + size, relaxed, relaxed)))
        ;
    if ((reserved.load(relaxed) > mem_max))
        throw;

    return exchanged;
}

void terminate(size_t lrq)
{
    // exception allocation needed because of malloc call during libc exception allocation
    throw TerminationRequest(std::max(tl_data->mem_allocated + lrq, std::max(tl_data->peak_mem_allocated, tl_data->mem_reserved)));
}
