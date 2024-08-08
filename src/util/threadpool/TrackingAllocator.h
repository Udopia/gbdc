#pragma once
#include "MemoryManagement.h"

template <typename T>
class TrackingAllocator : public std::allocator<T>
{
public:
    using Base = std::allocator<T>;
    using value_type = T;

    TrackingAllocator() : Base()
    {
        if (tl_id == UNTRACKED)
        {
            init_thread_datum();
        }
    }

    T *allocate(size_t n)
    {
        if (!has_time())
        {
            throw TimeLimitExceeded();
        }
        size_t bytes = n * sizeof(T);
        if (tl_data->memout || !can_fit(std::max(tl_data->mem_reserved, tl_data->mem_allocated + bytes)))
        {
            throw MemoryLimitExceeded();
        }
        size_t rmem_needed = tl_data->rmem_needed(bytes);
        while (!can_alloc(rmem_needed))
        {
            if (termination_ongoing.try_lock())
                terminate(bytes);
        }
        inc_allocated(bytes);
        return static_cast<T *>(::operator new(bytes));
    }

    void deallocate(T *p, size_t n) noexcept
    {
        std::size_t bytes = n * sizeof(T);
        dec_allocated(bytes);
        delete (p);
    }

    template <typename U>
    struct rebind
    {
        using other = TrackingAllocator<U>;
    };

private:
};
