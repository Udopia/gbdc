#pragma once
#include "MemoryManagement.h"

template <typename T>
class TrackingAllocator : public std::allocator<T>
{
public:
    using Base = std::allocator<T>;
    using value_type = T;

    TrackingAllocator() : Base() {
        if (tl_id == UNTRACKED){
            init_thread_datum();
        }
    }

    T *allocate(std::size_t n)
    {
        std::size_t bytes = n * sizeof(T);
        while (!can_alloc(tl_data->rmem_needed(bytes)))
        {
            if (termination_ongoing.try_lock())
                terminate(bytes);
        }
        inc_allocated(tl_id, bytes);
        return static_cast<T *>(::operator new(bytes));
    }

    void deallocate(T *p, std::size_t n) noexcept
    {
        std::size_t bytes = n * sizeof(T);
        dec_allocated(tl_id, bytes);
        delete (p);
    }

    template <typename U>
    struct rebind
    {
        using other = TrackingAllocator<U>;
    };

private:
};
