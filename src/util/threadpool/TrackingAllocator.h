#pragma once
#include "MemoryManagement.h"

template <typename T>
class TrackingAllocator : public std::allocator<T>
{
public:
    using Base = std::allocator<T>; // Alias for the base class
    using value_type = T;

    // Inherit constructors
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
        // std::cout << "Allocating " << bytes << " bytes. Total allocated: " << tl_data->mem_allocated << " bytes.\n";
        return static_cast<T *>(::operator new(bytes));
    }

    void deallocate(T *p, std::size_t n) noexcept
    {
        std::size_t bytes = n * sizeof(T);
        dec_allocated(tl_id, bytes);
        // std::cout << "Deallocating " << bytes << " bytes. Total outstanding: " << tl_data->mem_allocated << " bytes.\n";
        delete (p);
    }

    template <typename U>
    struct rebind
    {
        using other = TrackingAllocator<U>;
    };

private:
    static inline std::size_t totalAllocated = 0;
    static inline std::size_t totalOutstanding = 0;
};
