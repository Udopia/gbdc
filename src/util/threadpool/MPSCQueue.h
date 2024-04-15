#pragma once
#include <queue>
#include <mutex>
#include <list>
#include <atomic>
#include <condition_variable>

template <typename T>
class MPSCQueue
{
private:
    std::queue<T> q;
    std::atomic<uint16_t> producers;
    std::mutex m;
    std::condition_variable cv;

public:
    explicit MPSCQueue(std::uint16_t _producers) : q(), producers(_producers), m() {}
    void push(T &&t)
    {
        std::unique_lock<std::mutex> l(m);
        q.push(t);
    }

    void push(const T &t)
    {
        {
            std::unique_lock<std::mutex> l(m);
            q.push(t);
        }
        cv.notify_one();
    }

    T pop()
    {
        std::unique_lock<std::mutex> l(m);
        if (q.empty())
            return T();
        T t = q.front();
        q.pop();
        return t;
    }

    bool empty()
    {
        return q.empty();
    }

    void quit()
    {
        producers.fetch_sub(1, std::memory_order_relaxed);
    }

    bool done()
    {
        return producers.load(std::memory_order_relaxed) == 0 && empty();
    }
};
