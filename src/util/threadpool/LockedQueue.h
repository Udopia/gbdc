#pragma once
#include <queue>
#include <mutex>
#include <list>
#include <atomic>
#include <condition_variable>

template <typename T>
class LockedQueue
{
private:
    std::queue<T> q;
    std::mutex m;

public:
    explicit LockedQueue() : q(), m() {}

    void push(T &&t)
    {
        std::unique_lock<std::mutex> l(m);
        q.push(t);
    }

    void push(const T &t)
    {
        std::unique_lock<std::mutex> l(m);
        q.push(t);
    }

    bool try_pop(T &out)
    {
        std::unique_lock<std::mutex> l(m);
        if (empty())
            return false;
        out = q.front();
        q.pop();
        return true;
    }

    T pop()
    {
        std::unique_lock<std::mutex> l(m);
        T t = q.front();
        q.pop();
        return t;
    }

    bool empty()
    {
        return q.empty();
    }
};
