#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <chrono>
#include <iomanip>
#include <iostream>

class Stopwatch {
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

public:
    Stopwatch() { reset(); }
    void reset() { start_time = std::chrono::high_resolution_clock::now(); }
    double count() const
    {
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - start_time;
        return duration.count();
    }

    friend std::ostream &operator<<(std::ostream &os, const Stopwatch &rhs)
    {
        return os << std::fixed << std::setprecision(2) << rhs.count();
    }
};

#endif
