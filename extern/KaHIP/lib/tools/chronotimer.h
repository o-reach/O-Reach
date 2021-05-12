/******************************************************************************
 * chronotimer.h
 *****************************************************************************/

#ifndef CHRONOTIMER_H
#define CHRONOTIMER_H

#include <chrono>

class ChronoTimer {
        public:
                ChronoTimer() : start_time(std::chrono::steady_clock::now()) { }

                void restart() {
                    start_time = std::chrono::steady_clock::now();
                }

                template<typename time_unit = std::chrono::seconds>
                long long int elapsed_integral() {
                    static std::chrono::steady_clock::time_point stop_time;
                    stop_time = std::chrono::steady_clock::now();
                    return std::chrono::duration_cast<time_unit>(stop_time - start_time).count();
                }

                template<typename unit = std::ratio<1>>
                double elapsed() {
                    static std::chrono::steady_clock::time_point stop_time;
                    stop_time = std::chrono::steady_clock::now();
                    std::chrono::duration<double, unit> delta = stop_time - start_time;
                    return delta.count();
                }

        private:
                std::chrono::steady_clock::time_point start_time;
};

#endif /* CHRONOTIMER_H */
