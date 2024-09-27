#ifndef LIBDGP_TIMER_H
#define LIBDGP_TIMER_H
#include <chrono>

namespace dgp
{
    template <typename TimeT = std::chrono::milliseconds>
    class Timer {
    public:
        Timer() {
            start = std::chrono::system_clock::now();
        }

        size_t value() const {
            auto now = std::chrono::system_clock::now();
            auto duration = std::chrono::duration_cast<TimeT>(now - start);
            return (size_t) duration.count();
        }

        size_t reset() {
            auto now = std::chrono::system_clock::now();
            auto duration = std::chrono::duration_cast<TimeT>(now - start);
            start = now;
            return (size_t) duration.count();
        }

        void beginStage(const std::string &name) {
            reset();
            std::cout.flush();
        }

        void endStage(const std::string &str = "") {
            std::cout << "Done. took " << value() << "ms. " << std::endl;
        }
    private:
        std::chrono::system_clock::time_point start;
    };
}

#endif //LIBDGP_TIMER_H
