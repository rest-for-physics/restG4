#include "ProcessedEventCounter.h"

#include <atomic>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {
void CheckDrainBoundaries() {
    ProcessedEventCounter counter;
    counter.Increment();
    const auto before = counter.Take();
    counter.Increment();
    if (before + counter.Take() != 2 || counter.Take() != 0)
        throw std::runtime_error("An increment around a drain was lost or duplicated");
    counter.Increment();
    counter.Increment();
    if (counter.Take() != 2 || counter.Take() != 0)
        throw std::runtime_error("Drain did not consume exactly the pending events");
}

void CheckConcurrentAccounting() {
    // Exercise independent event producers and competing periodic/final drains.
    // Assert conservation of known work, without predicting a thread schedule.
    for (const int producers : {1, 4, 8}) {
        for (const int drains : {1, 2}) {
            ProcessedEventCounter counter;
            std::atomic<int> ready{0}, done{0};
            std::atomic<bool> start{false};
            std::atomic<long long> observed{0};
            constexpr int events = 250000;
            std::vector<std::thread> threads;
            for (int producer = 0; producer < producers; ++producer) {
                threads.emplace_back([&] {
                    ready.fetch_add(1);
                    while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                    for (int event = 0; event < events; ++event) {
                        counter.Increment();
                        if ((event & 1023) == 0) std::this_thread::yield();
                    }
                    done.fetch_add(1, std::memory_order_release);
                });
            }
            for (int drain = 0; drain < drains; ++drain) {
                threads.emplace_back([&] {
                    ready.fetch_add(1);
                    while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                    while (done.load(std::memory_order_acquire) != producers) {
                        observed.fetch_add(counter.Take(), std::memory_order_relaxed);
                        std::this_thread::yield();
                    }
                    observed.fetch_add(counter.Take(), std::memory_order_relaxed);
                });
            }
            while (ready.load() != producers + drains) std::this_thread::yield();
            start.store(true, std::memory_order_release);
            for (auto& thread : threads) thread.join();
            const auto actual = observed.load() + counter.Take();
            if (actual != static_cast<long long>(producers) * events || counter.Take() != 0)
                throw std::runtime_error("Concurrent event accounting lost or duplicated parents");
        }
    }
}
}  // namespace

#ifdef REST_EVENT_COUNTER_STANDALONE
int main() {
    CheckDrainBoundaries();
    CheckConcurrentAccounting();
    std::cout << "PASS drain boundaries and six concurrent conservation scenarios (6500000 increments)\n";
}
#else
#include <gtest/gtest.h>
TEST(ProcessedEventCounter, PreservesDrainBoundaries) { EXPECT_NO_THROW(CheckDrainBoundaries()); }
TEST(ProcessedEventCounter, ConservesConcurrentEvents) { EXPECT_NO_THROW(CheckConcurrentAccounting()); }
#endif
