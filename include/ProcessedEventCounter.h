#ifndef REST_PROCESSEDEVENTCOUNTER_H
#define REST_PROCESSEDEVENTCOUNTER_H

#include <atomic>

// Accounting only: no event payload is published through this counter.
class ProcessedEventCounter {
   public:
    void Increment() { fPending.fetch_add(1, std::memory_order_relaxed); }
    int Take() { return fPending.exchange(0, std::memory_order_relaxed); }

   private:
    std::atomic<int> fPending{0};
};

#endif
