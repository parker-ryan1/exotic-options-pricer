#pragma once

#include <memory>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <string>
#include <chrono>
#include <vector>
#include "logger.hpp"

/**
 * @file memory_profiler.hpp
 * @brief Memory profiling and leak detection utilities
 * @author Quantitative Finance Team
 * @version 1.0
 * @date 2025
 * 
 * This memory profiler provides:
 * - Memory allocation tracking
 * - Leak detection
 * - Memory usage statistics
 * - Peak memory monitoring
 * - Thread-safe operation
 * - Integration with custom allocators
 */

namespace Utils {

/**
 * @brief Memory allocation information
 */
struct AllocationInfo {
    size_t size;                    ///< Size of allocation in bytes
    std::string file;               ///< Source file where allocation occurred
    int line;                       ///< Line number where allocation occurred
    std::string function;           ///< Function name where allocation occurred
    std::chrono::system_clock::time_point timestamp;  ///< When allocation occurred
    std::thread::id thread_id;      ///< Thread that made the allocation
    
    AllocationInfo() : size(0), line(0), timestamp(std::chrono::system_clock::now()) {}
    
    AllocationInfo(size_t s, const std::string& f, int l, const std::string& func)
        : size(s), file(f), line(l), function(func), 
          timestamp(std::chrono::system_clock::now()),
          thread_id(std::this_thread::get_id()) {}
};

/**
 * @brief Memory usage statistics
 */
struct MemoryStats {
    std::atomic<size_t> total_allocated{0};     ///< Total bytes allocated
    std::atomic<size_t> total_deallocated{0};   ///< Total bytes deallocated
    std::atomic<size_t> current_usage{0};       ///< Current memory usage
    std::atomic<size_t> peak_usage{0};          ///< Peak memory usage
    std::atomic<size_t> allocation_count{0};    ///< Number of allocations
    std::atomic<size_t> deallocation_count{0};  ///< Number of deallocations
    std::atomic<size_t> active_allocations{0};  ///< Current active allocations
    
    /**
     * @brief Get memory efficiency ratio
     * @return Ratio of deallocated to allocated memory
     */
    double efficiency_ratio() const {
        size_t allocated = total_allocated.load();
        return allocated > 0 ? static_cast<double>(total_deallocated.load()) / allocated : 0.0;
    }
    
    /**
     * @brief Get average allocation size
     * @return Average size per allocation
     */
    double average_allocation_size() const {
        size_t count = allocation_count.load();
        return count > 0 ? static_cast<double>(total_allocated.load()) / count : 0.0;
    }
};

/**
 * @brief Thread-safe memory profiler
 * 
 * This class tracks memory allocations and deallocations to detect leaks
 * and provide memory usage statistics. It's designed to be used in
 * production environments with minimal performance overhead.
 */
class MemoryProfiler {
private:
    static std::unique_ptr<MemoryProfiler> instance_;
    static std::mutex instance_mutex_;
    
    mutable std::mutex allocations_mutex_;
    std::unordered_map<void*, AllocationInfo> active_allocations_;
    MemoryStats stats_;
    Logger logger_;
    bool enabled_;
    bool track_call_stacks_;
    size_t max_tracked_allocations_;
    
    // Private constructor for singleton
    MemoryProfiler();

public:
    /**
     * @brief Get singleton instance
     * @return Reference to memory profiler instance
     */
    static MemoryProfiler& getInstance();
    
    /**
     * @brief Initialize memory profiler
     * @param enabled Enable memory tracking
     * @param track_call_stacks Enable call stack tracking (expensive)
     * @param max_tracked_allocations Maximum number of allocations to track
     */
    static void initialize(bool enabled = true, bool track_call_stacks = false, 
                          size_t max_tracked_allocations = 100000);
    
    /**
     * @brief Enable or disable memory profiling
     * @param enabled True to enable profiling
     */
    void setEnabled(bool enabled) { enabled_ = enabled; }
    
    /**
     * @brief Check if profiling is enabled
     * @return True if profiling is enabled
     */
    bool isEnabled() const { return enabled_; }
    
    /**
     * @brief Record memory allocation
     * @param ptr Pointer to allocated memory
     * @param size Size of allocation
     * @param file Source file name
     * @param line Line number
     * @param function Function name
     */
    void recordAllocation(void* ptr, size_t size, const char* file = nullptr, 
                         int line = 0, const char* function = nullptr);
    
    /**
     * @brief Record memory deallocation
     * @param ptr Pointer to deallocated memory
     * @return Size of deallocated memory (0 if not found)
     */
    size_t recordDeallocation(void* ptr);
    
    /**
     * @brief Get current memory statistics
     * @return Memory usage statistics
     */
    MemoryStats getStats() const { return stats_; }
    
    /**
     * @brief Get number of active allocations
     * @return Number of unfreed allocations
     */
    size_t getActiveAllocationCount() const;
    
    /**
     * @brief Check for memory leaks
     * @return Vector of leak information
     */
    std::vector<AllocationInfo> detectLeaks() const;
    
    /**
     * @brief Print memory usage report
     * @param log_level Log level to use for output
     */
    void printReport(LogLevel log_level = LogLevel::INFO) const;
    
    /**
     * @brief Print detailed leak report
     * @param max_leaks Maximum number of leaks to report
     */
    void printLeakReport(size_t max_leaks = 10) const;
    
    /**
     * @brief Reset all statistics
     */
    void reset();
    
    /**
     * @brief Get memory usage in human-readable format
     * @param bytes Number of bytes
     * @return Formatted string (e.g., "1.5 MB")
     */
    static std::string formatBytes(size_t bytes);
    
    /**
     * @brief Check if memory usage exceeds threshold
     * @param threshold_bytes Memory threshold in bytes
     * @return True if current usage exceeds threshold
     */
    bool exceedsThreshold(size_t threshold_bytes) const {
        return stats_.current_usage.load() > threshold_bytes;
    }
    
    /**
     * @brief Get memory usage trend
     * @param window_size Number of recent samples to consider
     * @return Positive for increasing, negative for decreasing, 0 for stable
     */
    double getUsageTrend(size_t window_size = 10) const;

private:
    void updatePeakUsage();
    std::string getCallStack() const;
};

/**
 * @brief RAII memory tracker for scoped allocations
 */
class ScopedMemoryTracker {
private:
    std::string scope_name_;
    size_t initial_usage_;
    std::chrono::high_resolution_clock::time_point start_time_;
    Logger logger_;

public:
    /**
     * @brief Start tracking memory for a scope
     * @param scope_name Name of the scope being tracked
     */
    explicit ScopedMemoryTracker(const std::string& scope_name);
    
    /**
     * @brief Destructor - reports memory usage for scope
     */
    ~ScopedMemoryTracker();
    
    /**
     * @brief Get current memory usage for this scope
     * @return Memory used since construction
     */
    size_t getCurrentUsage() const;
};

/**
 * @brief Custom allocator with memory tracking
 * 
 * This allocator can be used with STL containers to track their memory usage.
 */
template<typename T>
class TrackedAllocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    
    template<typename U>
    struct rebind {
        using other = TrackedAllocator<U>;
    };
    
    TrackedAllocator() = default;
    
    template<typename U>
    TrackedAllocator(const TrackedAllocator<U>&) {}
    
    pointer allocate(size_type n) {
        size_t size = n * sizeof(T);
        pointer ptr = static_cast<pointer>(std::malloc(size));
        
        if (!ptr) {
            throw std::bad_alloc();
        }
        
        MemoryProfiler::getInstance().recordAllocation(ptr, size, __FILE__, __LINE__, __FUNCTION__);
        return ptr;
    }
    
    void deallocate(pointer ptr, size_type) {
        MemoryProfiler::getInstance().recordDeallocation(ptr);
        std::free(ptr);
    }
    
    template<typename U, typename... Args>
    void construct(U* ptr, Args&&... args) {
        new(ptr) U(std::forward<Args>(args)...);
    }
    
    template<typename U>
    void destroy(U* ptr) {
        ptr->~U();
    }
    
    bool operator==(const TrackedAllocator&) const { return true; }
    bool operator!=(const TrackedAllocator&) const { return false; }
};

// Convenience type aliases
template<typename T>
using TrackedVector = std::vector<T, TrackedAllocator<T>>;

template<typename Key, typename Value>
using TrackedMap = std::unordered_map<Key, Value, std::hash<Key>, std::equal_to<Key>, 
                                     TrackedAllocator<std::pair<const Key, Value>>>;

} // namespace Utils

// Memory tracking macros
#ifdef ENABLE_MEMORY_PROFILING

#define TRACK_MEMORY_SCOPE(name) \
    Utils::ScopedMemoryTracker _memory_tracker(name)

#define TRACK_ALLOCATION(ptr, size) \
    Utils::MemoryProfiler::getInstance().recordAllocation(ptr, size, __FILE__, __LINE__, __FUNCTION__)

#define TRACK_DEALLOCATION(ptr) \
    Utils::MemoryProfiler::getInstance().recordDeallocation(ptr)

// Override global new/delete operators for automatic tracking
void* operator new(size_t size);
void* operator new[](size_t size);
void operator delete(void* ptr) noexcept;
void operator delete[](void* ptr) noexcept;
void operator delete(void* ptr, size_t size) noexcept;
void operator delete[](void* ptr, size_t size) noexcept;

#else

#define TRACK_MEMORY_SCOPE(name)
#define TRACK_ALLOCATION(ptr, size)
#define TRACK_DEALLOCATION(ptr)

#endif // ENABLE_MEMORY_PROFILING