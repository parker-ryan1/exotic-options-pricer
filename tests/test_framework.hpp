#pragma once

#include <string>
#include <vector>
#include <functional>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include "../src/utils/logger.hpp"

/**
 * @file test_framework.hpp
 * @brief Testing framework for exotic options pricing system
 * @author Quantitative Finance Team
 * @version 1.0
 * @date 2025
 */

namespace Testing {

enum class TestStatus {
    PASSED,
    FAILED,
    SKIPPED,
    ERROR
};

std::string to_string(TestStatus status);

struct TestResult {
    std::string test_name;
    TestStatus status;
    std::string error_message;
    double execution_time_ms;
    size_t memory_used_bytes;
    
    TestResult(const std::string& name) 
        : test_name(name), status(TestStatus::PASSED), execution_time_ms(0.0), memory_used_bytes(0) {}
};

struct TestSuiteStats {
    int total_tests = 0;
    int passed_tests = 0;
    int failed_tests = 0;
    int skipped_tests = 0;
    int error_tests = 0;
    double total_time_ms = 0.0;
    size_t total_memory_bytes = 0;
    
    double success_rate() const {
        return total_tests > 0 ? (static_cast<double>(passed_tests) / total_tests) * 100.0 : 0.0;
    }
};

class AssertionFailure : public std::exception {
private:
    std::string message_;
    
public:
    explicit AssertionFailure(const std::string& message) : message_(message) {}
    const char* what() const noexcept override { return message_.c_str(); }
};

class TestFixture {
public:
    virtual ~TestFixture() = default;
    virtual void setUp() {}
    virtual void tearDown() {}
    
protected:
    Utils::Logger logger_{"TestFixture"};
};

using TestFunction = std::function<void()>;

template<typename FixtureType>
using TestMethodFunction = std::function<void(FixtureType&)>;

class TestCase {
private:
    std::string name_;
    TestFunction test_function_;
    bool enabled_;
    std::vector<std::string> tags_;
    
public:
    TestCase(const std::string& name, TestFunction func, bool enabled = true)
        : name_(name), test_function_(func), enabled_(enabled) {}
    
    const std::string& getName() const { return name_; }
    bool isEnabled() const { return enabled_; }
    void setEnabled(bool enabled) { enabled_ = enabled; }
    
    void addTag(const std::string& tag) { tags_.push_back(tag); }
    const std::vector<std::string>& getTags() const { return tags_; }
    bool hasTag(const std::string& tag) const {
        return std::find(tags_.begin(), tags_.end(), tag) != tags_.end();
    }
    
    TestResult run();
};

class TestSuite {
private:
    std::string name_;
    std::vector<std::unique_ptr<TestCase>> test_cases_;
    Utils::Logger logger_;
    bool verbose_output_;
    std::vector<std::string> enabled_tags_;
    std::vector<std::string> disabled_tags_;
    
public:
    explicit TestSuite(const std::string& name) 
        : name_(name), logger_("TestSuite::" + name), verbose_output_(false) {}
    
    void addTest(const std::string& name, TestFunction func, bool enabled = true);
    void addTest(const std::string& name, TestFunction func, 
                const std::vector<std::string>& tags, bool enabled = true);
    
    template<typename FixtureType>
    void addTestMethod(const std::string& name, TestMethodFunction<FixtureType> method, bool enabled = true) {
        auto wrapper = [method]() {
            FixtureType fixture;
            fixture.setUp();
            try {
                method(fixture);
                fixture.tearDown();
            } catch (...) {
                fixture.tearDown();
                throw;
            }
        };
        addTest(name, wrapper, enabled);
    }
    
    TestSuiteStats runAll();
    TestSuiteStats runWithTag(const std::string& tag);
    TestResult runTest(const std::string& test_name);
    
    void setVerbose(bool verbose) { verbose_output_ = verbose; }
    void enableTags(const std::vector<std::string>& tags) { enabled_tags_ = tags; }
    void disableTags(const std::vector<std::string>& tags) { disabled_tags_ = tags; }
    
    size_t getTestCount() const { return test_cases_.size(); }
    std::vector<std::string> getTestNames() const;
    const std::string& getName() const { return name_; }

private:
    bool shouldRunTest(const TestCase& test) const;
    void printTestResult(const TestResult& result) const;
    void printSummary(const TestSuiteStats& stats) const;
};

class TestRegistry {
private:
    std::vector<std::unique_ptr<TestSuite>> test_suites_;
    Utils::Logger logger_{"TestRegistry"};
    
    TestRegistry() = default;
    
public:
    static TestRegistry& getInstance();
    
    void registerSuite(std::unique_ptr<TestSuite> suite);
    TestSuiteStats runAllSuites();
    TestSuiteStats runSuite(const std::string& suite_name);
    std::vector<std::string> getSuiteNames() const;
    void printDiscovery() const;
};

class Benchmark {
private:
    std::string name_;
    std::chrono::high_resolution_clock::time_point start_time_;
    Utils::Logger logger_{"Benchmark"};
    
public:
    explicit Benchmark(const std::string& name) : name_(name) {
        start_time_ = std::chrono::high_resolution_clock::now();
        logger_.info("Starting benchmark: {}", name_);
    }
    
    ~Benchmark() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time_);
        double ms = duration.count() / 1000.0;
        logger_.info("Benchmark '{}' completed in {:.3f}ms", name_, ms);
    }
    
    double elapsed_ms() const {
        auto current_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time_);
        return duration.count() / 1000.0;
    }
};

// Assertion macros
#define ASSERT_TRUE(condition) \
    do { \
        if (!(condition)) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #condition " at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_FALSE(condition) \
    do { \
        if (condition) { \
            std::ostringstream oss; \
            oss << "Assertion failed: !(" #condition ") at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_EQ(expected, actual) \
    do { \
        if ((expected) != (actual)) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #expected " == " #actual \
                << " (expected: " << (expected) << ", actual: " << (actual) \
                << ") at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_NE(expected, actual) \
    do { \
        if ((expected) == (actual)) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #expected " != " #actual \
                << " (both values: " << (expected) \
                << ") at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_NEAR(expected, actual, tolerance) \
    do { \
        double diff = std::abs(static_cast<double>(expected) - static_cast<double>(actual)); \
        if (diff > (tolerance)) { \
            std::ostringstream oss; \
            oss << "Assertion failed: |" #expected " - " #actual "| <= " #tolerance \
                << " (expected: " << (expected) << ", actual: " << (actual) \
                << ", diff: " << diff << ", tolerance: " << (tolerance) \
                << ") at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_GT(left, right) \
    do { \
        if (!((left) > (right))) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #left " > " #right \
                << " (" << (left) << " > " << (right) \
                << ") at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_GE(left, right) \
    do { \
        if (!((left) >= (right))) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #left " >= " #right \
                << " (" << (left) << " >= " << (right) \
                << ") at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_LT(left, right) \
    do { \
        if (!((left) < (right))) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #left " < " #right \
                << " (" << (left) << " < " << (right) \
                << ") at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_LE(left, right) \
    do { \
        if (!((left) <= (right))) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #left " <= " #right \
                << " (" << (left) << " <= " << (right) \
                << ") at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_THROWS(statement, exception_type) \
    do { \
        bool caught = false; \
        try { \
            statement; \
        } catch (const exception_type&) { \
            caught = true; \
        } catch (...) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #statement " should throw " #exception_type \
                << " but threw different exception at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
        if (!caught) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #statement " should throw " #exception_type \
                << " but didn't throw at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define ASSERT_NO_THROW(statement) \
    do { \
        try { \
            statement; \
        } catch (...) { \
            std::ostringstream oss; \
            oss << "Assertion failed: " #statement " should not throw" \
                << " at " << __FILE__ << ":" << __LINE__; \
            throw Testing::AssertionFailure(oss.str()); \
        } \
    } while(0)

#define BENCHMARK(name) \
    Testing::Benchmark _benchmark(name)

} // namespace Testing