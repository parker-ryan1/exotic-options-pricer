#include "test_framework.hpp"
#include <algorithm>
#include <iomanip>

namespace Testing {

std::string to_string(TestStatus status) {
    switch (status) {
        case TestStatus::PASSED:  return "PASSED";
        case TestStatus::FAILED:  return "FAILED";
        case TestStatus::SKIPPED: return "SKIPPED";
        case TestStatus::ERROR:   return "ERROR";
        default:                  return "UNKNOWN";
    }
}

TestResult TestCase::run() {
    TestResult result(name_);
    
    if (!enabled_) {
        result.status = TestStatus::SKIPPED;
        return result;
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        test_function_();
        result.status = TestStatus::PASSED;
    } catch (const AssertionFailure& e) {
        result.status = TestStatus::FAILED;
        result.error_message = e.what();
    } catch (const std::exception& e) {
        result.status = TestStatus::ERROR;
        result.error_message = std::string("Unexpected exception: ") + e.what();
    } catch (...) {
        result.status = TestStatus::ERROR;
        result.error_message = "Unknown exception thrown";
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    result.execution_time_ms = duration.count() / 1000.0;
    
    return result;
}

void TestSuite::addTest(const std::string& name, TestFunction func, bool enabled) {
    auto test_case = std::make_unique<TestCase>(name, func, enabled);
    test_cases_.push_back(std::move(test_case));
}

void TestSuite::addTest(const std::string& name, TestFunction func, 
                       const std::vector<std::string>& tags, bool enabled) {
    auto test_case = std::make_unique<TestCase>(name, func, enabled);
    for (const auto& tag : tags) {
        test_case->addTag(tag);
    }
    test_cases_.push_back(std::move(test_case));
}

TestSuiteStats TestSuite::runAll() {
    TestSuiteStats stats;
    
    logger_.info("Running test suite: {} ({} tests)", name_, test_cases_.size());
    
    for (const auto& test_case : test_cases_) {
        if (!shouldRunTest(*test_case)) {
            continue;
        }
        
        TestResult result = test_case->run();
        
        stats.total_tests++;
        stats.total_time_ms += result.execution_time_ms;
        stats.total_memory_bytes += result.memory_used_bytes;
        
        switch (result.status) {
            case TestStatus::PASSED:
                stats.passed_tests++;
                break;
            case TestStatus::FAILED:
                stats.failed_tests++;
                break;
            case TestStatus::SKIPPED:
                stats.skipped_tests++;
                break;
            case TestStatus::ERROR:
                stats.error_tests++;
                break;
        }
        
        if (verbose_output_ || result.status != TestStatus::PASSED) {
            printTestResult(result);
        }
    }
    
    printSummary(stats);
    return stats;
}

TestSuiteStats TestSuite::runWithTag(const std::string& tag) {
    TestSuiteStats stats;
    
    logger_.info("Running tests with tag '{}' in suite: {}", tag, name_);
    
    for (const auto& test_case : test_cases_) {
        if (!test_case->hasTag(tag) || !shouldRunTest(*test_case)) {
            continue;
        }
        
        TestResult result = test_case->run();
        
        stats.total_tests++;
        stats.total_time_ms += result.execution_time_ms;
        
        switch (result.status) {
            case TestStatus::PASSED:
                stats.passed_tests++;
                break;
            case TestStatus::FAILED:
                stats.failed_tests++;
                break;
            case TestStatus::SKIPPED:
                stats.skipped_tests++;
                break;
            case TestStatus::ERROR:
                stats.error_tests++;
                break;
        }
        
        if (verbose_output_ || result.status != TestStatus::PASSED) {
            printTestResult(result);
        }
    }
    
    printSummary(stats);
    return stats;
}

TestResult TestSuite::runTest(const std::string& test_name) {
    for (const auto& test_case : test_cases_) {
        if (test_case->getName() == test_name) {
            TestResult result = test_case->run();
            printTestResult(result);
            return result;
        }
    }
    
    TestResult result(test_name);
    result.status = TestStatus::ERROR;
    result.error_message = "Test not found: " + test_name;
    return result;
}

std::vector<std::string> TestSuite::getTestNames() const {
    std::vector<std::string> names;
    names.reserve(test_cases_.size());
    
    for (const auto& test_case : test_cases_) {
        names.push_back(test_case->getName());
    }
    
    return names;
}

bool TestSuite::shouldRunTest(const TestCase& test) const {
    if (!test.isEnabled()) {
        return false;
    }
    
    for (const auto& disabled_tag : disabled_tags_) {
        if (test.hasTag(disabled_tag)) {
            return false;
        }
    }
    
    if (!enabled_tags_.empty()) {
        bool has_enabled_tag = false;
        for (const auto& enabled_tag : enabled_tags_) {
            if (test.hasTag(enabled_tag)) {
                has_enabled_tag = true;
                break;
            }
        }
        if (!has_enabled_tag) {
            return false;
        }
    }
    
    return true;
}

void TestSuite::printTestResult(const TestResult& result) const {
    std::string status_color;
    switch (result.status) {
        case TestStatus::PASSED:
            status_color = "\033[32m";  // Green
            break;
        case TestStatus::FAILED:
            status_color = "\033[31m";  // Red
            break;
        case TestStatus::SKIPPED:
            status_color = "\033[33m";  // Yellow
            break;
        case TestStatus::ERROR:
            status_color = "\033[35m";  // Magenta
            break;
    }
    
    std::cout << status_color << "[" << std::setw(7) << to_string(result.status) << "]\033[0m "
              << std::setw(40) << std::left << result.test_name
              << " (" << std::fixed << std::setprecision(2) << result.execution_time_ms << "ms)";
    
    if (!result.error_message.empty()) {
        std::cout << "\n    " << result.error_message;
    }
    
    std::cout << std::endl;
}

void TestSuite::printSummary(const TestSuiteStats& stats) const {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Test Suite: " << name_ << " - Summary" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    std::cout << "Total Tests:   " << stats.total_tests << std::endl;
    std::cout << "Passed:        " << "\033[32m" << stats.passed_tests << "\033[0m" << std::endl;
    std::cout << "Failed:        " << "\033[31m" << stats.failed_tests << "\033[0m" << std::endl;
    std::cout << "Errors:        " << "\033[35m" << stats.error_tests << "\033[0m" << std::endl;
    std::cout << "Skipped:       " << "\033[33m" << stats.skipped_tests << "\033[0m" << std::endl;
    std::cout << "Success Rate:  " << std::fixed << std::setprecision(1) << stats.success_rate() << "%" << std::endl;
    std::cout << "Total Time:    " << std::fixed << std::setprecision(2) << stats.total_time_ms << "ms" << std::endl;
    
    if (stats.total_memory_bytes > 0) {
        std::cout << "Memory Used:   " << (stats.total_memory_bytes / 1024) << " KB" << std::endl;
    }
    
    std::cout << std::string(60, '=') << std::endl;
}

TestRegistry& TestRegistry::getInstance() {
    static TestRegistry instance;
    return instance;
}

void TestRegistry::registerSuite(std::unique_ptr<TestSuite> suite) {
    logger_.info("Registering test suite: {}", suite->getName());
    test_suites_.push_back(std::move(suite));
}

TestSuiteStats TestRegistry::runAllSuites() {
    TestSuiteStats combined_stats;
    
    logger_.info("Running all test suites ({} suites)", test_suites_.size());
    
    for (const auto& suite : test_suites_) {
        TestSuiteStats suite_stats = suite->runAll();
        
        combined_stats.total_tests += suite_stats.total_tests;
        combined_stats.passed_tests += suite_stats.passed_tests;
        combined_stats.failed_tests += suite_stats.failed_tests;
        combined_stats.skipped_tests += suite_stats.skipped_tests;
        combined_stats.error_tests += suite_stats.error_tests;
        combined_stats.total_time_ms += suite_stats.total_time_ms;
        combined_stats.total_memory_bytes += suite_stats.total_memory_bytes;
    }
    
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "OVERALL TEST RESULTS" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    
    std::cout << "Test Suites:   " << test_suites_.size() << std::endl;
    std::cout << "Total Tests:   " << combined_stats.total_tests << std::endl;
    std::cout << "Passed:        " << "\033[32m" << combined_stats.passed_tests << "\033[0m" << std::endl;
    std::cout << "Failed:        " << "\033[31m" << combined_stats.failed_tests << "\033[0m" << std::endl;
    std::cout << "Errors:        " << "\033[35m" << combined_stats.error_tests << "\033[0m" << std::endl;
    std::cout << "Skipped:       " << "\033[33m" << combined_stats.skipped_tests << "\033[0m" << std::endl;
    std::cout << "Success Rate:  " << std::fixed << std::setprecision(1) << combined_stats.success_rate() << "%" << std::endl;
    std::cout << "Total Time:    " << std::fixed << std::setprecision(2) << combined_stats.total_time_ms << "ms" << std::endl;
    
    if (combined_stats.total_memory_bytes > 0) {
        std::cout << "Memory Used:   " << (combined_stats.total_memory_bytes / 1024) << " KB" << std::endl;
    }
    
    std::cout << std::string(80, '=') << std::endl;
    
    return combined_stats;
}

TestSuiteStats TestRegistry::runSuite(const std::string& suite_name) {
    for (const auto& suite : test_suites_) {
        if (suite->getName() == suite_name) {
            return suite->runAll();
        }
    }
    
    logger_.error("Test suite not found: {}", suite_name);
    return TestSuiteStats{};
}

std::vector<std::string> TestRegistry::getSuiteNames() const {
    std::vector<std::string> names;
    names.reserve(test_suites_.size());
    
    for (const auto& suite : test_suites_) {
        names.push_back(suite->getName());
    }
    
    return names;
}

void TestRegistry::printDiscovery() const {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "TEST DISCOVERY" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int total_tests = 0;
    
    for (const auto& suite : test_suites_) {
        std::cout << "\nSuite: " << suite->getName() << " (" << suite->getTestCount() << " tests)" << std::endl;
        
        auto test_names = suite->getTestNames();
        for (const auto& test_name : test_names) {
            std::cout << "  - " << test_name << std::endl;
            total_tests++;
        }
    }
    
    std::cout << "\nTotal: " << test_suites_.size() << " suites, " << total_tests << " tests" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
}

} // namespace Testing