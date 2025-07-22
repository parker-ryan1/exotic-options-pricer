#pragma once

#include <string>
#include <map>
#include <memory>
#include <mutex>
#include <fstream>
#include "../utils/logger.hpp"

/**
 * @file config.hpp
 * @brief Configuration management for exotic options pricing system
 * @author Quantitative Finance Team
 * @version 1.0
 * @date 2025
 */

namespace Config {

enum class ValueType {
    STRING,
    INTEGER,
    DOUBLE,
    BOOLEAN
};

class ConfigValue {
private:
    std::string string_value_;
    ValueType type_;
    
public:
    ConfigValue() : type_(ValueType::STRING) {}
    explicit ConfigValue(const std::string& value) : string_value_(value), type_(ValueType::STRING) {}
    explicit ConfigValue(int value) : string_value_(std::to_string(value)), type_(ValueType::INTEGER) {}
    explicit ConfigValue(double value) : string_value_(std::to_string(value)), type_(ValueType::DOUBLE) {}
    explicit ConfigValue(bool value) : string_value_(value ? "true" : "false"), type_(ValueType::BOOLEAN) {}
    
    operator std::string() const { return string_value_; }
    operator int() const { return std::stoi(string_value_); }
    operator double() const { return std::stod(string_value_); }
    operator bool() const { return string_value_ == "true" || string_value_ == "1"; }
    
    ValueType getType() const { return type_; }
    const std::string& getString() const { return string_value_; }
};

class ConfigManager {
private:
    mutable std::mutex mutex_;
    std::map<std::string, ConfigValue> config_map_;
    std::string config_file_path_;
    Utils::Logger logger_;
    
    ConfigManager();
    
    bool loadFromFile(const std::string& file_path);
    void loadEnvironmentOverrides();
    bool parseJsonContent(const std::string& json_content);
    void setDefaults();
    bool validateConfiguration();

public:
    static ConfigManager& getInstance();
    
    ConfigManager(const ConfigManager&) = delete;
    ConfigManager& operator=(const ConfigManager&) = delete;
    
    bool initialize(const std::string& config_file_path = "config.json");
    
    std::string getString(const std::string& key, const std::string& default_value = "") const;
    int getInt(const std::string& key, int default_value = 0) const;
    double getDouble(const std::string& key, double default_value = 0.0) const;
    bool getBool(const std::string& key, bool default_value = false) const;
    
    void set(const std::string& key, const ConfigValue& value);
    bool hasKey(const std::string& key) const;
    bool saveToFile(const std::string& file_path = "") const;
    bool reload();
    std::vector<std::string> getAllKeys() const;
    void printConfiguration(Utils::LogLevel log_level = Utils::LogLevel::INFO) const;
    
    // Convenience methods
    int getMonteCarloSimulations() const { return getInt("monte_carlo.simulations", 100000); }
    int getMonteCarloSteps() const { return getInt("monte_carlo.steps", 252); }
    bool getUseAntitheticVariates() const { return getBool("monte_carlo.use_antithetic", true); }
    int getRandomSeed() const { return getInt("monte_carlo.random_seed", 42); }
    
    double getConvergenceTolerance() const { return getDouble("numerical.tolerance", 1e-6); }
    int getMaxIterations() const { return getInt("numerical.max_iterations", 1000); }
    
    std::string getLogLevel() const { return getString("logging.level", "INFO"); }
    std::string getLogFile() const { return getString("logging.file", "exotic_options.log"); }
    bool getLogToConsole() const { return getBool("logging.console", true); }
    bool getLogToFile() const { return getBool("logging.file_output", true); }
    int getMaxLogFiles() const { return getInt("logging.max_files", 5); }
    size_t getMaxLogFileSize() const { return static_cast<size_t>(getInt("logging.max_file_size_mb", 10)) * 1024 * 1024; }
    
    bool getEnableThreadSafety() const { return getBool("threading.enable_safety", true); }
    int getMaxThreads() const { return getInt("threading.max_threads", std::thread::hardware_concurrency()); }
    bool getEnableParallelMC() const { return getBool("threading.enable_parallel_mc", true); }
    
    bool getEnableMemoryProfiling() const { return getBool("memory.enable_profiling", false); }
    size_t getMaxMemoryUsageMB() const { return static_cast<size_t>(getInt("memory.max_usage_mb", 1024)); }
    bool getEnableLeakDetection() const { return getBool("memory.enable_leak_detection", false); }
    
    bool getEnablePerformanceLogging() const { return getBool("performance.enable_logging", true); }
    bool getEnableProfiling() const { return getBool("performance.enable_profiling", false); }
    
    // Exotic options specific settings
    double getBarrierMonitoringFrequency() const { return getDouble("exotic.barrier_monitoring_frequency", 1.0); }
    bool getUseControlVariates() const { return getBool("exotic.use_control_variates", false); }
    int getAsianAveragingPoints() const { return getInt("exotic.asian_averaging_points", 252); }
    double getDigitalSmoothingParameter() const { return getDouble("exotic.digital_smoothing", 0.01); }
};

#define Config ConfigManager::getInstance()

} // namespace Config