#include "config.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <thread>

namespace Config {

ConfigManager::ConfigManager() : logger_("ConfigManager") {
    setDefaults();
}

ConfigManager& ConfigManager::getInstance() {
    static ConfigManager instance;
    return instance;
}

bool ConfigManager::initialize(const std::string& config_file_path) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    logger_.info("Initializing configuration system with file: {}", config_file_path);
    
    config_file_path_ = config_file_path;
    
    setDefaults();
    
    if (!config_file_path.empty()) {
        if (!loadFromFile(config_file_path)) {
            logger_.warning("Failed to load configuration file: {}, using defaults", config_file_path);
        }
    }
    
    loadEnvironmentOverrides();
    
    if (!validateConfiguration()) {
        logger_.error("Configuration validation failed");
        return false;
    }
    
    logger_.info("Configuration system initialized successfully");
    printConfiguration(Utils::LogLevel::DEBUG);
    
    return true;
}

void ConfigManager::setDefaults() {
    // Monte Carlo settings
    config_map_["monte_carlo.simulations"] = ConfigValue(100000);
    config_map_["monte_carlo.steps"] = ConfigValue(252);
    config_map_["monte_carlo.use_antithetic"] = ConfigValue(true);
    config_map_["monte_carlo.random_seed"] = ConfigValue(42);
    config_map_["monte_carlo.variance_reduction"] = ConfigValue(true);
    
    // Numerical settings
    config_map_["numerical.tolerance"] = ConfigValue(1e-6);
    config_map_["numerical.max_iterations"] = ConfigValue(1000);
    config_map_["numerical.use_high_precision"] = ConfigValue(false);
    
    // Logging settings
    config_map_["logging.level"] = ConfigValue("INFO");
    config_map_["logging.file"] = ConfigValue("exotic_options.log");
    config_map_["logging.console"] = ConfigValue(true);
    config_map_["logging.file_output"] = ConfigValue(true);
    config_map_["logging.max_files"] = ConfigValue(5);
    config_map_["logging.max_file_size_mb"] = ConfigValue(10);
    
    // Performance settings
    config_map_["performance.enable_logging"] = ConfigValue(true);
    config_map_["performance.enable_profiling"] = ConfigValue(false);
    config_map_["performance.profile_memory"] = ConfigValue(false);
    config_map_["performance.benchmark_mode"] = ConfigValue(false);
    
    // Threading settings
    config_map_["threading.enable_safety"] = ConfigValue(true);
    config_map_["threading.max_threads"] = ConfigValue(static_cast<int>(std::thread::hardware_concurrency()));
    config_map_["threading.enable_parallel_mc"] = ConfigValue(true);
    config_map_["threading.thread_pool_size"] = ConfigValue(4);
    
    // Memory management
    config_map_["memory.enable_profiling"] = ConfigValue(false);
    config_map_["memory.max_usage_mb"] = ConfigValue(1024);
    config_map_["memory.enable_leak_detection"] = ConfigValue(false);
    config_map_["memory.gc_frequency"] = ConfigValue(1000);
    
    // Exotic options specific settings
    config_map_["exotic.barrier_monitoring_frequency"] = ConfigValue(1.0);
    config_map_["exotic.use_control_variates"] = ConfigValue(false);
    config_map_["exotic.asian_averaging_points"] = ConfigValue(252);
    config_map_["exotic.digital_smoothing"] = ConfigValue(0.01);
    config_map_["exotic.lookback_monitoring_frequency"] = ConfigValue(1.0);
    config_map_["exotic.rainbow_correlation_threshold"] = ConfigValue(0.99);
    
    // Risk management
    config_map_["risk.var_confidence_95"] = ConfigValue(0.95);
    config_map_["risk.var_confidence_99"] = ConfigValue(0.99);
    config_map_["risk.enable_stress_testing"] = ConfigValue(true);
    config_map_["risk.max_position_size"] = ConfigValue(1000000.0);
    
    // Market data
    config_map_["market.default_risk_free_rate"] = ConfigValue(0.05);
    config_map_["market.default_dividend_yield"] = ConfigValue(0.0);
    config_map_["market.default_volatility"] = ConfigValue(0.2);
    config_map_["market.correlation_decay_factor"] = ConfigValue(0.95);
    
    // Validation settings
    config_map_["validation.enable_parameter_checks"] = ConfigValue(true);
    config_map_["validation.warn_extreme_values"] = ConfigValue(true);
    config_map_["validation.max_volatility"] = ConfigValue(5.0);
    config_map_["validation.max_time_to_expiry"] = ConfigValue(30.0);
    config_map_["validation.min_barrier_distance"] = ConfigValue(0.01);
}

bool ConfigManager::loadFromFile(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        logger_.error("Cannot open configuration file: {}", file_path);
        return false;
    }
    
    std::ostringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();
    
    if (content.empty()) {
        logger_.warning("Configuration file is empty: {}", file_path);
        return false;
    }
    
    return parseJsonContent(content);
}

bool ConfigManager::parseJsonContent(const std::string& json_content) {
    std::istringstream stream(json_content);
    std::string line;
    
    while (std::getline(stream, line)) {
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        
        if (line.empty() || line[0] == '#' || line.substr(0, 2) == "//") {
            continue;
        }
        
        if (line == "{" || line == "}" || line == "[" || line == "]") {
            continue;
        }
        
        if (!line.empty() && line.back() == ',') {
            line.pop_back();
        }
        
        size_t colon_pos = line.find(':');
        if (colon_pos == std::string::npos) {
            continue;
        }
        
        std::string key = line.substr(0, colon_pos);
        std::string value = line.substr(colon_pos + 1);
        
        if (key.front() == '"' && key.back() == '"') {
            key = key.substr(1, key.length() - 2);
        }
        if (value.front() == '"' && value.back() == '"') {
            value = value.substr(1, value.length() - 2);
        }
        
        if (value == "true" || value == "false") {
            config_map_[key] = ConfigValue(value == "true");
        } else if (value.find('.') != std::string::npos) {
            try {
                config_map_[key] = ConfigValue(std::stod(value));
            } catch (...) {
                config_map_[key] = ConfigValue(value);
            }
        } else {
            try {
                config_map_[key] = ConfigValue(std::stoi(value));
            } catch (...) {
                config_map_[key] = ConfigValue(value);
            }
        }
    }
    
    logger_.info("Loaded {} configuration values from JSON", config_map_.size());
    return true;
}

void ConfigManager::loadEnvironmentOverrides() {
    const char* env_vars[] = {
        "EXOTIC_OPTIONS_MONTE_CARLO_SIMULATIONS",
        "EXOTIC_OPTIONS_MONTE_CARLO_STEPS",
        "EXOTIC_OPTIONS_LOGGING_LEVEL",
        "EXOTIC_OPTIONS_LOGGING_FILE",
        "EXOTIC_OPTIONS_THREADING_MAX_THREADS",
        "EXOTIC_OPTIONS_MEMORY_MAX_USAGE_MB",
        nullptr
    };
    
    const char* config_keys[] = {
        "monte_carlo.simulations",
        "monte_carlo.steps",
        "logging.level",
        "logging.file",
        "threading.max_threads",
        "memory.max_usage_mb",
        nullptr
    };
    
    for (int i = 0; env_vars[i] != nullptr; ++i) {
        const char* env_value = std::getenv(env_vars[i]);
        if (env_value != nullptr) {
            logger_.info("Environment override: {} = {}", config_keys[i], env_value);
            
            std::string value_str(env_value);
            if (value_str == "true" || value_str == "false") {
                config_map_[config_keys[i]] = ConfigValue(value_str == "true");
            } else if (value_str.find('.') != std::string::npos) {
                try {
                    config_map_[config_keys[i]] = ConfigValue(std::stod(value_str));
                } catch (...) {
                    config_map_[config_keys[i]] = ConfigValue(value_str);
                }
            } else {
                try {
                    config_map_[config_keys[i]] = ConfigValue(std::stoi(value_str));
                } catch (...) {
                    config_map_[config_keys[i]] = ConfigValue(value_str);
                }
            }
        }
    }
}

bool ConfigManager::validateConfiguration() {
    bool is_valid = true;
    
    if (getInt("monte_carlo.simulations") <= 0) {
        logger_.error("Invalid monte_carlo.simulations: must be positive");
        is_valid = false;
    }
    
    if (getInt("monte_carlo.steps") <= 0) {
        logger_.error("Invalid monte_carlo.steps: must be positive");
        is_valid = false;
    }
    
    if (getDouble("numerical.tolerance") <= 0.0) {
        logger_.error("Invalid numerical.tolerance: must be positive");
        is_valid = false;
    }
    
    int max_threads = getInt("threading.max_threads");
    if (max_threads <= 0 || max_threads > 1000) {
        logger_.error("Invalid threading.max_threads: must be between 1 and 1000");
        is_valid = false;
    }
    
    if (getInt("memory.max_usage_mb") <= 0) {
        logger_.error("Invalid memory.max_usage_mb: must be positive");
        is_valid = false;
    }
    
    std::string log_level = getString("logging.level");
    if (log_level != "DEBUG" && log_level != "INFO" && log_level != "WARNING" && 
        log_level != "ERROR" && log_level != "CRITICAL") {
        logger_.error("Invalid logging.level: must be DEBUG, INFO, WARNING, ERROR, or CRITICAL");
        is_valid = false;
    }
    
    return is_valid;
}

std::string ConfigManager::getString(const std::string& key, const std::string& default_value) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = config_map_.find(key);
    if (it != config_map_.end()) {
        return it->second.getString();
    }
    
    return default_value;
}

int ConfigManager::getInt(const std::string& key, int default_value) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = config_map_.find(key);
    if (it != config_map_.end()) {
        try {
            return static_cast<int>(it->second);
        } catch (...) {
            logger_.warning("Failed to convert config value '{}' to int, using default", key);
        }
    }
    
    return default_value;
}

double ConfigManager::getDouble(const std::string& key, double default_value) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = config_map_.find(key);
    if (it != config_map_.end()) {
        try {
            return static_cast<double>(it->second);
        } catch (...) {
            logger_.warning("Failed to convert config value '{}' to double, using default", key);
        }
    }
    
    return default_value;
}

bool ConfigManager::getBool(const std::string& key, bool default_value) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = config_map_.find(key);
    if (it != config_map_.end()) {
        try {
            return static_cast<bool>(it->second);
        } catch (...) {
            logger_.warning("Failed to convert config value '{}' to bool, using default", key);
        }
    }
    
    return default_value;
}

void ConfigManager::set(const std::string& key, const ConfigValue& value) {
    std::lock_guard<std::mutex> lock(mutex_);
    config_map_[key] = value;
    logger_.debug("Configuration updated: {} = {}", key, value.getString());
}

bool ConfigManager::hasKey(const std::string& key) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return config_map_.find(key) != config_map_.end();
}

bool ConfigManager::saveToFile(const std::string& file_path) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::string output_path = file_path.empty() ? config_file_path_ : file_path;
    
    std::ofstream file(output_path);
    if (!file.is_open()) {
        logger_.error("Cannot open configuration file for writing: {}", output_path);
        return false;
    }
    
    file << "{\n";
    
    bool first = true;
    for (const auto& pair : config_map_) {
        if (!first) {
            file << ",\n";
        }
        
        file << "  \"" << pair.first << "\": ";
        
        if (pair.second.getType() == ValueType::STRING) {
            file << "\"" << pair.second.getString() << "\"";
        } else {
            file << pair.second.getString();
        }
        
        first = false;
    }
    
    file << "\n}\n";
    
    logger_.info("Configuration saved to: {}", output_path);
    return true;
}

bool ConfigManager::reload() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    logger_.info("Reloading configuration from: {}", config_file_path_);
    
    setDefaults();
    
    if (!config_file_path_.empty()) {
        loadFromFile(config_file_path_);
    }
    
    loadEnvironmentOverrides();
    
    return validateConfiguration();
}

std::vector<std::string> ConfigManager::getAllKeys() const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::vector<std::string> keys;
    keys.reserve(config_map_.size());
    
    for (const auto& pair : config_map_) {
        keys.push_back(pair.first);
    }
    
    return keys;
}

void ConfigManager::printConfiguration(Utils::LogLevel log_level) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!Utils::Logger::is_enabled(log_level)) {
        return;
    }
    
    logger_.info("Current Configuration ({} values):", config_map_.size());
    
    for (const auto& pair : config_map_) {
        if (log_level == Utils::LogLevel::DEBUG) {
            logger_.debug("  {} = {}", pair.first, pair.second.getString());
        } else if (log_level == Utils::LogLevel::INFO) {
            logger_.info("  {} = {}", pair.first, pair.second.getString());
        }
    }
}

} // namespace Config