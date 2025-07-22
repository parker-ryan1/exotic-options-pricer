#include "logger.hpp"
#include <iostream>
#include <filesystem>
#include <algorithm>

namespace Utils {

// Static member definitions
std::mutex Logger::global_mutex_;
std::ofstream Logger::log_file_;
LogLevel Logger::min_level_ = LogLevel::INFO;
bool Logger::console_output_ = true;
bool Logger::file_output_ = true;
std::string Logger::log_filename_ = "exotic_options.log";
size_t Logger::max_file_size_ = 10 * 1024 * 1024;
size_t Logger::current_file_size_ = 0;
int Logger::max_log_files_ = 5;

Utils::Logger g_logger("Global");

std::string to_string(LogLevel level) {
    switch (level) {
        case LogLevel::DEBUG:    return "DEBUG";
        case LogLevel::INFO:     return "INFO";
        case LogLevel::WARNING:  return "WARNING";
        case LogLevel::ERROR:    return "ERROR";
        case LogLevel::CRITICAL: return "CRITICAL";
        default:                 return "UNKNOWN";
    }
}

Logger::Logger(const std::string& component_name) 
    : component_name_(component_name) {
    static std::once_flag init_flag;
    std::call_once(init_flag, []() {
        configure();
    });
}

Logger::~Logger() {
    flush();
}

void Logger::configure(
    LogLevel min_level,
    bool console_output,
    bool file_output,
    const std::string& log_filename,
    size_t max_file_size,
    int max_log_files) {
    
    std::lock_guard<std::mutex> lock(global_mutex_);
    
    min_level_ = min_level;
    console_output_ = console_output;
    file_output_ = file_output;
    log_filename_ = log_filename;
    max_file_size_ = max_file_size;
    max_log_files_ = max_log_files;
    
    if (log_file_.is_open()) {
        log_file_.close();
    }
    
    if (file_output_) {
        log_file_.open(log_filename_, std::ios::app);
        if (log_file_.is_open()) {
            log_file_.seekp(0, std::ios::end);
            current_file_size_ = log_file_.tellp();
            
            std::string config_msg = "Logger configured - Level: " + to_string(min_level_) +
                                   ", Console: " + (console_output_ ? "ON" : "OFF") +
                                   ", File: " + (file_output_ ? "ON" : "OFF");
            
            log_file_ << get_timestamp() << " [INFO] [Logger] " << config_msg << std::endl;
            current_file_size_ += config_msg.length() + 50;
        }
    }
}

std::string Logger::get_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    
    return oss.str();
}

std::string Logger::get_thread_id() {
    std::ostringstream oss;
    oss << std::this_thread::get_id();
    return oss.str();
}

void Logger::rotate_log_files() {
    if (!file_output_ || !log_file_.is_open()) {
        return;
    }
    
    log_file_.close();
    
    for (int i = max_log_files_ - 1; i > 0; --i) {
        std::string old_name = log_filename_ + "." + std::to_string(i);
        std::string new_name = log_filename_ + "." + std::to_string(i + 1);
        
        if (std::filesystem::exists(old_name)) {
            if (i == max_log_files_ - 1) {
                std::filesystem::remove(new_name);
            }
            std::filesystem::rename(old_name, new_name);
        }
    }
    
    if (std::filesystem::exists(log_filename_)) {
        std::filesystem::rename(log_filename_, log_filename_ + ".1");
    }
    
    log_file_.open(log_filename_, std::ios::out);
    current_file_size_ = 0;
    
    if (log_file_.is_open()) {
        std::string rotation_msg = "Log file rotated";
        log_file_ << get_timestamp() << " [INFO] [Logger] " << rotation_msg << std::endl;
        current_file_size_ += rotation_msg.length() + 50;
    }
}

void Logger::write_log(LogLevel level, const std::string& message) {
    std::lock_guard<std::mutex> lock(global_mutex_);
    
    std::string timestamp = get_timestamp();
    std::string thread_id = get_thread_id();
    std::string level_str = to_string(level);
    
    std::ostringstream log_line;
    log_line << timestamp << " [" << level_str << "] [" << component_name_ 
             << "] [T:" << thread_id << "] " << message;
    
    std::string full_message = log_line.str();
    
    if (console_output_) {
        if (level >= LogLevel::ERROR) {
            std::cerr << full_message << std::endl;
        } else {
            std::cout << full_message << std::endl;
        }
    }
    
    if (file_output_ && log_file_.is_open()) {
        log_file_ << full_message << std::endl;
        current_file_size_ += full_message.length() + 1;
        
        if (current_file_size_ > max_file_size_) {
            rotate_log_files();
        }
    }
}

void Logger::format_string(std::ostringstream& oss, const std::string& format) {
    oss << format;
}

void Logger::flush() {
    std::lock_guard<std::mutex> lock(global_mutex_);
    
    if (console_output_) {
        std::cout.flush();
        std::cerr.flush();
    }
    
    if (file_output_ && log_file_.is_open()) {
        log_file_.flush();
    }
}

PerformanceTimer::PerformanceTimer(Logger& logger, const std::string& operation_name, LogLevel log_level)
    : logger_(logger), operation_name_(operation_name), log_level_(log_level) {
    start_time_ = std::chrono::high_resolution_clock::now();
    
    if (logger_.is_enabled(log_level_)) {
        logger_.debug("Starting operation: {}", operation_name_);
    }
}

PerformanceTimer::~PerformanceTimer() {
    double elapsed = elapsed_ms();
    
    if (logger_.is_enabled(log_level_)) {
        if (log_level_ == LogLevel::DEBUG) {
            logger_.debug("Operation '{}' completed in {:.3f}ms", operation_name_, elapsed);
        } else if (log_level_ == LogLevel::INFO) {
            logger_.info("Operation '{}' completed in {:.3f}ms", operation_name_, elapsed);
        } else if (log_level_ == LogLevel::WARNING) {
            logger_.warning("Operation '{}' completed in {:.3f}ms", operation_name_, elapsed);
        } else if (log_level_ == LogLevel::ERROR) {
            logger_.error("Operation '{}' completed in {:.3f}ms", operation_name_, elapsed);
        } else if (log_level_ == LogLevel::CRITICAL) {
            logger_.critical("Operation '{}' completed in {:.3f}ms", operation_name_, elapsed);
        }
    }
}

double PerformanceTimer::elapsed_ms() const {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time_);
    return duration.count() / 1000.0;
}

} // namespace Utils