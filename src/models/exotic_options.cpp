#include "exotic_options.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <chrono>

namespace ExoticOptions {

// Thread-local RNG initialization
thread_local std::mt19937 ThreadSafeRNG::generator_;
thread_local std::normal_distribution<double> ThreadSafeRNG::normal_dist_(0.0, 1.0);
thread_local std::uniform_real_distribution<double> ThreadSafeRNG::uniform_dist_(0.0, 1.0);
thread_local bool ThreadSafeRNG::initialized_ = false;

void ThreadSafeRNG::initialize_thread_local() {
    if (!initialized_) {
        generator_.seed(std::random_device{}());
        initialized_ = true;
    }
}

double ThreadSafeRNG::normal() {
    initialize_thread_local();
    return normal_dist_(generator_);
}

double ThreadSafeRNG::uniform() {
    initialize_thread_local();
    return uniform_dist_(generator_);
}

std::pair<double, double> ThreadSafeRNG::correlated_normal(double correlation) {
    initialize_thread_local();
    
    double z1 = normal_dist_(generator_);
    double z2 = normal_dist_(generator_);
    
    double w1 = z1;
    double w2 = correlation * z1 + std::sqrt(1.0 - correlation * correlation) * z2;
    
    return {w1, w2};
}

void ThreadSafeRNG::set_seed(unsigned int seed) {
    generator_.seed(seed);
    initialized_ = true;
}

// BaseParameters implementation
bool BaseParameters::is_valid() const noexcept {
    return spot_price > 0.0 && 
           strike_price > 0.0 && 
           time_to_expiry > 0.0 && 
           risk_free_rate >= 0.0 && 
           volatility > 0.0 && 
           dividend_yield >= 0.0;
}

std::string BaseParameters::validation_error() const noexcept {
    std::ostringstream oss;
    
    if (spot_price <= 0.0) oss << "Spot price must be positive. ";
    if (strike_price <= 0.0) oss << "Strike price must be positive. ";
    if (time_to_expiry <= 0.0) oss << "Time to expiry must be positive. ";
    if (risk_free_rate < 0.0) oss << "Risk-free rate cannot be negative. ";
    if (volatility <= 0.0) oss << "Volatility must be positive. ";
    if (dividend_yield < 0.0) oss << "Dividend yield cannot be negative. ";
    
    return oss.str();
}

// ExoticOption base class implementation
std::vector<double> ExoticOption::generate_gbm_path(int num_steps) const {
    std::vector<double> path;
    path.reserve(num_steps + 1);
    path.push_back(params_.spot_price);
    
    const double dt = params_.time_to_expiry / num_steps;
    const double drift = (params_.risk_free_rate - params_.dividend_yield - 0.5 * params_.volatility * params_.volatility) * dt;
    const double diffusion = params_.volatility * std::sqrt(dt);
    
    double current_price = params_.spot_price;
    
    for (int i = 0; i < num_steps; ++i) {
        double dW = ThreadSafeRNG::normal();
        current_price *= std::exp(drift + diffusion * dW);
        path.push_back(current_price);
    }
    
    return path;
}

PricingResult ExoticOption::run_monte_carlo(int num_simulations, int num_steps) const {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Use config defaults if not specified
    if (num_simulations == 0) {
        num_simulations = Config::ConfigManager::getInstance().getMonteCarloSimulations();
    }
    if (num_steps == 0) {
        num_steps = Config::ConfigManager::getInstance().getMonteCarloSteps();
    }
    
    logger_.info("Starting Monte Carlo simulation: {} simulations, {} steps", 
                 num_simulations, num_steps);
    
    PricingResult result;
    std::vector<double> payoffs;
    payoffs.reserve(num_simulations);
    
    try {
        for (int sim = 0; sim < num_simulations; ++sim) {
            auto path = generate_gbm_path(num_steps);
            double payoff = calculate_payoff(path);
            payoffs.push_back(payoff);
        }
        
        // Calculate statistics
        double sum = std::accumulate(payoffs.begin(), payoffs.end(), 0.0);
        result.price = std::exp(-params_.risk_free_rate * params_.time_to_expiry) * sum / num_simulations;
        
        // Calculate standard error
        double sum_squared = 0.0;
        for (double payoff : payoffs) {
            double discounted_payoff = std::exp(-params_.risk_free_rate * params_.time_to_expiry) * payoff;
            sum_squared += (discounted_payoff - result.price) * (discounted_payoff - result.price);
        }
        
        result.standard_error = std::sqrt(sum_squared / (num_simulations * (num_simulations - 1)));
        result.confidence_interval = 1.96 * result.standard_error;  // 95% CI
        result.simulations_used = num_simulations;
        result.is_valid = true;
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        result.execution_time_ms = duration.count() / 1000.0;
        
        logger_.info("Monte Carlo completed: price=${:.4f} ± ${:.4f} ({:.2f}ms)", 
                     result.price, result.confidence_interval, result.execution_time_ms);
        
    } catch (const std::exception& e) {
        result.is_valid = false;
        result.error_message = e.what();
        logger_.error("Monte Carlo simulation failed: {}", e.what());
    }
    
    return result;
}

// AsianOption implementation
PricingResult AsianOption::price(int num_simulations, int num_steps) {
    logger_.info("Pricing {} Asian option", is_call_ ? "call" : "put");
    return run_monte_carlo(num_simulations, num_steps);
}

std::string AsianOption::get_name() const {
    std::ostringstream oss;
    oss << (arithmetic_average_ ? "Arithmetic" : "Geometric") << " Asian "
        << (is_call_ ? "Call" : "Put") << " Option";
    return oss.str();
}

double AsianOption::calculate_payoff(const std::vector<double>& path) const {
    double average_price;
    
    if (arithmetic_average_) {
        // Arithmetic average
        double sum = std::accumulate(path.begin(), path.end(), 0.0);
        average_price = sum / path.size();
    } else {
        // Geometric average
        double log_sum = 0.0;
        for (double price : path) {
            log_sum += std::log(price);
        }
        average_price = std::exp(log_sum / path.size());
    }
    
    if (is_call_) {
        return std::max(average_price - params_.strike_price, 0.0);
    } else {
        return std::max(params_.strike_price - average_price, 0.0);
    }
}

// BarrierOption implementation
PricingResult BarrierOption::price(int num_simulations, int num_steps) {
    logger_.info("Pricing {} barrier option with barrier at {}", 
                 Utils::barrier_type_to_string(barrier_type_), barrier_level_);
    return run_monte_carlo(num_simulations, num_steps);
}

std::string BarrierOption::get_name() const {
    std::ostringstream oss;
    oss << Utils::barrier_type_to_string(barrier_type_) << " "
        << (is_call_ ? "Call" : "Put") << " Barrier Option";
    return oss.str();
}

double BarrierOption::calculate_payoff(const std::vector<double>& path) const {
    bool barrier_hit = check_barrier_hit(path);
    double final_price = path.back();
    
    bool should_payout = false;
    double payoff = 0.0;
    
    switch (barrier_type_) {
        case BarrierType::UP_AND_OUT:
        case BarrierType::DOWN_AND_OUT:
            // Knock-out: pay if barrier NOT hit
            should_payout = !barrier_hit;
            if (!should_payout && rebate_ > 0.0) {
                payoff = rebate_;  // Pay rebate if knocked out
            }
            break;
            
        case BarrierType::UP_AND_IN:
        case BarrierType::DOWN_AND_IN:
            // Knock-in: pay if barrier hit
            should_payout = barrier_hit;
            break;
    }
    
    if (should_payout) {
        if (is_call_) {
            payoff = std::max(final_price - params_.strike_price, 0.0);
        } else {
            payoff = std::max(params_.strike_price - final_price, 0.0);
        }
    }
    
    return payoff;
}

bool BarrierOption::check_barrier_hit(const std::vector<double>& path) const {
    for (double price : path) {
        switch (barrier_type_) {
            case BarrierType::UP_AND_OUT:
            case BarrierType::UP_AND_IN:
                if (price >= barrier_level_) return true;
                break;
                
            case BarrierType::DOWN_AND_OUT:
            case BarrierType::DOWN_AND_IN:
                if (price <= barrier_level_) return true;
                break;
        }
    }
    return false;
}

// LookbackOption implementation
PricingResult LookbackOption::price(int num_simulations, int num_steps) {
    logger_.info("Pricing {} {} lookback option", 
                 fixed_strike_ ? "fixed strike" : "floating strike",
                 is_call_ ? "call" : "put");
    return run_monte_carlo(num_simulations, num_steps);
}

std::string LookbackOption::get_name() const {
    std::ostringstream oss;
    oss << (fixed_strike_ ? "Fixed Strike" : "Floating Strike") << " Lookback "
        << (is_call_ ? "Call" : "Put") << " Option";
    return oss.str();
}

double LookbackOption::calculate_payoff(const std::vector<double>& path) const {
    double max_price = *std::max_element(path.begin(), path.end());
    double min_price = *std::min_element(path.begin(), path.end());
    double final_price = path.back();
    
    if (fixed_strike_) {
        // Fixed strike lookback
        if (is_call_) {
            return std::max(max_price - params_.strike_price, 0.0);
        } else {
            return std::max(params_.strike_price - min_price, 0.0);
        }
    } else {
        // Floating strike lookback
        if (is_call_) {
            return final_price - min_price;  // Always non-negative
        } else {
            return max_price - final_price;  // Always non-negative
        }
    }
}

// DigitalOption implementation
PricingResult DigitalOption::price(int num_simulations, int num_steps) {
    logger_.info("Pricing {} digital option with cash amount {}", 
                 is_call_ ? "call" : "put", cash_amount_);
    return run_monte_carlo(num_simulations, num_steps);
}

std::string DigitalOption::get_name() const {
    std::ostringstream oss;
    oss << "Digital " << (is_call_ ? "Call" : "Put") << " Option";
    return oss.str();
}

double DigitalOption::calculate_payoff(const std::vector<double>& path) const {
    double final_price = path.back();
    
    bool in_the_money = is_call_ ? (final_price > params_.strike_price) 
                                 : (final_price < params_.strike_price);
    
    return in_the_money ? cash_amount_ : 0.0;
}

// MultiAssetParameters implementation
bool MultiAssetParameters::is_valid() const noexcept {
    if (spot_prices.empty() || volatilities.empty()) return false;
    if (spot_prices.size() != volatilities.size()) return false;
    if (correlation_matrix.size() != spot_prices.size()) return false;
    
    // Check positive prices and volatilities
    for (double price : spot_prices) {
        if (price <= 0.0) return false;
    }
    for (double vol : volatilities) {
        if (vol <= 0.0) return false;
    }
    
    // Check correlation matrix
    if (!Utils::is_valid_correlation_matrix(correlation_matrix)) return false;
    
    // Check other parameters
    return risk_free_rate >= 0.0 && dividend_yield >= 0.0 && time_to_expiry > 0.0;
}

std::string MultiAssetParameters::validation_error() const noexcept {
    std::ostringstream oss;
    
    if (spot_prices.empty()) oss << "No spot prices provided. ";
    if (volatilities.empty()) oss << "No volatilities provided. ";
    if (spot_prices.size() != volatilities.size()) oss << "Spot prices and volatilities size mismatch. ";
    if (correlation_matrix.size() != spot_prices.size()) oss << "Correlation matrix size mismatch. ";
    
    for (size_t i = 0; i < spot_prices.size(); ++i) {
        if (spot_prices[i] <= 0.0) {
            oss << "Spot price " << i << " must be positive. ";
        }
    }
    
    for (size_t i = 0; i < volatilities.size(); ++i) {
        if (volatilities[i] <= 0.0) {
            oss << "Volatility " << i << " must be positive. ";
        }
    }
    
    if (!Utils::is_valid_correlation_matrix(correlation_matrix)) {
        oss << "Invalid correlation matrix. ";
    }
    
    if (risk_free_rate < 0.0) oss << "Risk-free rate cannot be negative. ";
    if (dividend_yield < 0.0) oss << "Dividend yield cannot be negative. ";
    if (time_to_expiry <= 0.0) oss << "Time to expiry must be positive. ";
    
    return oss.str();
}

// RainbowOption implementation
PricingResult RainbowOption::price(int num_simulations, int num_steps) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Use config defaults if not specified
    if (num_simulations == 0) {
        num_simulations = Config::ConfigManager::getInstance().getMonteCarloSimulations();
    }
    if (num_steps == 0) {
        num_steps = Config::ConfigManager::getInstance().getMonteCarloSteps();
    }
    
    logger_.info("Pricing {} rainbow option on {} assets", 
                 is_best_of_ ? "best-of" : "worst-of", params_.num_assets());
    
    PricingResult result;
    std::vector<double> payoffs;
    payoffs.reserve(num_simulations);
    
    try {
        for (int sim = 0; sim < num_simulations; ++sim) {
            auto paths = generate_correlated_paths(num_steps);
            double payoff = calculate_payoff(paths);
            payoffs.push_back(payoff);
        }
        
        // Calculate statistics
        double sum = std::accumulate(payoffs.begin(), payoffs.end(), 0.0);
        result.price = std::exp(-params_.risk_free_rate * params_.time_to_expiry) * sum / num_simulations;
        
        // Calculate standard error
        double sum_squared = 0.0;
        for (double payoff : payoffs) {
            double discounted_payoff = std::exp(-params_.risk_free_rate * params_.time_to_expiry) * payoff;
            sum_squared += (discounted_payoff - result.price) * (discounted_payoff - result.price);
        }
        
        result.standard_error = std::sqrt(sum_squared / (num_simulations * (num_simulations - 1)));
        result.confidence_interval = 1.96 * result.standard_error;
        result.simulations_used = num_simulations;
        result.is_valid = true;
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        result.execution_time_ms = duration.count() / 1000.0;
        
        logger_.info("Rainbow option priced: ${:.4f} ± ${:.4f} ({:.2f}ms)", 
                     result.price, result.confidence_interval, result.execution_time_ms);
        
    } catch (const std::exception& e) {
        result.is_valid = false;
        result.error_message = e.what();
        logger_.error("Rainbow option pricing failed: {}", e.what());
    }
    
    return result;
}

std::string RainbowOption::get_name() const {
    std::ostringstream oss;
    oss << (is_best_of_ ? "Best of " : "Worst of ") << params_.num_assets()
        << " Rainbow " << (is_call_ ? "Call" : "Put") << " Option";
    return oss.str();
}

double RainbowOption::calculate_payoff(const std::vector<std::vector<double>>& paths) const {
    std::vector<double> final_prices;
    final_prices.reserve(paths.size());
    
    for (const auto& path : paths) {
        final_prices.push_back(path.back());
    }
    
    double selected_price;
    if (is_best_of_) {
        selected_price = *std::max_element(final_prices.begin(), final_prices.end());
    } else {
        selected_price = *std::min_element(final_prices.begin(), final_prices.end());
    }
    
    if (is_call_) {
        return std::max(selected_price - strike_price_, 0.0);
    } else {
        return std::max(strike_price_ - selected_price, 0.0);
    }
}

std::vector<std::vector<double>> RainbowOption::generate_correlated_paths(int num_steps) const {
    const size_t num_assets = params_.num_assets();
    std::vector<std::vector<double>> paths(num_assets);
    
    // Initialize paths with spot prices
    for (size_t i = 0; i < num_assets; ++i) {
        paths[i].reserve(num_steps + 1);
        paths[i].push_back(params_.spot_prices[i]);
    }
    
    // Cholesky decomposition for correlation
    auto chol_matrix = cholesky_decomposition(params_.correlation_matrix);
    
    const double dt = params_.time_to_expiry / num_steps;
    
    for (int step = 0; step < num_steps; ++step) {
        // Generate independent normal random numbers
        std::vector<double> independent_normals(num_assets);
        for (size_t i = 0; i < num_assets; ++i) {
            independent_normals[i] = ThreadSafeRNG::normal();
        }
        
        // Apply Cholesky decomposition to get correlated normals
        std::vector<double> correlated_normals(num_assets, 0.0);
        for (size_t i = 0; i < num_assets; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                correlated_normals[i] += chol_matrix[i * num_assets + j] * independent_normals[j];
            }
        }
        
        // Update each asset price
        for (size_t i = 0; i < num_assets; ++i) {
            double current_price = paths[i].back();
            double drift = (params_.risk_free_rate - params_.dividend_yield - 
                           0.5 * params_.volatilities[i] * params_.volatilities[i]) * dt;
            double diffusion = params_.volatilities[i] * std::sqrt(dt) * correlated_normals[i];
            
            double new_price = current_price * std::exp(drift + diffusion);
            paths[i].push_back(new_price);
        }
    }
    
    return paths;
}

std::vector<double> RainbowOption::cholesky_decomposition(const std::vector<std::vector<double>>& matrix) const {
    const size_t n = matrix.size();
    std::vector<double> L(n * n, 0.0);
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            if (i == j) {
                double sum = 0.0;
                for (size_t k = 0; k < j; ++k) {
                    sum += L[i * n + k] * L[i * n + k];
                }
                L[i * n + j] = std::sqrt(matrix[i][j] - sum);
            } else {
                double sum = 0.0;
                for (size_t k = 0; k < j; ++k) {
                    sum += L[i * n + k] * L[j * n + k];
                }
                L[i * n + j] = (matrix[i][j] - sum) / L[j * n + j];
            }
        }
    }
    
    return L;
}

// ExoticPortfolio implementation
void ExoticPortfolio::add_position(std::unique_ptr<ExoticOption> option, double quantity, 
                                  const std::string& name) {
    std::lock_guard<std::mutex> lock(portfolio_mutex_);
    
    std::string position_name = name.empty() ? generate_position_name(*option) : name;
    positions_.emplace_back(std::move(option), quantity, position_name);
    
    logger_.info("Added position: {} x{}", position_name, quantity);
}

double ExoticPortfolio::calculate_total_value(int num_simulations) {
    std::lock_guard<std::mutex> lock(portfolio_mutex_);
    
    double total_value = 0.0;
    
    logger_.info("Calculating portfolio value for {} positions", positions_.size());
    
    for (const auto& position : positions_) {
        auto result = position.option->price(num_simulations);
        if (result.is_valid) {
            double position_value = position.quantity * result.price;
            total_value += position_value;
            
            logger_.info("Position {}: ${:.4f} x {} = ${:.4f}", 
                        position.name, result.price, position.quantity, position_value);
        } else {
            logger_.error("Failed to price position {}: {}", position.name, result.error_message);
        }
    }
    
    logger_.info("Total portfolio value: ${:.4f}", total_value);
    return total_value;
}

size_t ExoticPortfolio::size() const {
    std::lock_guard<std::mutex> lock(portfolio_mutex_);
    return positions_.size();
}

void ExoticPortfolio::print_summary() const {
    std::lock_guard<std::mutex> lock(portfolio_mutex_);
    
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "                    EXOTIC OPTIONS PORTFOLIO" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    
    for (size_t i = 0; i < positions_.size(); ++i) {
        const auto& position = positions_[i];
        std::cout << std::left << std::setw(40) << position.name
                  << " | Quantity: " << std::setw(8) << std::fixed << std::setprecision(2) 
                  << position.quantity << std::endl;
    }
    
    std::cout << std::string(80, '=') << std::endl;
}

void ExoticPortfolio::clear() {
    std::lock_guard<std::mutex> lock(portfolio_mutex_);
    positions_.clear();
    logger_.info("Portfolio cleared");
}

std::string ExoticPortfolio::generate_position_name(const ExoticOption& option) const {
    return option.get_name();
}

// Utility functions
namespace Utils {

std::string barrier_type_to_string(BarrierType type) {
    switch (type) {
        case BarrierType::UP_AND_OUT:   return "Up-and-Out";
        case BarrierType::UP_AND_IN:    return "Up-and-In";
        case BarrierType::DOWN_AND_OUT: return "Down-and-Out";
        case BarrierType::DOWN_AND_IN:  return "Down-and-In";
        default:                        return "Unknown";
    }
}

double geometric_asian_price(const BaseParameters& params, bool is_call) {
    // Adjusted parameters for geometric Asian option
    double adj_vol = params.volatility / std::sqrt(3.0);
    double adj_rate = 0.5 * (params.risk_free_rate - params.dividend_yield + params.volatility * params.volatility / 6.0);
    
    // Use Black-Scholes formula with adjusted parameters
    double d1 = (std::log(params.spot_price / params.strike_price) + 
                (adj_rate + 0.5 * adj_vol * adj_vol) * params.time_to_expiry) / 
               (adj_vol * std::sqrt(params.time_to_expiry));
    double d2 = d1 - adj_vol * std::sqrt(params.time_to_expiry);
    
    // Normal CDF approximation
    auto norm_cdf = [](double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); };
    
    if (is_call) {
        return params.spot_price * std::exp(-params.dividend_yield * params.time_to_expiry) * norm_cdf(d1) -
               params.strike_price * std::exp(-params.risk_free_rate * params.time_to_expiry) * norm_cdf(d2);
    } else {
        return params.strike_price * std::exp(-params.risk_free_rate * params.time_to_expiry) * norm_cdf(-d2) -
               params.spot_price * std::exp(-params.dividend_yield * params.time_to_expiry) * norm_cdf(-d1);
    }
}

bool is_valid_correlation_matrix(const std::vector<std::vector<double>>& matrix) {
    const size_t n = matrix.size();
    
    // Check if square matrix
    for (const auto& row : matrix) {
        if (row.size() != n) return false;
    }
    
    // Check diagonal elements are 1
    for (size_t i = 0; i < n; ++i) {
        if (std::abs(matrix[i][i] - 1.0) > 1e-10) return false;
    }
    
    // Check symmetry and bounds
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (std::abs(matrix[i][j] - matrix[j][i]) > 1e-10) return false;
            if (matrix[i][j] < -1.0 || matrix[i][j] > 1.0) return false;
        }
    }
    
    // Check positive semi-definiteness (simplified check)
    // For a more rigorous check, we would compute eigenvalues
    return true;
}

std::vector<double> generate_antithetic_variates(const std::vector<double>& normal_variates) {
    std::vector<double> antithetic;
    antithetic.reserve(normal_variates.size());
    
    for (double variate : normal_variates) {
        antithetic.push_back(-variate);
    }
    
    return antithetic;
}

} // namespace Utils

} // namespace ExoticOptions