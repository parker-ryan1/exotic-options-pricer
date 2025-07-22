#pragma once

#include <vector>
#include <memory>
#include <random>
#include <string>
#include <functional>
#include <atomic>
#include <mutex>
#include "../utils/logger.hpp"
#include "../config/config.hpp"

/**
 * @file exotic_options.hpp
 * @brief Exotic options pricing models with Monte Carlo simulation
 * @author Quantitative Finance Team
 * @version 1.0
 * @date 2025
 * 
 * This file contains implementations of various exotic option types:
 * - Asian Options (arithmetic and geometric averaging)
 * - Barrier Options (knock-in/knock-out, up/down)
 * - Lookback Options (fixed and floating strike)
 * - Digital/Binary Options (cash-or-nothing)
 * - Rainbow Options (multi-asset derivatives)
 * 
 * All implementations are thread-safe and optimized for performance.
 */

namespace ExoticOptions {

/**
 * @brief Base parameters for all exotic options
 */
struct BaseParameters {
    double spot_price;      ///< Current underlying price (S > 0)
    double strike_price;    ///< Strike price (K > 0)
    double time_to_expiry;  ///< Time to expiration in years (T > 0)
    double risk_free_rate;  ///< Risk-free interest rate (r >= 0)
    double volatility;      ///< Annualized volatility (σ > 0)
    double dividend_yield;  ///< Continuous dividend yield (q >= 0)
    
    BaseParameters(double S = 100.0, double K = 100.0, double T = 1.0,
                   double r = 0.05, double sigma = 0.20, double q = 0.0)
        : spot_price(S), strike_price(K), time_to_expiry(T),
          risk_free_rate(r), volatility(sigma), dividend_yield(q) {}
    
    /**
     * @brief Validate parameters
     * @return true if all parameters are valid
     */
    bool is_valid() const noexcept;
    
    /**
     * @brief Get validation error message
     * @return Error message if invalid, empty if valid
     */
    std::string validation_error() const noexcept;
};

/**
 * @brief Pricing result with statistics
 */
struct PricingResult {
    double price;               ///< Option price
    double standard_error;      ///< Monte Carlo standard error
    double confidence_interval; ///< 95% confidence interval
    int simulations_used;       ///< Number of simulations performed
    double execution_time_ms;   ///< Execution time in milliseconds
    bool is_valid;              ///< Whether calculation was successful
    std::string error_message;  ///< Error message if calculation failed
    
    PricingResult() : price(0.0), standard_error(0.0), confidence_interval(0.0),
                     simulations_used(0), execution_time_ms(0.0), is_valid(false) {}
};

/**
 * @brief Thread-safe random number generator
 */
class ThreadSafeRNG {
private:
    thread_local static std::mt19937 generator_;
    thread_local static std::normal_distribution<double> normal_dist_;
    thread_local static std::uniform_real_distribution<double> uniform_dist_;
    thread_local static bool initialized_;
    
    static void initialize_thread_local();

public:
    /**
     * @brief Generate standard normal random number
     * @return N(0,1) random variable
     */
    static double normal();
    
    /**
     * @brief Generate uniform random number
     * @return U(0,1) random variable
     */
    static double uniform();
    
    /**
     * @brief Generate correlated normal random numbers
     * @param correlation Correlation coefficient (-1 <= ρ <= 1)
     * @return Pair of correlated N(0,1) random variables
     */
    static std::pair<double, double> correlated_normal(double correlation);
    
    /**
     * @brief Set random seed for current thread
     * @param seed Random seed
     */
    static void set_seed(unsigned int seed);
};

/**
 * @brief Base class for all exotic options
 */
class ExoticOption {
protected:
    BaseParameters params_;
    mutable Utils::Logger logger_;
    
public:
    explicit ExoticOption(const BaseParameters& params)
        : params_(params), logger_("ExoticOption") {
        if (!params_.is_valid()) {
            throw std::invalid_argument(params_.validation_error());
        }
    }
    
    virtual ~ExoticOption() = default;
    
    /**
     * @brief Price the option using Monte Carlo simulation
     * @param num_simulations Number of Monte Carlo simulations
     * @param num_steps Number of time steps per simulation
     * @return Pricing result with statistics
     */
    virtual PricingResult price(int num_simulations = 0, int num_steps = 0) = 0;
    
    /**
     * @brief Get option name/description
     * @return Option name
     */
    virtual std::string get_name() const = 0;
    
    /**
     * @brief Get option type (call/put)
     * @return true for call, false for put
     */
    virtual bool is_call() const = 0;
    
    /**
     * @brief Calculate payoff for given path
     * @param path Price path
     * @return Option payoff
     */
    virtual double calculate_payoff(const std::vector<double>& path) const = 0;
    
    /**
     * @brief Get underlying parameters
     * @return Base parameters
     */
    const BaseParameters& get_parameters() const { return params_; }

protected:
    /**
     * @brief Generate geometric Brownian motion path
     * @param num_steps Number of time steps
     * @return Price path
     */
    std::vector<double> generate_gbm_path(int num_steps) const;
    
    /**
     * @brief Run Monte Carlo simulation
     * @param num_simulations Number of simulations
     * @param num_steps Number of time steps
     * @return Pricing result
     */
    PricingResult run_monte_carlo(int num_simulations, int num_steps) const;
};

/**
 * @brief Asian option with arithmetic averaging
 */
class AsianOption : public ExoticOption {
private:
    bool is_call_;
    bool arithmetic_average_;  ///< true for arithmetic, false for geometric
    
public:
    AsianOption(const BaseParameters& params, bool is_call = true, 
                bool arithmetic_average = true)
        : ExoticOption(params), is_call_(is_call), arithmetic_average_(arithmetic_average) {}
    
    PricingResult price(int num_simulations = 0, int num_steps = 0) override;
    std::string get_name() const override;
    bool is_call() const override { return is_call_; }
    double calculate_payoff(const std::vector<double>& path) const override;
    
    /**
     * @brief Get averaging type
     * @return true for arithmetic, false for geometric
     */
    bool is_arithmetic() const { return arithmetic_average_; }
};

/**
 * @brief Barrier option types
 */
enum class BarrierType {
    UP_AND_OUT,     ///< Knock-out when price goes above barrier
    UP_AND_IN,      ///< Knock-in when price goes above barrier
    DOWN_AND_OUT,   ///< Knock-out when price goes below barrier
    DOWN_AND_IN     ///< Knock-in when price goes below barrier
};

/**
 * @brief Barrier option implementation
 */
class BarrierOption : public ExoticOption {
private:
    bool is_call_;
    BarrierType barrier_type_;
    double barrier_level_;
    double rebate_;  ///< Rebate paid if barrier is hit (for knock-out options)
    
public:
    BarrierOption(const BaseParameters& params, double barrier_level,
                  BarrierType barrier_type, bool is_call = true, double rebate = 0.0)
        : ExoticOption(params), is_call_(is_call), barrier_type_(barrier_type),
          barrier_level_(barrier_level), rebate_(rebate) {
        if (barrier_level <= 0.0) {
            throw std::invalid_argument("Barrier level must be positive");
        }
    }
    
    PricingResult price(int num_simulations = 0, int num_steps = 0) override;
    std::string get_name() const override;
    bool is_call() const override { return is_call_; }
    double calculate_payoff(const std::vector<double>& path) const override;
    
    /**
     * @brief Get barrier level
     * @return Barrier level
     */
    double get_barrier_level() const { return barrier_level_; }
    
    /**
     * @brief Get barrier type
     * @return Barrier type
     */
    BarrierType get_barrier_type() const { return barrier_type_; }

private:
    bool check_barrier_hit(const std::vector<double>& path) const;
};

/**
 * @brief Lookback option implementation
 */
class LookbackOption : public ExoticOption {
private:
    bool is_call_;
    bool fixed_strike_;  ///< true for fixed strike, false for floating strike
    
public:
    LookbackOption(const BaseParameters& params, bool is_call = true, 
                   bool fixed_strike = true)
        : ExoticOption(params), is_call_(is_call), fixed_strike_(fixed_strike) {}
    
    PricingResult price(int num_simulations = 0, int num_steps = 0) override;
    std::string get_name() const override;
    bool is_call() const override { return is_call_; }
    double calculate_payoff(const std::vector<double>& path) const override;
    
    /**
     * @brief Check if option has fixed strike
     * @return true for fixed strike, false for floating
     */
    bool is_fixed_strike() const { return fixed_strike_; }
};

/**
 * @brief Digital/Binary option implementation
 */
class DigitalOption : public ExoticOption {
private:
    bool is_call_;
    double cash_amount_;  ///< Cash amount paid if option finishes in-the-money
    
public:
    DigitalOption(const BaseParameters& params, bool is_call = true, 
                  double cash_amount = 1.0)
        : ExoticOption(params), is_call_(is_call), cash_amount_(cash_amount) {
        if (cash_amount <= 0.0) {
            throw std::invalid_argument("Cash amount must be positive");
        }
    }
    
    PricingResult price(int num_simulations = 0, int num_steps = 0) override;
    std::string get_name() const override;
    bool is_call() const override { return is_call_; }
    double calculate_payoff(const std::vector<double>& path) const override;
    
    /**
     * @brief Get cash amount
     * @return Cash amount paid if ITM
     */
    double get_cash_amount() const { return cash_amount_; }
};

/**
 * @brief Multi-asset parameters for rainbow options
 */
struct MultiAssetParameters {
    std::vector<double> spot_prices;     ///< Initial prices for each asset
    std::vector<double> volatilities;    ///< Volatilities for each asset
    std::vector<std::vector<double>> correlation_matrix;  ///< Correlation matrix
    double risk_free_rate;               ///< Risk-free rate
    double dividend_yield;               ///< Dividend yield (assumed same for all assets)
    double time_to_expiry;               ///< Time to expiration
    
    MultiAssetParameters(const std::vector<double>& spots,
                        const std::vector<double>& vols,
                        const std::vector<std::vector<double>>& corr_matrix,
                        double r = 0.05, double q = 0.0, double T = 1.0)
        : spot_prices(spots), volatilities(vols), correlation_matrix(corr_matrix),
          risk_free_rate(r), dividend_yield(q), time_to_expiry(T) {}
    
    bool is_valid() const noexcept;
    std::string validation_error() const noexcept;
    size_t num_assets() const { return spot_prices.size(); }
};

/**
 * @brief Rainbow option on multiple assets
 */
class RainbowOption {
private:
    MultiAssetParameters params_;
    double strike_price_;
    bool is_call_;
    bool is_best_of_;  ///< true for best-of, false for worst-of
    mutable Utils::Logger logger_;
    
public:
    RainbowOption(const MultiAssetParameters& params, double strike,
                  bool is_call = true, bool is_best_of = true)
        : params_(params), strike_price_(strike), is_call_(is_call),
          is_best_of_(is_best_of), logger_("RainbowOption") {
        if (!params_.is_valid()) {
            throw std::invalid_argument(params_.validation_error());
        }
        if (strike <= 0.0) {
            throw std::invalid_argument("Strike price must be positive");
        }
    }
    
    /**
     * @brief Price the rainbow option
     * @param num_simulations Number of Monte Carlo simulations
     * @param num_steps Number of time steps per simulation
     * @return Pricing result
     */
    PricingResult price(int num_simulations = 0, int num_steps = 0);
    
    std::string get_name() const;
    bool is_call() const { return is_call_; }
    bool is_best_of() const { return is_best_of_; }
    
    /**
     * @brief Calculate payoff for given multi-asset path
     * @param paths Price paths for all assets
     * @return Option payoff
     */
    double calculate_payoff(const std::vector<std::vector<double>>& paths) const;

private:
    std::vector<std::vector<double>> generate_correlated_paths(int num_steps) const;
    std::vector<double> cholesky_decomposition(const std::vector<std::vector<double>>& matrix) const;
};

/**
 * @brief Portfolio of exotic options
 */
class ExoticPortfolio {
private:
    struct PortfolioPosition {
        std::unique_ptr<ExoticOption> option;
        double quantity;
        std::string name;
        
        PortfolioPosition(std::unique_ptr<ExoticOption> opt, double qty, const std::string& n)
            : option(std::move(opt)), quantity(qty), name(n) {}
    };
    
    std::vector<PortfolioPosition> positions_;
    mutable Utils::Logger logger_;
    mutable std::mutex portfolio_mutex_;

public:
    ExoticPortfolio() : logger_("ExoticPortfolio") {}
    
    /**
     * @brief Add option to portfolio
     * @param option Exotic option to add
     * @param quantity Position size (positive for long, negative for short)
     * @param name Position name
     */
    void add_position(std::unique_ptr<ExoticOption> option, double quantity, 
                     const std::string& name = "");
    
    /**
     * @brief Calculate total portfolio value
     * @param num_simulations Number of simulations per option
     * @return Total portfolio value
     */
    double calculate_total_value(int num_simulations = 0);
    
    /**
     * @brief Get number of positions
     * @return Number of positions in portfolio
     */
    size_t size() const;
    
    /**
     * @brief Print portfolio summary
     */
    void print_summary() const;
    
    /**
     * @brief Clear all positions
     */
    void clear();

private:
    std::string generate_position_name(const ExoticOption& option) const;
};

/**
 * @brief Utility functions for exotic options
 */
namespace Utils {
    /**
     * @brief Convert barrier type to string
     * @param type Barrier type
     * @return String representation
     */
    std::string barrier_type_to_string(BarrierType type);
    
    /**
     * @brief Calculate theoretical Asian option price (geometric average)
     * @param params Base parameters
     * @param is_call true for call, false for put
     * @return Theoretical price
     */
    double geometric_asian_price(const BaseParameters& params, bool is_call);
    
    /**
     * @brief Validate correlation matrix
     * @param matrix Correlation matrix
     * @return true if valid correlation matrix
     */
    bool is_valid_correlation_matrix(const std::vector<std::vector<double>>& matrix);
    
    /**
     * @brief Generate antithetic variates for variance reduction
     * @param normal_variates Original normal random numbers
     * @return Antithetic variates
     */
    std::vector<double> generate_antithetic_variates(const std::vector<double>& normal_variates);
}

} // namespace ExoticOptions