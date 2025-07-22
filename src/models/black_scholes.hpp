#pragma once

#include <cmath>
#include <stdexcept>
#include "../utils/logger.hpp"
#include "../config/config.hpp"

/**
 * @file black_scholes.hpp
 * @brief Black-Scholes option pricing model implementation
 * @author Quantitative Finance Team
 * @version 1.0
 * @date 2025
 * 
 * This file contains the implementation of the Black-Scholes option pricing model
 * for European options. The model assumes constant volatility and risk-free rate.
 * 
 * Mathematical Model:
 * - Call Price: S*N(d1) - K*exp(-r*T)*N(d2)
 * - Put Price: K*exp(-r*T)*N(-d2) - S*N(-d1)
 * - d1 = (ln(S/K) + (r + σ²/2)*T) / (σ*√T)
 * - d2 = d1 - σ*√T
 * 
 * Where:
 * - S: Current stock price
 * - K: Strike price
 * - T: Time to expiration
 * - r: Risk-free rate
 * - σ: Volatility
 * - N(): Standard normal cumulative distribution function
 */

namespace BlackScholes {

/**
 * @brief Parameters for Black-Scholes option pricing
 * 
 * This structure holds all necessary parameters for Black-Scholes calculations.
 * All parameters are validated upon construction to ensure mathematical validity.
 */
struct Parameters {
    double spot_price;      ///< Current stock price (S > 0)
    double strike_price;    ///< Strike price (K > 0)
    double time_to_expiry;  ///< Time to expiration in years (T > 0)
    double risk_free_rate;  ///< Risk-free interest rate (r >= 0)
    double volatility;      ///< Annualized volatility (σ > 0)
    double dividend_yield;  ///< Continuous dividend yield (q >= 0)
    
    /**
     * @brief Construct Black-Scholes parameters with validation
     * @param S Current stock price
     * @param K Strike price
     * @param T Time to expiration in years
     * @param r Risk-free rate
     * @param sigma Volatility
     * @param q Dividend yield (default: 0.0)
     * @throws std::invalid_argument if any parameter is invalid
     */
    Parameters(double S, double K, double T, double r, double sigma, double q = 0.0);
    
    /**
     * @brief Validate all parameters
     * @return true if all parameters are valid
     */
    bool is_valid() const noexcept;
    
    /**
     * @brief Get parameter validation error message
     * @return Error message if invalid, empty string if valid
     */
    std::string validation_error() const noexcept;
};

/**
 * @brief Greeks (risk sensitivities) for options
 * 
 * The Greeks measure the sensitivity of option prices to various factors:
 * - Delta: Price sensitivity to underlying price changes
 * - Gamma: Delta sensitivity to underlying price changes
 * - Theta: Price sensitivity to time decay
 * - Vega: Price sensitivity to volatility changes
 * - Rho: Price sensitivity to interest rate changes
 */
struct Greeks {
    double delta;   ///< ∂V/∂S - Price sensitivity to spot price
    double gamma;   ///< ∂²V/∂S² - Delta sensitivity to spot price
    double theta;   ///< ∂V/∂T - Price sensitivity to time (per day)
    double vega;    ///< ∂V/∂σ - Price sensitivity to volatility (per 1%)
    double rho;     ///< ∂V/∂r - Price sensitivity to interest rate (per 1%)
    
    Greeks() = default;
    Greeks(double d, double g, double t, double v, double r) 
        : delta(d), gamma(g), theta(t), vega(v), rho(r) {}
};

/**
 * @brief Option pricing result containing price and Greeks
 */
struct PricingResult {
    double price;           ///< Option price
    Greeks greeks;          ///< Risk sensitivities
    double implied_vol;     ///< Implied volatility (if calculated)
    bool is_valid;          ///< Whether calculation was successful
    std::string error_msg;  ///< Error message if calculation failed
    
    PricingResult() : price(0.0), implied_vol(0.0), is_valid(false) {}
};

/**
 * @brief Black-Scholes option pricer class
 * 
 * Thread-safe implementation of Black-Scholes option pricing model.
 * Supports both call and put options with full Greeks calculation.
 */
class OptionPricer {
private:
    static thread_local Logger logger_;  ///< Thread-local logger instance
    
    /**
     * @brief Calculate d1 parameter for Black-Scholes formula
     * @param params Black-Scholes parameters
     * @return d1 value
     */
    static double calculate_d1(const Parameters& params) noexcept;
    
    /**
     * @brief Calculate d2 parameter for Black-Scholes formula
     * @param params Black-Scholes parameters
     * @param d1 Pre-calculated d1 value
     * @return d2 value
     */
    static double calculate_d2(const Parameters& params, double d1) noexcept;

public:
    /**
     * @brief Price a European call option
     * @param params Black-Scholes parameters
     * @return Pricing result with price and Greeks
     * @throws std::invalid_argument if parameters are invalid
     */
    static PricingResult price_call(const Parameters& params);
    
    /**
     * @brief Price a European put option
     * @param params Black-Scholes parameters
     * @return Pricing result with price and Greeks
     * @throws std::invalid_argument if parameters are invalid
     */
    static PricingResult price_put(const Parameters& params);
    
    /**
     * @brief Calculate Greeks for a call option
     * @param params Black-Scholes parameters
     * @return Greeks structure
     */
    static Greeks calculate_call_greeks(const Parameters& params);
    
    /**
     * @brief Calculate Greeks for a put option
     * @param params Black-Scholes parameters
     * @return Greeks structure
     */
    static Greeks calculate_put_greeks(const Parameters& params);
    
    /**
     * @brief Calculate implied volatility using Newton-Raphson method
     * @param market_price Observed market price
     * @param params Black-Scholes parameters (volatility will be ignored)
     * @param is_call true for call option, false for put
     * @param max_iterations Maximum number of iterations (default from config)
     * @param tolerance Convergence tolerance (default from config)
     * @return Implied volatility or NaN if convergence failed
     */
    static double calculate_implied_volatility(
        double market_price,
        const Parameters& params,
        bool is_call,
        int max_iterations = 0,  // 0 means use config default
        double tolerance = 0.0   // 0 means use config default
    );
    
    /**
     * @brief Validate Black-Scholes assumptions
     * @param params Parameters to validate
     * @return Validation warnings (empty if no issues)
     */
    static std::vector<std::string> validate_assumptions(const Parameters& params);
};

/**
 * @brief Utility functions for normal distribution
 */
namespace MathUtils {
    /**
     * @brief Standard normal cumulative distribution function
     * @param x Input value
     * @return N(x) - probability that standard normal variable ≤ x
     */
    double normal_cdf(double x) noexcept;
    
    /**
     * @brief Standard normal probability density function
     * @param x Input value
     * @return φ(x) - standard normal density at x
     */
    double normal_pdf(double x) noexcept;
    
    /**
     * @brief Inverse normal cumulative distribution function
     * @param p Probability (0 < p < 1)
     * @return x such that N(x) = p
     * @throws std::invalid_argument if p is not in (0,1)
     */
    double normal_inv_cdf(double p);
}

} // namespace BlackScholes