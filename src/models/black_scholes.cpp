#include "black_scholes.hpp"
#include <cmath>
#include <sstream>
#include <algorithm>

namespace BlackScholes {

// Thread-local logger initialization
thread_local Logger OptionPricer::logger_("BlackScholes::OptionPricer");

// Parameters implementation
Parameters::Parameters(double S, double K, double T, double r, double sigma, double q)
    : spot_price(S), strike_price(K), time_to_expiry(T), 
      risk_free_rate(r), volatility(sigma), dividend_yield(q) {
    
    if (!is_valid()) {
        logger_.error("Invalid Black-Scholes parameters: {}", validation_error());
        throw std::invalid_argument(validation_error());
    }
    
    logger_.debug("Created Black-Scholes parameters: S={}, K={}, T={}, r={}, σ={}, q={}", 
                  S, K, T, r, sigma, q);
}

bool Parameters::is_valid() const noexcept {
    return spot_price > 0.0 && 
           strike_price > 0.0 && 
           time_to_expiry > 0.0 && 
           risk_free_rate >= 0.0 && 
           volatility > 0.0 && 
           dividend_yield >= 0.0;
}

std::string Parameters::validation_error() const noexcept {
    std::ostringstream oss;
    
    if (spot_price <= 0.0) oss << "Spot price must be positive. ";
    if (strike_price <= 0.0) oss << "Strike price must be positive. ";
    if (time_to_expiry <= 0.0) oss << "Time to expiry must be positive. ";
    if (risk_free_rate < 0.0) oss << "Risk-free rate cannot be negative. ";
    if (volatility <= 0.0) oss << "Volatility must be positive. ";
    if (dividend_yield < 0.0) oss << "Dividend yield cannot be negative. ";
    
    return oss.str();
}

// OptionPricer implementation
double OptionPricer::calculate_d1(const Parameters& params) noexcept {
    const double numerator = std::log(params.spot_price / params.strike_price) + 
                           (params.risk_free_rate - params.dividend_yield + 
                            0.5 * params.volatility * params.volatility) * params.time_to_expiry;
    const double denominator = params.volatility * std::sqrt(params.time_to_expiry);
    
    return numerator / denominator;
}

double OptionPricer::calculate_d2(const Parameters& params, double d1) noexcept {
    return d1 - params.volatility * std::sqrt(params.time_to_expiry);
}

PricingResult OptionPricer::price_call(const Parameters& params) {
    logger_.debug("Pricing call option with S={}, K={}, T={}", 
                  params.spot_price, params.strike_price, params.time_to_expiry);
    
    PricingResult result;
    
    try {
        const double d1 = calculate_d1(params);
        const double d2 = calculate_d2(params, d1);
        
        const double N_d1 = MathUtils::normal_cdf(d1);
        const double N_d2 = MathUtils::normal_cdf(d2);
        
        const double discount_factor = std::exp(-params.risk_free_rate * params.time_to_expiry);
        const double dividend_factor = std::exp(-params.dividend_yield * params.time_to_expiry);
        
        result.price = params.spot_price * dividend_factor * N_d1 - 
                      params.strike_price * discount_factor * N_d2;
        
        result.greeks = calculate_call_greeks(params);
        result.is_valid = true;
        
        logger_.info("Call option priced successfully: ${:.4f}", result.price);
        
    } catch (const std::exception& e) {
        result.is_valid = false;
        result.error_msg = e.what();
        logger_.error("Failed to price call option: {}", e.what());
    }
    
    return result;
}

PricingResult OptionPricer::price_put(const Parameters& params) {
    logger_.debug("Pricing put option with S={}, K={}, T={}", 
                  params.spot_price, params.strike_price, params.time_to_expiry);
    
    PricingResult result;
    
    try {
        const double d1 = calculate_d1(params);
        const double d2 = calculate_d2(params, d1);
        
        const double N_minus_d1 = MathUtils::normal_cdf(-d1);
        const double N_minus_d2 = MathUtils::normal_cdf(-d2);
        
        const double discount_factor = std::exp(-params.risk_free_rate * params.time_to_expiry);
        const double dividend_factor = std::exp(-params.dividend_yield * params.time_to_expiry);
        
        result.price = params.strike_price * discount_factor * N_minus_d2 - 
                      params.spot_price * dividend_factor * N_minus_d1;
        
        result.greeks = calculate_put_greeks(params);
        result.is_valid = true;
        
        logger_.info("Put option priced successfully: ${:.4f}", result.price);
        
    } catch (const std::exception& e) {
        result.is_valid = false;
        result.error_msg = e.what();
        logger_.error("Failed to price put option: {}", e.what());
    }
    
    return result;
}

Greeks OptionPricer::calculate_call_greeks(const Parameters& params) {
    const double d1 = calculate_d1(params);
    const double d2 = calculate_d2(params, d1);
    
    const double N_d1 = MathUtils::normal_cdf(d1);
    const double N_d2 = MathUtils::normal_cdf(d2);
    const double phi_d1 = MathUtils::normal_pdf(d1);
    
    const double sqrt_T = std::sqrt(params.time_to_expiry);
    const double discount_factor = std::exp(-params.risk_free_rate * params.time_to_expiry);
    const double dividend_factor = std::exp(-params.dividend_yield * params.time_to_expiry);
    
    Greeks greeks;
    
    // Delta: ∂V/∂S
    greeks.delta = dividend_factor * N_d1;
    
    // Gamma: ∂²V/∂S²
    greeks.gamma = dividend_factor * phi_d1 / (params.spot_price * params.volatility * sqrt_T);
    
    // Theta: ∂V/∂T (per day)
    const double theta_term1 = -params.spot_price * dividend_factor * phi_d1 * params.volatility / (2 * sqrt_T);
    const double theta_term2 = params.dividend_yield * params.spot_price * dividend_factor * N_d1;
    const double theta_term3 = -params.risk_free_rate * params.strike_price * discount_factor * N_d2;
    greeks.theta = (theta_term1 + theta_term2 + theta_term3) / 365.0;  // Per day
    
    // Vega: ∂V/∂σ (per 1%)
    greeks.vega = params.spot_price * dividend_factor * phi_d1 * sqrt_T / 100.0;
    
    // Rho: ∂V/∂r (per 1%)
    greeks.rho = params.strike_price * params.time_to_expiry * discount_factor * N_d2 / 100.0;
    
    return greeks;
}

Greeks OptionPricer::calculate_put_greeks(const Parameters& params) {
    const double d1 = calculate_d1(params);
    const double d2 = calculate_d2(params, d1);
    
    const double N_minus_d1 = MathUtils::normal_cdf(-d1);
    const double N_minus_d2 = MathUtils::normal_cdf(-d2);
    const double phi_d1 = MathUtils::normal_pdf(d1);
    
    const double sqrt_T = std::sqrt(params.time_to_expiry);
    const double discount_factor = std::exp(-params.risk_free_rate * params.time_to_expiry);
    const double dividend_factor = std::exp(-params.dividend_yield * params.time_to_expiry);
    
    Greeks greeks;
    
    // Delta: ∂V/∂S
    greeks.delta = -dividend_factor * N_minus_d1;
    
    // Gamma: ∂²V/∂S² (same as call)
    greeks.gamma = dividend_factor * phi_d1 / (params.spot_price * params.volatility * sqrt_T);
    
    // Theta: ∂V/∂T (per day)
    const double theta_term1 = -params.spot_price * dividend_factor * phi_d1 * params.volatility / (2 * sqrt_T);
    const double theta_term2 = -params.dividend_yield * params.spot_price * dividend_factor * N_minus_d1;
    const double theta_term3 = params.risk_free_rate * params.strike_price * discount_factor * N_minus_d2;
    greeks.theta = (theta_term1 + theta_term2 + theta_term3) / 365.0;  // Per day
    
    // Vega: ∂V/∂σ (per 1%, same as call)
    greeks.vega = params.spot_price * dividend_factor * phi_d1 * sqrt_T / 100.0;
    
    // Rho: ∂V/∂r (per 1%)
    greeks.rho = -params.strike_price * params.time_to_expiry * discount_factor * N_minus_d2 / 100.0;
    
    return greeks;
}

double OptionPricer::calculate_implied_volatility(
    double market_price,
    const Parameters& params,
    bool is_call,
    int max_iterations,
    double tolerance) {
    
    // Use config defaults if not specified
    if (max_iterations == 0) {
        max_iterations = Config::getInstance().getImpliedVolMaxIterations();
    }
    if (tolerance == 0.0) {
        tolerance = Config::getInstance().getImpliedVolTolerance();
    }
    
    logger_.debug("Calculating implied volatility for market price ${:.4f}", market_price);
    
    // Initial guess
    double vol = 0.2;
    
    for (int i = 0; i < max_iterations; ++i) {
        // Create temporary parameters with current volatility guess
        Parameters temp_params = params;
        temp_params.volatility = vol;
        
        // Calculate theoretical price and vega
        PricingResult result = is_call ? price_call(temp_params) : price_put(temp_params);
        
        if (!result.is_valid) {
            logger_.error("Failed to calculate theoretical price during IV calculation");
            return std::numeric_limits<double>::quiet_NaN();
        }
        
        const double price_diff = result.price - market_price;
        const double vega = result.greeks.vega * 100.0;  // Convert back to per unit
        
        // Check convergence
        if (std::abs(price_diff) < tolerance) {
            logger_.info("Implied volatility converged: {:.4f} after {} iterations", vol, i + 1);
            return vol;
        }
        
        // Check for zero vega
        if (std::abs(vega) < 1e-10) {
            logger_.warning("Zero vega encountered in IV calculation");
            break;
        }
        
        // Newton-Raphson update
        vol -= price_diff / vega;
        vol = std::max(vol, 0.001);  // Ensure positive volatility
    }
    
    logger_.warning("Implied volatility failed to converge after {} iterations", max_iterations);
    return std::numeric_limits<double>::quiet_NaN();
}

std::vector<std::string> OptionPricer::validate_assumptions(const Parameters& params) {
    std::vector<std::string> warnings;
    
    // Check for extreme parameters
    if (params.volatility > 2.0) {
        warnings.push_back("Very high volatility (>200%) may indicate model breakdown");
    }
    
    if (params.time_to_expiry > 10.0) {
        warnings.push_back("Very long time to expiry (>10 years) may reduce model accuracy");
    }
    
    if (params.risk_free_rate > 0.20) {
        warnings.push_back("Very high risk-free rate (>20%) is unusual");
    }
    
    // Check moneyness
    const double moneyness = params.spot_price / params.strike_price;
    if (moneyness < 0.5 || moneyness > 2.0) {
        warnings.push_back("Extreme moneyness may reduce model accuracy");
    }
    
    // Check for dividend yield vs risk-free rate
    if (params.dividend_yield > params.risk_free_rate + 0.10) {
        warnings.push_back("Dividend yield significantly higher than risk-free rate");
    }
    
    return warnings;
}

// MathUtils implementation
namespace MathUtils {

double normal_cdf(double x) noexcept {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

double normal_pdf(double x) noexcept {
    static const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

double normal_inv_cdf(double p) {
    if (p <= 0.0 || p >= 1.0) {
        throw std::invalid_argument("Probability must be in (0,1)");
    }
    
    // Beasley-Springer-Moro algorithm
    static const double a[4] = {
        2.50662823884,
        -18.61500062529,
        41.39119773534,
        -25.44106049637
    };
    
    static const double b[4] = {
        -8.47351093090,
        23.08336743743,
        -21.06224101826,
        3.13082909833
    };
    
    static const double c[9] = {
        0.3374754822726147,
        0.9761690190917186,
        0.1607979714918209,
        0.0276438810333863,
        0.0038405729373609,
        0.0003951896511919,
        0.0000321767881768,
        0.0000002888167364,
        0.0000003960315187
    };
    
    double x = p - 0.5;
    
    if (std::abs(x) < 0.42) {
        double r = x * x;
        return x * (((a[3] * r + a[2]) * r + a[1]) * r + a[0]) /
               ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
    }
    
    double r = p;
    if (x > 0.0) r = 1.0 - p;
    r = std::log(-std::log(r));
    
    double result = c[0];
    for (int i = 1; i < 9; ++i) {
        result += c[i] * std::pow(r, i);
    }
    
    if (x < 0.0) result = -result;
    
    return result;
}

} // namespace MathUtils

} // namespace BlackScholes