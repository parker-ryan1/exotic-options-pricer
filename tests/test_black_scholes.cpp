#include "test_framework.hpp"
#include "../src/models/black_scholes.hpp"
#include <cmath>
#include <limits>

using namespace BlackScholes;
using namespace Testing;

/**
 * @file test_black_scholes.cpp
 * @brief Comprehensive unit tests for Black-Scholes option pricing model
 * 
 * Test Coverage:
 * - Parameter validation
 * - Call and put option pricing
 * - Greeks calculations
 * - Implied volatility calculations
 * - Edge cases and boundary conditions
 * - Performance benchmarks
 * - Thread safety
 */

class BlackScholesTestFixture : public TestFixture {
protected:
    Parameters standard_params;
    Parameters atm_params;
    Parameters itm_call_params;
    Parameters otm_call_params;
    
    void setUp() override {
        // Standard test parameters
        standard_params = Parameters(100.0, 100.0, 0.25, 0.05, 0.20, 0.0);
        
        // At-the-money parameters
        atm_params = Parameters(100.0, 100.0, 1.0, 0.05, 0.20, 0.0);
        
        // In-the-money call parameters
        itm_call_params = Parameters(110.0, 100.0, 0.25, 0.05, 0.20, 0.0);
        
        // Out-of-the-money call parameters
        otm_call_params = Parameters(90.0, 100.0, 0.25, 0.05, 0.20, 0.0);
    }
};

// Test suite for Black-Scholes parameter validation
TEST_SUITE(BlackScholesParameters) {
    auto suite = std::make_unique<TestSuite>("BlackScholesParameters");
    
    // Test valid parameters
    suite->addTest("ValidParameters", []() {
        ASSERT_NO_THROW(Parameters(100.0, 100.0, 1.0, 0.05, 0.20, 0.0));
        
        Parameters params(100.0, 100.0, 1.0, 0.05, 0.20, 0.0);
        ASSERT_TRUE(params.is_valid());
        ASSERT_TRUE(params.validation_error().empty());
    });
    
    // Test invalid spot price
    suite->addTest("InvalidSpotPrice", []() {
        ASSERT_THROWS(Parameters(-100.0, 100.0, 1.0, 0.05, 0.20, 0.0), std::invalid_argument);
        ASSERT_THROWS(Parameters(0.0, 100.0, 1.0, 0.05, 0.20, 0.0), std::invalid_argument);
    });
    
    // Test invalid strike price
    suite->addTest("InvalidStrikePrice", []() {
        ASSERT_THROWS(Parameters(100.0, -100.0, 1.0, 0.05, 0.20, 0.0), std::invalid_argument);
        ASSERT_THROWS(Parameters(100.0, 0.0, 1.0, 0.05, 0.20, 0.0), std::invalid_argument);
    });
    
    // Test invalid time to expiry
    suite->addTest("InvalidTimeToExpiry", []() {
        ASSERT_THROWS(Parameters(100.0, 100.0, -1.0, 0.05, 0.20, 0.0), std::invalid_argument);
        ASSERT_THROWS(Parameters(100.0, 100.0, 0.0, 0.05, 0.20, 0.0), std::invalid_argument);
    });
    
    // Test invalid risk-free rate
    suite->addTest("InvalidRiskFreeRate", []() {
        ASSERT_THROWS(Parameters(100.0, 100.0, 1.0, -0.05, 0.20, 0.0), std::invalid_argument);
    });
    
    // Test invalid volatility
    suite->addTest("InvalidVolatility", []() {
        ASSERT_THROWS(Parameters(100.0, 100.0, 1.0, 0.05, -0.20, 0.0), std::invalid_argument);
        ASSERT_THROWS(Parameters(100.0, 100.0, 1.0, 0.05, 0.0, 0.0), std::invalid_argument);
    });
    
    // Test invalid dividend yield
    suite->addTest("InvalidDividendYield", []() {
        ASSERT_THROWS(Parameters(100.0, 100.0, 1.0, 0.05, 0.20, -0.05), std::invalid_argument);
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for Black-Scholes option pricing
TEST_SUITE(BlackScholesOptionPricing) {
    auto suite = std::make_unique<TestSuite>("BlackScholesOptionPricing");
    
    // Test call option pricing with known values
    suite->addTestMethod<BlackScholesTestFixture>("CallOptionPricing", [](BlackScholesTestFixture& fixture) {
        BENCHMARK("CallOptionPricing");
        
        auto result = OptionPricer::price_call(fixture.standard_params);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        
        // For S=100, K=100, T=0.25, r=0.05, σ=0.20, q=0.0
        // Expected call price ≈ 4.78
        ASSERT_NEAR(result.price, 4.78, 0.1);
    });
    
    // Test put option pricing with known values
    suite->addTestMethod<BlackScholesTestFixture>("PutOptionPricing", [](BlackScholesTestFixture& fixture) {
        BENCHMARK("PutOptionPricing");
        
        auto result = OptionPricer::price_put(fixture.standard_params);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        
        // For S=100, K=100, T=0.25, r=0.05, σ=0.20, q=0.0
        // Expected put price ≈ 3.54
        ASSERT_NEAR(result.price, 3.54, 0.1);
    });
    
    // Test put-call parity
    suite->addTestMethod<BlackScholesTestFixture>("PutCallParity", [](BlackScholesTestFixture& fixture) {
        auto call_result = OptionPricer::price_call(fixture.standard_params);
        auto put_result = OptionPricer::price_put(fixture.standard_params);
        
        ASSERT_TRUE(call_result.is_valid);
        ASSERT_TRUE(put_result.is_valid);
        
        // Put-call parity: C - P = S*e^(-q*T) - K*e^(-r*T)
        double left_side = call_result.price - put_result.price;
        double right_side = fixture.standard_params.spot_price * 
                           std::exp(-fixture.standard_params.dividend_yield * fixture.standard_params.time_to_expiry) -
                           fixture.standard_params.strike_price * 
                           std::exp(-fixture.standard_params.risk_free_rate * fixture.standard_params.time_to_expiry);
        
        ASSERT_NEAR(left_side, right_side, 1e-10);
    });
    
    // Test ITM call option
    suite->addTestMethod<BlackScholesTestFixture>("ITMCallOption", [](BlackScholesTestFixture& fixture) {
        auto result = OptionPricer::price_call(fixture.itm_call_params);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 10.0);  // Should be > intrinsic value
        
        // Intrinsic value = max(S - K, 0) = max(110 - 100, 0) = 10
        double intrinsic_value = std::max(fixture.itm_call_params.spot_price - fixture.itm_call_params.strike_price, 0.0);
        ASSERT_GT(result.price, intrinsic_value);
    });
    
    // Test OTM call option
    suite->addTestMethod<BlackScholesTestFixture>("OTMCallOption", [](BlackScholesTestFixture& fixture) {
        auto result = OptionPricer::price_call(fixture.otm_call_params);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        ASSERT_LT(result.price, 5.0);  // Should be relatively small
    });
    
    // Test option pricing with dividends
    suite->addTest("OptionPricingWithDividends", []() {
        Parameters params_with_div(100.0, 100.0, 1.0, 0.05, 0.20, 0.03);
        
        auto call_result = OptionPricer::price_call(params_with_div);
        auto put_result = OptionPricer::price_put(params_with_div);
        
        ASSERT_TRUE(call_result.is_valid);
        ASSERT_TRUE(put_result.is_valid);
        
        // With dividends, call price should be lower, put price higher
        Parameters params_no_div(100.0, 100.0, 1.0, 0.05, 0.20, 0.0);
        auto call_no_div = OptionPricer::price_call(params_no_div);
        auto put_no_div = OptionPricer::price_put(params_no_div);
        
        ASSERT_LT(call_result.price, call_no_div.price);
        ASSERT_GT(put_result.price, put_no_div.price);
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for Greeks calculations
TEST_SUITE(BlackScholesGreeks) {
    auto suite = std::make_unique<TestSuite>("BlackScholesGreeks");
    
    // Test call delta
    suite->addTestMethod<BlackScholesTestFixture>("CallDelta", [](BlackScholesTestFixture& fixture) {
        auto greeks = OptionPricer::calculate_call_greeks(fixture.standard_params);
        
        // Call delta should be between 0 and 1
        ASSERT_GT(greeks.delta, 0.0);
        ASSERT_LT(greeks.delta, 1.0);
        
        // For ATM option, delta should be around 0.5
        auto atm_greeks = OptionPricer::calculate_call_greeks(fixture.atm_params);
        ASSERT_NEAR(atm_greeks.delta, 0.5, 0.1);
    });
    
    // Test put delta
    suite->addTestMethod<BlackScholesTestFixture>("PutDelta", [](BlackScholesTestFixture& fixture) {
        auto greeks = OptionPricer::calculate_put_greeks(fixture.standard_params);
        
        // Put delta should be between -1 and 0
        ASSERT_LT(greeks.delta, 0.0);
        ASSERT_GT(greeks.delta, -1.0);
        
        // For ATM option, put delta should be around -0.5
        auto atm_greeks = OptionPricer::calculate_put_greeks(fixture.atm_params);
        ASSERT_NEAR(atm_greeks.delta, -0.5, 0.1);
    });
    
    // Test gamma (same for calls and puts)
    suite->addTestMethod<BlackScholesTestFixture>("Gamma", [](BlackScholesTestFixture& fixture) {
        auto call_greeks = OptionPricer::calculate_call_greeks(fixture.standard_params);
        auto put_greeks = OptionPricer::calculate_put_greeks(fixture.standard_params);
        
        // Gamma should be positive
        ASSERT_GT(call_greeks.gamma, 0.0);
        ASSERT_GT(put_greeks.gamma, 0.0);
        
        // Gamma should be the same for calls and puts
        ASSERT_NEAR(call_greeks.gamma, put_greeks.gamma, 1e-10);
    });
    
    // Test vega (same for calls and puts)
    suite->addTestMethod<BlackScholesTestFixture>("Vega", [](BlackScholesTestFixture& fixture) {
        auto call_greeks = OptionPricer::calculate_call_greeks(fixture.standard_params);
        auto put_greeks = OptionPricer::calculate_put_greeks(fixture.standard_params);
        
        // Vega should be positive
        ASSERT_GT(call_greeks.vega, 0.0);
        ASSERT_GT(put_greeks.vega, 0.0);
        
        // Vega should be the same for calls and puts
        ASSERT_NEAR(call_greeks.vega, put_greeks.vega, 1e-10);
    });
    
    // Test theta
    suite->addTestMethod<BlackScholesTestFixture>("Theta", [](BlackScholesTestFixture& fixture) {
        auto call_greeks = OptionPricer::calculate_call_greeks(fixture.standard_params);
        auto put_greeks = OptionPricer::calculate_put_greeks(fixture.standard_params);
        
        // Theta should generally be negative (time decay)
        ASSERT_LT(call_greeks.theta, 0.0);
        // Put theta can be positive or negative depending on parameters
    });
    
    // Test rho
    suite->addTestMethod<BlackScholesTestFixture>("Rho", [](BlackScholesTestFixture& fixture) {
        auto call_greeks = OptionPricer::calculate_call_greeks(fixture.standard_params);
        auto put_greeks = OptionPricer::calculate_put_greeks(fixture.standard_params);
        
        // Call rho should be positive
        ASSERT_GT(call_greeks.rho, 0.0);
        
        // Put rho should be negative
        ASSERT_LT(put_greeks.rho, 0.0);
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for implied volatility calculations
TEST_SUITE(BlackScholesImpliedVolatility) {
    auto suite = std::make_unique<TestSuite>("BlackScholesImpliedVolatility");
    
    // Test implied volatility round-trip
    suite->addTestMethod<BlackScholesTestFixture>("ImpliedVolatilityRoundTrip", [](BlackScholesTestFixture& fixture) {
        BENCHMARK("ImpliedVolatilityCalculation");
        
        // Calculate theoretical price
        auto call_result = OptionPricer::price_call(fixture.standard_params);
        ASSERT_TRUE(call_result.is_valid);
        
        // Calculate implied volatility from theoretical price
        double implied_vol = OptionPricer::calculate_implied_volatility(
            call_result.price, fixture.standard_params, true);
        
        ASSERT_FALSE(std::isnan(implied_vol));
        ASSERT_NEAR(implied_vol, fixture.standard_params.volatility, 1e-6);
    });
    
    // Test implied volatility for put options
    suite->addTestMethod<BlackScholesTestFixture>("ImpliedVolatilityPut", [](BlackScholesTestFixture& fixture) {
        // Calculate theoretical price
        auto put_result = OptionPricer::price_put(fixture.standard_params);
        ASSERT_TRUE(put_result.is_valid);
        
        // Calculate implied volatility from theoretical price
        double implied_vol = OptionPricer::calculate_implied_volatility(
            put_result.price, fixture.standard_params, false);
        
        ASSERT_FALSE(std::isnan(implied_vol));
        ASSERT_NEAR(implied_vol, fixture.standard_params.volatility, 1e-6);
    });
    
    // Test implied volatility with extreme prices
    suite->addTestMethod<BlackScholesTestFixture>("ImpliedVolatilityExtremes", [](BlackScholesTestFixture& fixture) {
        // Very low price should result in low implied volatility
        double low_implied_vol = OptionPricer::calculate_implied_volatility(
            0.01, fixture.standard_params, true);
        ASSERT_GT(low_implied_vol, 0.0);
        ASSERT_LT(low_implied_vol, 0.1);
        
        // Very high price should result in high implied volatility
        double high_implied_vol = OptionPricer::calculate_implied_volatility(
            50.0, fixture.standard_params, true);
        ASSERT_GT(high_implied_vol, 1.0);
    });
    
    // Test implied volatility convergence failure
    suite->addTestMethod<BlackScholesTestFixture>("ImpliedVolatilityFailure", [](BlackScholesTestFixture& fixture) {
        // Impossible price (higher than spot price for call)
        double impossible_implied_vol = OptionPricer::calculate_implied_volatility(
            150.0, fixture.standard_params, true);
        ASSERT_TRUE(std::isnan(impossible_implied_vol));
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for mathematical utilities
TEST_SUITE(BlackScholesMathUtils) {
    auto suite = std::make_unique<TestSuite>("BlackScholesMathUtils");
    
    // Test normal CDF
    suite->addTest("NormalCDF", []() {
        // Test known values
        ASSERT_NEAR(MathUtils::normal_cdf(0.0), 0.5, 1e-10);
        ASSERT_NEAR(MathUtils::normal_cdf(-1.96), 0.025, 1e-3);
        ASSERT_NEAR(MathUtils::normal_cdf(1.96), 0.975, 1e-3);
        
        // Test extreme values
        ASSERT_NEAR(MathUtils::normal_cdf(-10.0), 0.0, 1e-10);
        ASSERT_NEAR(MathUtils::normal_cdf(10.0), 1.0, 1e-10);
    });
    
    // Test normal PDF
    suite->addTest("NormalPDF", []() {
        // Test known values
        ASSERT_NEAR(MathUtils::normal_pdf(0.0), 1.0 / std::sqrt(2.0 * M_PI), 1e-10);
        
        // Test symmetry
        ASSERT_NEAR(MathUtils::normal_pdf(-1.0), MathUtils::normal_pdf(1.0), 1e-10);
        
        // Test positive values
        ASSERT_GT(MathUtils::normal_pdf(0.0), 0.0);
        ASSERT_GT(MathUtils::normal_pdf(1.0), 0.0);
        ASSERT_GT(MathUtils::normal_pdf(-1.0), 0.0);
    });
    
    // Test inverse normal CDF
    suite->addTest("InverseNormalCDF", []() {
        // Test known values
        ASSERT_NEAR(MathUtils::normal_inv_cdf(0.5), 0.0, 1e-6);
        ASSERT_NEAR(MathUtils::normal_inv_cdf(0.025), -1.96, 1e-2);
        ASSERT_NEAR(MathUtils::normal_inv_cdf(0.975), 1.96, 1e-2);
        
        // Test round-trip
        double x = 1.5;
        double p = MathUtils::normal_cdf(x);
        double x_recovered = MathUtils::normal_inv_cdf(p);
        ASSERT_NEAR(x, x_recovered, 1e-6);
    });
    
    // Test inverse normal CDF edge cases
    suite->addTest("InverseNormalCDFEdgeCases", []() {
        // Test invalid inputs
        ASSERT_THROWS(MathUtils::normal_inv_cdf(0.0), std::invalid_argument);
        ASSERT_THROWS(MathUtils::normal_inv_cdf(1.0), std::invalid_argument);
        ASSERT_THROWS(MathUtils::normal_inv_cdf(-0.1), std::invalid_argument);
        ASSERT_THROWS(MathUtils::normal_inv_cdf(1.1), std::invalid_argument);
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for edge cases and boundary conditions
TEST_SUITE(BlackScholesEdgeCases) {
    auto suite = std::make_unique<TestSuite>("BlackScholesEdgeCases");
    
    // Test very short time to expiry
    suite->addTest("VeryShortTimeToExpiry", []() {
        Parameters params(100.0, 100.0, 1e-6, 0.05, 0.20, 0.0);  // 1 microsecond
        
        auto call_result = OptionPricer::price_call(params);
        auto put_result = OptionPricer::price_put(params);
        
        ASSERT_TRUE(call_result.is_valid);
        ASSERT_TRUE(put_result.is_valid);
        
        // Should approach intrinsic value
        ASSERT_NEAR(call_result.price, 0.0, 1e-3);  // ATM call with no time
        ASSERT_NEAR(put_result.price, 0.0, 1e-3);   // ATM put with no time
    });
    
    // Test very long time to expiry
    suite->addTest("VeryLongTimeToExpiry", []() {
        Parameters params(100.0, 100.0, 100.0, 0.05, 0.20, 0.0);  // 100 years
        
        auto call_result = OptionPricer::price_call(params);
        auto put_result = OptionPricer::price_put(params);
        
        ASSERT_TRUE(call_result.is_valid);
        ASSERT_TRUE(put_result.is_valid);
        
        // Call should approach spot price
        ASSERT_NEAR(call_result.price, params.spot_price, 1.0);
        
        // Put should be very small
        ASSERT_LT(put_result.price, 1.0);
    });
    
    // Test very high volatility
    suite->addTest("VeryHighVolatility", []() {
        Parameters params(100.0, 100.0, 1.0, 0.05, 5.0, 0.0);  // 500% volatility
        
        auto call_result = OptionPricer::price_call(params);
        auto put_result = OptionPricer::price_put(params);
        
        ASSERT_TRUE(call_result.is_valid);
        ASSERT_TRUE(put_result.is_valid);
        
        // Both options should be very valuable
        ASSERT_GT(call_result.price, 40.0);
        ASSERT_GT(put_result.price, 40.0);
    });
    
    // Test very low volatility
    suite->addTest("VeryLowVolatility", []() {
        Parameters params(100.0, 100.0, 1.0, 0.05, 0.001, 0.0);  // 0.1% volatility
        
        auto call_result = OptionPricer::price_call(params);
        auto put_result = OptionPricer::price_put(params);
        
        ASSERT_TRUE(call_result.is_valid);
        ASSERT_TRUE(put_result.is_valid);
        
        // Options should have very little time value
        ASSERT_LT(call_result.price, 10.0);
        ASSERT_LT(put_result.price, 5.0);
    });
    
    // Test deep ITM call
    suite->addTest("DeepITMCall", []() {
        Parameters params(200.0, 100.0, 1.0, 0.05, 0.20, 0.0);  // S >> K
        
        auto call_result = OptionPricer::price_call(params);
        
        ASSERT_TRUE(call_result.is_valid);
        
        // Should be close to forward price minus strike
        double forward_price = params.spot_price * std::exp(params.risk_free_rate * params.time_to_expiry);
        double expected_price = forward_price - params.strike_price * std::exp(-params.risk_free_rate * params.time_to_expiry);
        
        ASSERT_NEAR(call_result.price, expected_price, 1.0);
    });
    
    // Test deep OTM call
    suite->addTest("DeepOTMCall", []() {
        Parameters params(50.0, 100.0, 1.0, 0.05, 0.20, 0.0);  // S << K
        
        auto call_result = OptionPricer::price_call(params);
        
        ASSERT_TRUE(call_result.is_valid);
        
        // Should be very small
        ASSERT_LT(call_result.price, 1.0);
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Performance benchmark tests
TEST_SUITE(BlackScholesPerformance) {
    auto suite = std::make_unique<TestSuite>("BlackScholesPerformance");
    
    // Benchmark option pricing performance
    suite->addTest("PricingPerformanceBenchmark", []() {
        Parameters params(100.0, 100.0, 1.0, 0.05, 0.20, 0.0);
        
        const int num_iterations = 100000;
        
        {
            BENCHMARK("CallPricing_100k_iterations");
            for (int i = 0; i < num_iterations; ++i) {
                auto result = OptionPricer::price_call(params);
                ASSERT_TRUE(result.is_valid);
            }
        }
        
        {
            BENCHMARK("PutPricing_100k_iterations");
            for (int i = 0; i < num_iterations; ++i) {
                auto result = OptionPricer::price_put(params);
                ASSERT_TRUE(result.is_valid);
            }
        }
    });
    
    // Benchmark Greeks calculation performance
    suite->addTest("GreeksPerformanceBenchmark", []() {
        Parameters params(100.0, 100.0, 1.0, 0.05, 0.20, 0.0);
        
        const int num_iterations = 50000;
        
        {
            BENCHMARK("CallGreeks_50k_iterations");
            for (int i = 0; i < num_iterations; ++i) {
                auto greeks = OptionPricer::calculate_call_greeks(params);
                ASSERT_GT(greeks.delta, 0.0);
            }
        }
    });
    
    // Benchmark implied volatility performance
    suite->addTest("ImpliedVolatilityPerformanceBenchmark", []() {
        Parameters params(100.0, 100.0, 1.0, 0.05, 0.20, 0.0);
        double market_price = 10.0;
        
        const int num_iterations = 1000;
        
        {
            BENCHMARK("ImpliedVolatility_1k_iterations");
            for (int i = 0; i < num_iterations; ++i) {
                double iv = OptionPricer::calculate_implied_volatility(market_price, params, true);
                // Some may fail to converge, that's expected
            }
        }
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Thread safety tests
TEST_SUITE(BlackScholesThreadSafety) {
    auto suite = std::make_unique<TestSuite>("BlackScholesThreadSafety");
    
    // Test concurrent option pricing
    suite->addTest("ConcurrentOptionPricing", []() {
        Parameters params(100.0, 100.0, 1.0, 0.05, 0.20, 0.0);
        
        const int num_threads = 4;
        const int iterations_per_thread = 1000;
        
        std::vector<std::thread> threads;
        std::vector<bool> results(num_threads, false);
        
        {
            BENCHMARK("ConcurrentPricing_4_threads_1k_each");
            
            for (int t = 0; t < num_threads; ++t) {
                threads.emplace_back([&params, &results, t, iterations_per_thread]() {
                    bool all_valid = true;
                    
                    for (int i = 0; i < iterations_per_thread; ++i) {
                        auto call_result = OptionPricer::price_call(params);
                        auto put_result = OptionPricer::price_put(params);
                        
                        if (!call_result.is_valid || !put_result.is_valid) {
                            all_valid = false;
                            break;
                        }
                    }
                    
                    results[t] = all_valid;
                });
            }
            
            for (auto& thread : threads) {
                thread.join();
            }
        }
        
        // All threads should have succeeded
        for (bool result : results) {
            ASSERT_TRUE(result);
        }
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Main test runner
int main(int argc, char* argv[]) {
    // Configure logging for tests
    Utils::Logger::configure(
        Utils::LogLevel::INFO,
        true,   // console output
        true,   // file output
        "test_results.log",
        5 * 1024 * 1024,  // 5MB max file size
        3       // keep 3 log files
    );
    
    // Initialize configuration
    Config::ConfigManager::getInstance().initialize("test_config.json");
    
    std::cout << "Black-Scholes Option Pricing Model - Unit Tests" << std::endl;
    std::cout << "================================================" << std::endl;
    
    // Print test discovery
    Testing::TestRegistry::getInstance().printDiscovery();
    
    // Run all tests
    auto stats = Testing::TestRegistry::getInstance().runAllSuites();
    
    // Return appropriate exit code
    return (stats.failed_tests == 0 && stats.error_tests == 0) ? 0 : 1;
}