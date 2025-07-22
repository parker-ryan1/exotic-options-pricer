#include "test_framework.hpp"
#include "../src/models/exotic_options.hpp"
#include <cmath>
#include <limits>

using namespace ExoticOptions;
using namespace Testing;

/**
 * @file test_exotic_options.cpp
 * @brief Comprehensive unit tests for exotic options pricing models
 * 
 * Test Coverage:
 * - Parameter validation for all option types
 * - Asian options (arithmetic and geometric)
 * - Barrier options (all four types)
 * - Lookback options (fixed and floating strike)
 * - Digital options
 * - Rainbow options (multi-asset)
 * - Portfolio management
 * - Performance benchmarks
 * - Thread safety
 */

class ExoticOptionsTestFixture : public TestFixture {
protected:
    BaseParameters standard_params;
    BaseParameters high_vol_params;
    BaseParameters low_vol_params;
    
    void setUp() override {
        // Standard test parameters
        standard_params = BaseParameters(100.0, 100.0, 0.25, 0.05, 0.20, 0.0);
        
        // High volatility parameters
        high_vol_params = BaseParameters(100.0, 100.0, 0.25, 0.05, 0.50, 0.0);
        
        // Low volatility parameters
        low_vol_params = BaseParameters(100.0, 100.0, 0.25, 0.05, 0.05, 0.0);
    }
};

// Test suite for parameter validation
TEST_SUITE(ExoticOptionsParameters) {
    auto suite = std::make_unique<TestSuite>("ExoticOptionsParameters");
    
    // Test valid base parameters
    suite->addTest("ValidBaseParameters", []() {
        ASSERT_NO_THROW(BaseParameters(100.0, 100.0, 1.0, 0.05, 0.20, 0.0));
        
        BaseParameters params(100.0, 100.0, 1.0, 0.05, 0.20, 0.0);
        ASSERT_TRUE(params.is_valid());
        ASSERT_TRUE(params.validation_error().empty());
    });
    
    // Test invalid parameters
    suite->addTest("InvalidSpotPrice", []() {
        BaseParameters params(-100.0, 100.0, 1.0, 0.05, 0.20, 0.0);
        ASSERT_FALSE(params.is_valid());
        ASSERT_FALSE(params.validation_error().empty());
    });
    
    suite->addTest("InvalidStrikePrice", []() {
        BaseParameters params(100.0, -100.0, 1.0, 0.05, 0.20, 0.0);
        ASSERT_FALSE(params.is_valid());
    });
    
    suite->addTest("InvalidTimeToExpiry", []() {
        BaseParameters params(100.0, 100.0, -1.0, 0.05, 0.20, 0.0);
        ASSERT_FALSE(params.is_valid());
    });
    
    suite->addTest("InvalidVolatility", []() {
        BaseParameters params(100.0, 100.0, 1.0, 0.05, -0.20, 0.0);
        ASSERT_FALSE(params.is_valid());
    });
    
    // Test multi-asset parameters
    suite->addTest("ValidMultiAssetParameters", []() {
        std::vector<double> spots = {100.0, 110.0};
        std::vector<double> vols = {0.20, 0.25};
        std::vector<std::vector<double>> corr = {{1.0, 0.5}, {0.5, 1.0}};
        
        MultiAssetParameters params(spots, vols, corr);
        ASSERT_TRUE(params.is_valid());
    });
    
    suite->addTest("InvalidCorrelationMatrix", []() {
        std::vector<double> spots = {100.0, 110.0};
        std::vector<double> vols = {0.20, 0.25};
        std::vector<std::vector<double>> corr = {{1.0, 1.5}, {1.5, 1.0}};  // Invalid correlation > 1
        
        MultiAssetParameters params(spots, vols, corr);
        ASSERT_FALSE(params.is_valid());
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for Asian options
TEST_SUITE(AsianOptions) {
    auto suite = std::make_unique<TestSuite>("AsianOptions");
    
    // Test arithmetic Asian call pricing
    suite->addTestMethod<ExoticOptionsTestFixture>("ArithmeticAsianCall", [](ExoticOptionsTestFixture& fixture) {
        BENCHMARK("ArithmeticAsianCall");
        
        AsianOption option(fixture.standard_params, true, true);  // call, arithmetic
        auto result = option.price(10000, 100);  // Reduced simulations for testing
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        ASSERT_LT(result.price, fixture.standard_params.spot_price);  // Should be less than spot
        ASSERT_GT(result.standard_error, 0.0);
        ASSERT_EQ(result.simulations_used, 10000);
    });
    
    // Test arithmetic Asian put pricing
    suite->addTestMethod<ExoticOptionsTestFixture>("ArithmeticAsianPut", [](ExoticOptionsTestFixture& fixture) {
        AsianOption option(fixture.standard_params, false, true);  // put, arithmetic
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        ASSERT_LT(result.price, fixture.standard_params.strike_price);
    });
    
    // Test geometric Asian option
    suite->addTestMethod<ExoticOptionsTestFixture>("GeometricAsianCall", [](ExoticOptionsTestFixture& fixture) {
        AsianOption option(fixture.standard_params, true, false);  // call, geometric
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        
        // Geometric average should be less than arithmetic average
        AsianOption arith_option(fixture.standard_params, true, true);
        auto arith_result = arith_option.price(10000, 100);
        
        ASSERT_TRUE(arith_result.is_valid);
        ASSERT_LT(result.price, arith_result.price);
    });
    
    // Test Asian option with high volatility
    suite->addTestMethod<ExoticOptionsTestFixture>("AsianHighVolatility", [](ExoticOptionsTestFixture& fixture) {
        AsianOption option(fixture.high_vol_params, true, true);
        auto result = option.price(5000, 50);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        
        // Compare with standard volatility
        AsianOption std_option(fixture.standard_params, true, true);
        auto std_result = std_option.price(5000, 50);
        
        ASSERT_TRUE(std_result.is_valid);
        ASSERT_GT(result.price, std_result.price);  // Higher vol should give higher price
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for Barrier options
TEST_SUITE(BarrierOptions) {
    auto suite = std::make_unique<TestSuite>("BarrierOptions");
    
    // Test up-and-out barrier call
    suite->addTestMethod<ExoticOptionsTestFixture>("UpAndOutCall", [](ExoticOptionsTestFixture& fixture) {
        BENCHMARK("UpAndOutCall");
        
        double barrier = 120.0;  // Above current spot
        BarrierOption option(fixture.standard_params, barrier, BarrierType::UP_AND_OUT, true);
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GE(result.price, 0.0);
        ASSERT_EQ(option.get_barrier_level(), barrier);
        ASSERT_EQ(option.get_barrier_type(), BarrierType::UP_AND_OUT);
    });
    
    // Test down-and-in barrier put
    suite->addTestMethod<ExoticOptionsTestFixture>("DownAndInPut", [](ExoticOptionsTestFixture& fixture) {
        double barrier = 80.0;  // Below current spot
        BarrierOption option(fixture.standard_params, barrier, BarrierType::DOWN_AND_IN, false);
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GE(result.price, 0.0);
    });
    
    // Test barrier option with rebate
    suite->addTestMethod<ExoticOptionsTestFixture>("BarrierWithRebate", [](ExoticOptionsTestFixture& fixture) {
        double barrier = 120.0;
        double rebate = 5.0;
        BarrierOption option(fixture.standard_params, barrier, BarrierType::UP_AND_OUT, true, rebate);
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GE(result.price, 0.0);
        
        // Compare with no rebate
        BarrierOption no_rebate_option(fixture.standard_params, barrier, BarrierType::UP_AND_OUT, true, 0.0);
        auto no_rebate_result = no_rebate_option.price(10000, 100);
        
        ASSERT_TRUE(no_rebate_result.is_valid);
        ASSERT_GE(result.price, no_rebate_result.price);  // Rebate should increase value
    });
    
    // Test all barrier types
    suite->addTestMethod<ExoticOptionsTestFixture>("AllBarrierTypes", [](ExoticOptionsTestFixture& fixture) {
        std::vector<BarrierType> types = {
            BarrierType::UP_AND_OUT, BarrierType::UP_AND_IN,
            BarrierType::DOWN_AND_OUT, BarrierType::DOWN_AND_IN
        };
        
        for (auto type : types) {
            double barrier = (type == BarrierType::UP_AND_OUT || type == BarrierType::UP_AND_IN) ? 120.0 : 80.0;
            BarrierOption option(fixture.standard_params, barrier, type, true);
            auto result = option.price(5000, 50);
            
            ASSERT_TRUE(result.is_valid);
            ASSERT_GE(result.price, 0.0);
        }
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for Lookback options
TEST_SUITE(LookbackOptions) {
    auto suite = std::make_unique<TestSuite>("LookbackOptions");
    
    // Test fixed strike lookback call
    suite->addTestMethod<ExoticOptionsTestFixture>("FixedStrikeLookbackCall", [](ExoticOptionsTestFixture& fixture) {
        BENCHMARK("FixedStrikeLookbackCall");
        
        LookbackOption option(fixture.standard_params, true, true);  // call, fixed strike
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        ASSERT_TRUE(option.is_fixed_strike());
        
        // Fixed strike lookback call should be more valuable than regular call
        // (since it pays max(S_max - K, 0) instead of max(S_T - K, 0))
    });
    
    // Test floating strike lookback call
    suite->addTestMethod<ExoticOptionsTestFixture>("FloatingStrikeLookbackCall", [](ExoticOptionsTestFixture& fixture) {
        LookbackOption option(fixture.standard_params, true, false);  // call, floating strike
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        ASSERT_FALSE(option.is_fixed_strike());
        
        // Floating strike lookback call pays S_T - S_min, which is always >= 0
    });
    
    // Test lookback put options
    suite->addTestMethod<ExoticOptionsTestFixture>("LookbackPutOptions", [](ExoticOptionsTestFixture& fixture) {
        // Fixed strike put
        LookbackOption fixed_put(fixture.standard_params, false, true);
        auto fixed_result = fixed_put.price(5000, 50);
        
        ASSERT_TRUE(fixed_result.is_valid);
        ASSERT_GT(fixed_result.price, 0.0);
        
        // Floating strike put
        LookbackOption floating_put(fixture.standard_params, false, false);
        auto floating_result = floating_put.price(5000, 50);
        
        ASSERT_TRUE(floating_result.is_valid);
        ASSERT_GT(floating_result.price, 0.0);
    });
    
    // Test lookback option with high volatility
    suite->addTestMethod<ExoticOptionsTestFixture>("LookbackHighVolatility", [](ExoticOptionsTestFixture& fixture) {
        LookbackOption option(fixture.high_vol_params, true, true);
        auto result = option.price(5000, 50);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        
        // Compare with standard volatility
        LookbackOption std_option(fixture.standard_params, true, true);
        auto std_result = std_option.price(5000, 50);
        
        ASSERT_TRUE(std_result.is_valid);
        ASSERT_GT(result.price, std_result.price);  // Higher vol should give higher price
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for Digital options
TEST_SUITE(DigitalOptions) {
    auto suite = std::make_unique<TestSuite>("DigitalOptions");
    
    // Test digital call option
    suite->addTestMethod<ExoticOptionsTestFixture>("DigitalCall", [](ExoticOptionsTestFixture& fixture) {
        BENCHMARK("DigitalCall");
        
        double cash_amount = 10.0;
        DigitalOption option(fixture.standard_params, true, cash_amount);
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GE(result.price, 0.0);
        ASSERT_LE(result.price, cash_amount);  // Price should not exceed cash amount
        ASSERT_EQ(option.get_cash_amount(), cash_amount);
    });
    
    // Test digital put option
    suite->addTestMethod<ExoticOptionsTestFixture>("DigitalPut", [](ExoticOptionsTestFixture& fixture) {
        double cash_amount = 10.0;
        DigitalOption option(fixture.standard_params, false, cash_amount);
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GE(result.price, 0.0);
        ASSERT_LE(result.price, cash_amount);
    });
    
    // Test digital option with different cash amounts
    suite->addTestMethod<ExoticOptionsTestFixture>("DigitalDifferentCashAmounts", [](ExoticOptionsTestFixture& fixture) {
        std::vector<double> cash_amounts = {1.0, 5.0, 10.0, 20.0};
        
        for (double cash : cash_amounts) {
            DigitalOption option(fixture.standard_params, true, cash);
            auto result = option.price(5000, 50);
            
            ASSERT_TRUE(result.is_valid);
            ASSERT_GE(result.price, 0.0);
            ASSERT_LE(result.price, cash);
        }
    });
    
    // Test digital option pricing relationship
    suite->addTestMethod<ExoticOptionsTestFixture>("DigitalPricingRelationship", [](ExoticOptionsTestFixture& fixture) {
        double cash_amount = 10.0;
        
        // ATM digital call and put should sum to discounted cash amount
        DigitalOption call_option(fixture.standard_params, true, cash_amount);
        DigitalOption put_option(fixture.standard_params, false, cash_amount);
        
        auto call_result = call_option.price(10000, 100);
        auto put_result = put_option.price(10000, 100);
        
        ASSERT_TRUE(call_result.is_valid);
        ASSERT_TRUE(put_result.is_valid);
        
        double total_price = call_result.price + put_result.price;
        double discounted_cash = cash_amount * std::exp(-fixture.standard_params.risk_free_rate * 
                                                       fixture.standard_params.time_to_expiry);
        
        // Allow for Monte Carlo error
        ASSERT_NEAR(total_price, discounted_cash, 0.5);
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for Rainbow options
TEST_SUITE(RainbowOptions) {
    auto suite = std::make_unique<TestSuite>("RainbowOptions");
    
    // Test best-of-two rainbow call
    suite->addTest("BestOfTwoRainbowCall", []() {
        BENCHMARK("BestOfTwoRainbowCall");
        
        std::vector<double> spots = {100.0, 110.0};
        std::vector<double> vols = {0.20, 0.25};
        std::vector<std::vector<double>> corr = {{1.0, 0.5}, {0.5, 1.0}};
        
        MultiAssetParameters params(spots, vols, corr, 0.05, 0.0, 0.25);
        RainbowOption option(params, 100.0, true, true);  // call, best-of
        
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        ASSERT_TRUE(option.is_best_of());
    });
    
    // Test worst-of-two rainbow put
    suite->addTest("WorstOfTwoRainbowPut", []() {
        std::vector<double> spots = {100.0, 110.0};
        std::vector<double> vols = {0.20, 0.25};
        std::vector<std::vector<double>> corr = {{1.0, 0.3}, {0.3, 1.0}};
        
        MultiAssetParameters params(spots, vols, corr, 0.05, 0.0, 0.25);
        RainbowOption option(params, 105.0, false, false);  // put, worst-of
        
        auto result = option.price(10000, 100);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
        ASSERT_FALSE(option.is_best_of());
    });
    
    // Test correlation effect on rainbow options
    suite->addTest("RainbowCorrelationEffect", []() {
        std::vector<double> spots = {100.0, 100.0};
        std::vector<double> vols = {0.20, 0.20};
        
        // High correlation
        std::vector<std::vector<double>> high_corr = {{1.0, 0.9}, {0.9, 1.0}};
        MultiAssetParameters high_params(spots, vols, high_corr, 0.05, 0.0, 0.25);
        RainbowOption high_option(high_params, 100.0, true, true);
        
        // Low correlation
        std::vector<std::vector<double>> low_corr = {{1.0, 0.1}, {0.1, 1.0}};
        MultiAssetParameters low_params(spots, vols, low_corr, 0.05, 0.0, 0.25);
        RainbowOption low_option(low_params, 100.0, true, true);
        
        auto high_result = high_option.price(5000, 50);
        auto low_result = low_option.price(5000, 50);
        
        ASSERT_TRUE(high_result.is_valid);
        ASSERT_TRUE(low_result.is_valid);
        
        // For best-of options, lower correlation should give higher price
        ASSERT_GT(low_result.price, high_result.price);
    });
    
    // Test three-asset rainbow option
    suite->addTest("ThreeAssetRainbow", []() {
        std::vector<double> spots = {100.0, 110.0, 90.0};
        std::vector<double> vols = {0.20, 0.25, 0.18};
        std::vector<std::vector<double>> corr = {
            {1.0, 0.5, 0.3},
            {0.5, 1.0, 0.4},
            {0.3, 0.4, 1.0}
        };
        
        MultiAssetParameters params(spots, vols, corr, 0.05, 0.0, 0.25);
        RainbowOption option(params, 100.0, true, true);
        
        auto result = option.price(5000, 50);
        
        ASSERT_TRUE(result.is_valid);
        ASSERT_GT(result.price, 0.0);
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for Portfolio management
TEST_SUITE(ExoticPortfolio) {
    auto suite = std::make_unique<TestSuite>("ExoticPortfolio");
    
    // Test portfolio creation and management
    suite->addTestMethod<ExoticOptionsTestFixture>("PortfolioManagement", [](ExoticOptionsTestFixture& fixture) {
        BENCHMARK("PortfolioManagement");
        
        ExoticPortfolio portfolio;
        
        // Add Asian call
        auto asian_call = std::make_unique<AsianOption>(fixture.standard_params, true, true);
        portfolio.add_position(std::move(asian_call), 10.0, "Asian Call Long");
        
        // Add Barrier put
        auto barrier_put = std::make_unique<BarrierOption>(
            fixture.standard_params, 80.0, BarrierType::DOWN_AND_IN, false);
        portfolio.add_position(std::move(barrier_put), -5.0, "Barrier Put Short");
        
        // Add Digital call
        auto digital_call = std::make_unique<DigitalOption>(fixture.standard_params, true, 10.0);
        portfolio.add_position(std::move(digital_call), 20.0, "Digital Call Long");
        
        ASSERT_EQ(portfolio.size(), 3);
        
        // Calculate portfolio value
        double total_value = portfolio.calculate_total_value(5000);
        
        // Portfolio value should be finite
        ASSERT_TRUE(std::isfinite(total_value));
        
        // Print summary (for manual verification)
        portfolio.print_summary();
        
        // Clear portfolio
        portfolio.clear();
        ASSERT_EQ(portfolio.size(), 0);
    });
    
    // Test portfolio with mixed positions
    suite->addTestMethod<ExoticOptionsTestFixture>("MixedPortfolio", [](ExoticOptionsTestFixture& fixture) {
        ExoticPortfolio portfolio;
        
        // Long lookback call
        auto lookback_call = std::make_unique<LookbackOption>(fixture.standard_params, true, true);
        portfolio.add_position(std::move(lookback_call), 5.0);
        
        // Short lookback put
        auto lookback_put = std::make_unique<LookbackOption>(fixture.standard_params, false, true);
        portfolio.add_position(std::move(lookback_put), -3.0);
        
        ASSERT_EQ(portfolio.size(), 2);
        
        double total_value = portfolio.calculate_total_value(3000);
        ASSERT_TRUE(std::isfinite(total_value));
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Test suite for utility functions
TEST_SUITE(ExoticOptionsUtils) {
    auto suite = std::make_unique<TestSuite>("ExoticOptionsUtils");
    
    // Test barrier type string conversion
    suite->addTest("BarrierTypeToString", []() {
        ASSERT_EQ(Utils::barrier_type_to_string(BarrierType::UP_AND_OUT), "Up-and-Out");
        ASSERT_EQ(Utils::barrier_type_to_string(BarrierType::UP_AND_IN), "Up-and-In");
        ASSERT_EQ(Utils::barrier_type_to_string(BarrierType::DOWN_AND_OUT), "Down-and-Out");
        ASSERT_EQ(Utils::barrier_type_to_string(BarrierType::DOWN_AND_IN), "Down-and-In");
    });
    
    // Test correlation matrix validation
    suite->addTest("CorrelationMatrixValidation", []() {
        // Valid correlation matrix
        std::vector<std::vector<double>> valid_corr = {{1.0, 0.5}, {0.5, 1.0}};
        ASSERT_TRUE(Utils::is_valid_correlation_matrix(valid_corr));
        
        // Invalid correlation matrix (not symmetric)
        std::vector<std::vector<double>> invalid_corr = {{1.0, 0.5}, {0.3, 1.0}};
        ASSERT_FALSE(Utils::is_valid_correlation_matrix(invalid_corr));
        
        // Invalid correlation matrix (diagonal not 1)
        std::vector<std::vector<double>> invalid_diag = {{0.9, 0.5}, {0.5, 1.0}};
        ASSERT_FALSE(Utils::is_valid_correlation_matrix(invalid_diag));
    });
    
    // Test geometric Asian pricing
    suite->addTestMethod<ExoticOptionsTestFixture>("GeometricAsianPricing", [](ExoticOptionsTestFixture& fixture) {
        double call_price = Utils::geometric_asian_price(fixture.standard_params, true);
        double put_price = Utils::geometric_asian_price(fixture.standard_params, false);
        
        ASSERT_GT(call_price, 0.0);
        ASSERT_GT(put_price, 0.0);
        ASSERT_TRUE(std::isfinite(call_price));
        ASSERT_TRUE(std::isfinite(put_price));
    });
    
    // Test antithetic variates generation
    suite->addTest("AntitheticVariates", []() {
        std::vector<double> original = {1.0, -0.5, 2.0, -1.5};
        auto antithetic = Utils::generate_antithetic_variates(original);
        
        ASSERT_EQ(antithetic.size(), original.size());
        
        for (size_t i = 0; i < original.size(); ++i) {
            ASSERT_NEAR(antithetic[i], -original[i], 1e-10);
        }
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Performance benchmark tests
TEST_SUITE(ExoticOptionsPerformance) {
    auto suite = std::make_unique<TestSuite>("ExoticOptionsPerformance");
    
    // Benchmark Asian option pricing
    suite->addTestMethod<ExoticOptionsTestFixture>("AsianOptionBenchmark", [](ExoticOptionsTestFixture& fixture) {
        AsianOption option(fixture.standard_params, true, true);
        
        const int num_runs = 10;
        
        {
            BENCHMARK("AsianOption_10_runs_10k_sims");
            for (int i = 0; i < num_runs; ++i) {
                auto result = option.price(10000, 100);
                ASSERT_TRUE(result.is_valid);
            }
        }
    });
    
    // Benchmark barrier option pricing
    suite->addTestMethod<ExoticOptionsTestFixture>("BarrierOptionBenchmark", [](ExoticOptionsTestFixture& fixture) {
        BarrierOption option(fixture.standard_params, 120.0, BarrierType::UP_AND_OUT, true);
        
        const int num_runs = 10;
        
        {
            BENCHMARK("BarrierOption_10_runs_10k_sims");
            for (int i = 0; i < num_runs; ++i) {
                auto result = option.price(10000, 100);
                ASSERT_TRUE(result.is_valid);
            }
        }
    });
    
    // Benchmark rainbow option pricing
    suite->addTest("RainbowOptionBenchmark", []() {
        std::vector<double> spots = {100.0, 110.0};
        std::vector<double> vols = {0.20, 0.25};
        std::vector<std::vector<double>> corr = {{1.0, 0.5}, {0.5, 1.0}};
        
        MultiAssetParameters params(spots, vols, corr, 0.05, 0.0, 0.25);
        RainbowOption option(params, 100.0, true, true);
        
        const int num_runs = 5;
        
        {
            BENCHMARK("RainbowOption_5_runs_5k_sims");
            for (int i = 0; i < num_runs; ++i) {
                auto result = option.price(5000, 50);
                ASSERT_TRUE(result.is_valid);
            }
        }
    });
    
    TestRegistry::getInstance().registerSuite(std::move(suite));
}

// Thread safety tests
TEST_SUITE(ExoticOptionsThreadSafety) {
    auto suite = std::make_unique<TestSuite>("ExoticOptionsThreadSafety");
    
    // Test concurrent option pricing
    suite->addTestMethod<ExoticOptionsTestFixture>("ConcurrentPricing", [](ExoticOptionsTestFixture& fixture) {
        const int num_threads = 4;
        const int iterations_per_thread = 100;
        
        std::vector<std::thread> threads;
        std::vector<bool> results(num_threads, false);
        
        {
            BENCHMARK("ConcurrentPricing_4_threads_100_each");
            
            for (int t = 0; t < num_threads; ++t) {
                threads.emplace_back([&fixture, &results, t, iterations_per_thread]() {
                    bool all_valid = true;
                    
                    for (int i = 0; i < iterations_per_thread; ++i) {
                        AsianOption option(fixture.standard_params, true, true);
                        auto result = option.price(1000, 20);  // Smaller simulations for speed
                        
                        if (!result.is_valid) {
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
    
    // Test concurrent portfolio management
    suite->addTestMethod<ExoticOptionsTestFixture>("ConcurrentPortfolio", [](ExoticOptionsTestFixture& fixture) {
        ExoticPortfolio portfolio;
        
        const int num_threads = 3;
        std::vector<std::thread> threads;
        
        // Add options concurrently
        for (int t = 0; t < num_threads; ++t) {
            threads.emplace_back([&portfolio, &fixture, t]() {
                auto option = std::make_unique<AsianOption>(fixture.standard_params, true, true);
                portfolio.add_position(std::move(option), 1.0, "Thread_" + std::to_string(t));
            });
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
        
        ASSERT_EQ(portfolio.size(), num_threads);
        
        // Calculate portfolio value (should be thread-safe)
        double total_value = portfolio.calculate_total_value(1000);
        ASSERT_TRUE(std::isfinite(total_value));
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
        "exotic_options_test.log",
        5 * 1024 * 1024,  // 5MB max file size
        3       // keep 3 log files
    );
    
    // Initialize configuration
    Config::ConfigManager::getInstance().initialize("test_config.json");
    
    std::cout << "Exotic Options Pricing Models - Unit Tests" << std::endl;
    std::cout << "===========================================" << std::endl;
    
    // Print test discovery
    Testing::TestRegistry::getInstance().printDiscovery();
    
    // Run all tests
    auto stats = Testing::TestRegistry::getInstance().runAllSuites();
    
    // Return appropriate exit code
    return (stats.failed_tests == 0 && stats.error_tests == 0) ? 0 : 1;
}