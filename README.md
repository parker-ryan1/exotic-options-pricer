# Exotic Options Derivative Pricer (C++)

A comprehensive C++ application for pricing exotic derivative instruments using Monte Carlo simulation and analytical methods.

## Features

### Supported Exotic Options:
1. **Asian Options** - Arithmetic average price options
2. **Barrier Options** - Knock-in/knock-out with up/down barriers
3. **Lookback Options** - Fixed and floating strike variants
4. **Digital/Binary Options** - Cash-or-nothing payoffs
5. **Rainbow Options** - Best/worst of two correlated assets

### Advanced Capabilities:
- **Monte Carlo Simulation** - High-performance pricing engine
- **Portfolio Management** - Build and analyze complex portfolios
- **Risk Analytics** - Position sizing and portfolio valuation
- **Interactive CLI** - User-friendly command-line interface
- **Optimized Performance** - Compiled C++ with O3 optimization

## Mathematical Models

### Asian Options
- **Payoff**: max(S_avg - K, 0) for calls, max(K - S_avg, 0) for puts
- **Method**: Monte Carlo with arithmetic averaging
- **Applications**: Commodity derivatives, currency options

### Barrier Options
- **Types**: Up-and-Out, Up-and-In, Down-and-Out, Down-and-In
- **Payoff**: Standard European payoff if barrier condition met
- **Method**: Path-dependent Monte Carlo simulation

### Lookback Options
- **Fixed Strike**: max(S_max - K, 0) or max(K - S_min, 0)
- **Floating Strike**: S_T - S_min or S_max - S_T
- **Method**: Monte Carlo with extrema tracking

### Digital Options
- **Payoff**: Fixed cash amount if S_T > K (call) or S_T < K (put)
- **Method**: Analytical Black-Scholes formula
- **Use Cases**: Binary bets, structured products

### Rainbow Options
- **Payoff**: Based on best/worst performing asset
- **Method**: Correlated Monte Carlo simulation
- **Applications**: Multi-asset derivatives, basket options

## Build Instructions

### Prerequisites:
- C++17 compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- Make utility (optional but recommended)

### Compilation:

#### Using Make (Recommended):
```bash
# Standard build
make

# Optimized release build
make release

# Debug build
make debug

# Clean build files
make clean

# Build and run
make run
```

#### Manual Compilation:
```bash
# Windows (MinGW/MSYS2)
g++ -std=c++17 -O3 -Wall -Wextra -o exotic_options.exe main.cpp

# Linux/macOS
g++ -std=c++17 -O3 -Wall -Wextra -march=native -o exotic_options main.cpp
```

## Usage

### Running the Application:
```bash
# Windows
./exotic_options.exe

# Linux/macOS
./exotic_options
```

### Interactive Menu:
1. Select option type from the menu
2. Enter market parameters (spot price, strike, volatility, etc.)
3. Choose position size (positive for long, negative for short)
4. Add to portfolio or calculate individual prices
5. View portfolio summary and total value

### Example Session:
```
Select Option Type:
1. Asian Option (Arithmetic Average)
2. Barrier Option (Knock-in/Knock-out)
3. Lookback Option (Fixed/Floating Strike)
4. Digital/Binary Option
5. Rainbow Option (Best/Worst of Two)
6. View Portfolio
7. Clear Portfolio
0. Exit

Choice: 1

--- Asian Option Parameters ---
Current Stock Price: $100
Strike Price: $105
Time to Maturity (years): 0.25
Risk-free Rate (decimal): 0.05
Volatility (decimal): 0.20
Call option? (y/n): y
Position size (negative for short): 10

Option Price: $2.3456
Add to portfolio? (y/n): y
Added to portfolio!
```

## Performance Optimizations

- **Vectorized Operations**: Efficient random number generation
- **Memory Management**: Smart pointers and RAII principles
- **Compiler Optimizations**: -O3 flag with native architecture targeting
- **Monte Carlo Efficiency**: Optimized simulation loops

## Risk Management Features

- **Portfolio Aggregation**: Sum individual option values
- **Position Tracking**: Long/short position management
- **Real-time Pricing**: Instant recalculation capabilities
- **Scenario Analysis**: Multiple option types in single portfolio

## Technical Implementation

### Random Number Generation:
- Mersenne Twister MT19937 generator
- Box-Muller transformation for normal variates
- Correlation handling for multi-asset options

### Numerical Methods:
- Monte Carlo simulation (100,000 paths default)
- Analytical solutions where available
- Error function approximations for normal CDF

### Code Architecture:
- Object-oriented design with inheritance
- Template-based generic programming
- RAII memory management
- Exception-safe implementations

## Extending the System

### Adding New Option Types:
1. Inherit from the `Option` base class
2. Implement `price()` and `getName()` methods
3. Add to the menu system in `ExoticOptionsCalculator`
4. Update the Makefile if needed

### Example New Option:
```cpp
class CustomOption : public Option {
public:
    CustomOption(/* parameters */) : Option(/* base params */) {}
    
    double price() override {
        // Your pricing logic here
        return calculatedPrice;
    }
    
    std::string getName() const override {
        return "Custom Option";
    }
};
```

## Validation and Testing

The pricing models have been validated against:
- Known analytical solutions (where available)
- Industry-standard pricing libraries
- Academic literature benchmarks
- Monte Carlo convergence analysis

## License

This project is provided for educational and research purposes. Use at your own risk for production trading systems.