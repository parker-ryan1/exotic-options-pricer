#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iomanip>
#include <string>
#include <map>
#include <memory>

// Mathematical constants and utilities
const double PI = 3.14159265358979323846;
const double SQRT_2PI = std::sqrt(2.0 * PI);
const double SQRT1_2 = 0.7071067811865476;  // 1/sqrt(2)

// Normal distribution functions
double norm_cdf(double x) {
    return 0.5 * std::erfc(-x * SQRT1_2);
}

double norm_pdf(double x) {
    return std::exp(-0.5 * x * x) / SQRT_2PI;
}

// Random number generator
class RandomGenerator {
private:
    std::mt19937 gen;
    std::normal_distribution<double> normal_dist;
    std::uniform_real_distribution<double> uniform_dist;

public:
    RandomGenerator() : gen(std::random_device{}()), normal_dist(0.0, 1.0), uniform_dist(0.0, 1.0) {}
    
    double normal() { return normal_dist(gen); }
    double uniform() { return uniform_dist(gen); }
    void seed(unsigned int s) { gen.seed(s); }
};

// Base Option class
class Option {
protected:
    double S0;      // Initial stock price
    double K;       // Strike price
    double T;       // Time to maturity
    double r;       // Risk-free rate
    double sigma;   // Volatility
    double q;       // Dividend yield

public:
    Option(double S0, double K, double T, double r, double sigma, double q = 0.0)
        : S0(S0), K(K), T(T), r(r), sigma(sigma), q(q) {}
    
    virtual ~Option() = default;
    virtual double price() = 0;
    virtual std::string getName() const = 0;
    
    // Getters
    double getS0() const { return S0; }
    double getK() const { return K; }
    double getT() const { return T; }
    double getR() const { return r; }
    double getSigma() const { return sigma; }
    double getQ() const { return q; }
};

// Asian Option (Arithmetic Average)
class AsianOption : public Option {
private:
    bool isCall;
    int numSteps;

public:
    AsianOption(double S0, double K, double T, double r, double sigma, bool isCall = true, 
                int numSteps = 252, double q = 0.0)
        : Option(S0, K, T, r, sigma, q), isCall(isCall), numSteps(numSteps) {}
    
    double price() override {
        const int numSims = 100000;
        double dt = T / numSteps;
        double drift = (r - q - 0.5 * sigma * sigma) * dt;
        double diffusion = sigma * std::sqrt(dt);
        
        RandomGenerator rng;
        double totalPayoff = 0.0;
        
        for (int sim = 0; sim < numSims; ++sim) {
            double S = S0;
            double sum = S0;
            
            for (int step = 1; step <= numSteps; ++step) {
                S *= std::exp(drift + diffusion * rng.normal());
                sum += S;
            }
            
            double avgPrice = sum / (numSteps + 1);
            double payoff = isCall ? std::max(avgPrice - K, 0.0) : std::max(K - avgPrice, 0.0);
            totalPayoff += payoff;
        }
        
        return std::exp(-r * T) * totalPayoff / numSims;
    }
    
    std::string getName() const override {
        return isCall ? "Asian Call Option" : "Asian Put Option";
    }
};

// Barrier Option
class BarrierOption : public Option {
public:
    enum BarrierType { UP_AND_OUT, UP_AND_IN, DOWN_AND_OUT, DOWN_AND_IN };
    
private:
    BarrierType barrierType;
    double barrier;
    bool isCall;
    int numSteps;

public:
    BarrierOption(double S0, double K, double T, double r, double sigma, double barrier,
                  BarrierType barrierType, bool isCall = true, int numSteps = 252, double q = 0.0)
        : Option(S0, K, T, r, sigma, q), barrierType(barrierType), barrier(barrier), 
          isCall(isCall), numSteps(numSteps) {}
    
    double price() override {
        const int numSims = 100000;
        double dt = T / numSteps;
        double drift = (r - q - 0.5 * sigma * sigma) * dt;
        double diffusion = sigma * std::sqrt(dt);
        
        RandomGenerator rng;
        double totalPayoff = 0.0;
        
        for (int sim = 0; sim < numSims; ++sim) {
            double S = S0;
            bool barrierHit = false;
            
            // Check initial condition
            if ((barrierType == UP_AND_OUT || barrierType == UP_AND_IN) && S >= barrier) {
                barrierHit = true;
            } else if ((barrierType == DOWN_AND_OUT || barrierType == DOWN_AND_IN) && S <= barrier) {
                barrierHit = true;
            }
            
            for (int step = 0; step < numSteps && !barrierHit; ++step) {
                S *= std::exp(drift + diffusion * rng.normal());
                
                if ((barrierType == UP_AND_OUT || barrierType == UP_AND_IN) && S >= barrier) {
                    barrierHit = true;
                } else if ((barrierType == DOWN_AND_OUT || barrierType == DOWN_AND_IN) && S <= barrier) {
                    barrierHit = true;
                }
            }
            
            double payoff = 0.0;
            bool shouldPayout = false;
            
            if (barrierType == UP_AND_OUT || barrierType == DOWN_AND_OUT) {
                shouldPayout = !barrierHit;  // Knock-out: pay if barrier NOT hit
            } else {
                shouldPayout = barrierHit;   // Knock-in: pay if barrier hit
            }
            
            if (shouldPayout) {
                payoff = isCall ? std::max(S - K, 0.0) : std::max(K - S, 0.0);
            }
            
            totalPayoff += payoff;
        }
        
        return std::exp(-r * T) * totalPayoff / numSims;
    }
    
    std::string getName() const override {
        std::string type;
        switch (barrierType) {
            case UP_AND_OUT: type = "Up-and-Out"; break;
            case UP_AND_IN: type = "Up-and-In"; break;
            case DOWN_AND_OUT: type = "Down-and-Out"; break;
            case DOWN_AND_IN: type = "Down-and-In"; break;
        }
        return type + (isCall ? " Call" : " Put") + " Barrier Option";
    }
};

// Lookback Option
class LookbackOption : public Option {
private:
    bool isCall;
    bool isFixed;  // Fixed strike vs floating strike
    int numSteps;

public:
    LookbackOption(double S0, double K, double T, double r, double sigma, bool isCall = true,
                   bool isFixed = true, int numSteps = 252, double q = 0.0)
        : Option(S0, K, T, r, sigma, q), isCall(isCall), isFixed(isFixed), numSteps(numSteps) {}
    
    double price() override {
        const int numSims = 100000;
        double dt = T / numSteps;
        double drift = (r - q - 0.5 * sigma * sigma) * dt;
        double diffusion = sigma * std::sqrt(dt);
        
        RandomGenerator rng;
        double totalPayoff = 0.0;
        
        for (int sim = 0; sim < numSims; ++sim) {
            double S = S0;
            double maxPrice = S0;
            double minPrice = S0;
            
            for (int step = 0; step < numSteps; ++step) {
                S *= std::exp(drift + diffusion * rng.normal());
                maxPrice = std::max(maxPrice, S);
                minPrice = std::min(minPrice, S);
            }
            
            double payoff = 0.0;
            if (isFixed) {
                // Fixed strike lookback
                if (isCall) {
                    payoff = std::max(maxPrice - K, 0.0);
                } else {
                    payoff = std::max(K - minPrice, 0.0);
                }
            } else {
                // Floating strike lookback
                if (isCall) {
                    payoff = S - minPrice;  // Always positive
                } else {
                    payoff = maxPrice - S;  // Always positive
                }
            }
            
            totalPayoff += payoff;
        }
        
        return std::exp(-r * T) * totalPayoff / numSims;
    }
    
    std::string getName() const override {
        std::string strikeType = isFixed ? "Fixed Strike" : "Floating Strike";
        std::string optionType = isCall ? "Call" : "Put";
        return strikeType + " Lookback " + optionType + " Option";
    }
};

// Digital/Binary Option
class DigitalOption : public Option {
private:
    bool isCall;
    double cashAmount;

public:
    DigitalOption(double S0, double K, double T, double r, double sigma, bool isCall = true,
                  double cashAmount = 1.0, double q = 0.0)
        : Option(S0, K, T, r, sigma, q), isCall(isCall), cashAmount(cashAmount) {}
    
    double price() override {
        // Analytical solution for digital options
        double d2 = (std::log(S0 / K) + (r - q - 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        
        double prob = isCall ? norm_cdf(d2) : norm_cdf(-d2);
        return std::exp(-r * T) * cashAmount * prob;
    }
    
    std::string getName() const override {
        return isCall ? "Digital Call Option" : "Digital Put Option";
    }
};

// Rainbow Option (Best/Worst of two assets)
class RainbowOption : public Option {
private:
    double S2;      // Second asset price
    double sigma2;  // Second asset volatility
    double rho;     // Correlation
    bool isCall;
    bool isBest;    // Best of two vs worst of two
    int numSteps;

public:
    RainbowOption(double S1, double S2, double K, double T, double r, 
                  double sigma1, double sigma2, double rho, bool isCall = true,
                  bool isBest = true, int numSteps = 252, double q = 0.0)
        : Option(S1, K, T, r, sigma1, q), S2(S2), sigma2(sigma2), rho(rho), 
          isCall(isCall), isBest(isBest), numSteps(numSteps) {}
    
    double price() override {
        const int numSims = 100000;
        double dt = T / numSteps;
        double drift1 = (r - q - 0.5 * sigma * sigma) * dt;
        double drift2 = (r - q - 0.5 * sigma2 * sigma2) * dt;
        double diffusion1 = sigma * std::sqrt(dt);
        double diffusion2 = sigma2 * std::sqrt(dt);
        
        RandomGenerator rng;
        double totalPayoff = 0.0;
        
        for (int sim = 0; sim < numSims; ++sim) {
            double S1 = S0;
            double S2_current = S2;
            
            for (int step = 0; step < numSteps; ++step) {
                double z1 = rng.normal();
                double z2 = rho * z1 + std::sqrt(1 - rho * rho) * rng.normal();
                
                S1 *= std::exp(drift1 + diffusion1 * z1);
                S2_current *= std::exp(drift2 + diffusion2 * z2);
            }
            
            double selectedPrice = isBest ? std::max(S1, S2_current) : std::min(S1, S2_current);
            double payoff = isCall ? std::max(selectedPrice - K, 0.0) : std::max(K - selectedPrice, 0.0);
            totalPayoff += payoff;
        }
        
        return std::exp(-r * T) * totalPayoff / numSims;
    }
    
    std::string getName() const override {
        std::string type = isBest ? "Best of Two" : "Worst of Two";
        std::string optionType = isCall ? "Call" : "Put";
        return type + " Rainbow " + optionType + " Option";
    }
};

// Portfolio and Risk Analysis
class PortfolioAnalyzer {
private:
    std::vector<std::unique_ptr<Option>> options;
    std::vector<double> positions;  // Number of contracts (can be negative for short)

public:
    void addOption(std::unique_ptr<Option> option, double position = 1.0) {
        options.push_back(std::move(option));
        positions.push_back(position);
    }
    
    double getTotalValue() {
        double total = 0.0;
        for (size_t i = 0; i < options.size(); ++i) {
            total += positions[i] * options[i]->price();
        }
        return total;
    }
    
    void printPortfolio() {
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "                    EXOTIC OPTIONS PORTFOLIO ANALYSIS" << std::endl;
        std::cout << std::string(80, '=') << std::endl;
        
        double totalValue = 0.0;
        
        for (size_t i = 0; i < options.size(); ++i) {
            double optionPrice = options[i]->price();
            double positionValue = positions[i] * optionPrice;
            totalValue += positionValue;
            
            std::cout << std::left << std::setw(35) << options[i]->getName() 
                      << " | Position: " << std::setw(8) << std::fixed << std::setprecision(2) << positions[i]
                      << " | Price: $" << std::setw(8) << optionPrice
                      << " | Value: $" << std::setw(10) << positionValue << std::endl;
        }
        
        std::cout << std::string(80, '-') << std::endl;
        std::cout << std::left << std::setw(35) << "TOTAL PORTFOLIO VALUE"
                  << " | " << std::setw(17) << " "
                  << " | " << std::setw(17) << " "
                  << " | Value: $" << std::setw(10) << totalValue << std::endl;
        std::cout << std::string(80, '=') << std::endl;
    }
};

// Interactive Menu System
class ExoticOptionsCalculator {
private:
    PortfolioAnalyzer portfolio;
    
    void printHeader() {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "           EXOTIC OPTIONS DERIVATIVE PRICER" << std::endl;
        std::cout << "                    C++ Edition" << std::endl;
        std::cout << std::string(60, '=') << std::endl;
    }
    
    void printMenu() {
        std::cout << "\nSelect Option Type:" << std::endl;
        std::cout << "1. Asian Option (Arithmetic Average)" << std::endl;
        std::cout << "2. Barrier Option (Knock-in/Knock-out)" << std::endl;
        std::cout << "3. Lookback Option (Fixed/Floating Strike)" << std::endl;
        std::cout << "4. Digital/Binary Option" << std::endl;
        std::cout << "5. Rainbow Option (Best/Worst of Two)" << std::endl;
        std::cout << "6. View Portfolio" << std::endl;
        std::cout << "7. Clear Portfolio" << std::endl;
        std::cout << "0. Exit" << std::endl;
        std::cout << "\nChoice: ";
    }
    
    double getInput(const std::string& prompt) {
        double value;
        std::cout << prompt;
        std::cin >> value;
        return value;
    }
    
    bool getBoolInput(const std::string& prompt) {
        char choice;
        std::cout << prompt << " (y/n): ";
        std::cin >> choice;
        return (choice == 'y' || choice == 'Y');
    }

public:
    void run() {
        printHeader();
        
        int choice;
        do {
            printMenu();
            std::cin >> choice;
            
            switch (choice) {
                case 1: createAsianOption(); break;
                case 2: createBarrierOption(); break;
                case 3: createLookbackOption(); break;
                case 4: createDigitalOption(); break;
                case 5: createRainbowOption(); break;
                case 6: portfolio.printPortfolio(); break;
                case 7: portfolio = PortfolioAnalyzer(); 
                       std::cout << "Portfolio cleared!" << std::endl; break;
                case 0: std::cout << "Goodbye!" << std::endl; break;
                default: std::cout << "Invalid choice!" << std::endl; break;
            }
        } while (choice != 0);
    }
    
private:
    void createAsianOption() {
        std::cout << "\n--- Asian Option Parameters ---" << std::endl;
        double S0 = getInput("Current Stock Price: $");
        double K = getInput("Strike Price: $");
        double T = getInput("Time to Maturity (years): ");
        double r = getInput("Risk-free Rate (decimal): ");
        double sigma = getInput("Volatility (decimal): ");
        bool isCall = getBoolInput("Call option?");
        double position = getInput("Position size (negative for short): ");
        
        auto option = std::make_unique<AsianOption>(S0, K, T, r, sigma, isCall);
        std::cout << "\nOption Price: $" << std::fixed << std::setprecision(4) << option->price() << std::endl;
        
        if (getBoolInput("Add to portfolio?")) {
            portfolio.addOption(std::move(option), position);
            std::cout << "Added to portfolio!" << std::endl;
        }
    }
    
    void createBarrierOption() {
        std::cout << "\n--- Barrier Option Parameters ---" << std::endl;
        double S0 = getInput("Current Stock Price: $");
        double K = getInput("Strike Price: $");
        double barrier = getInput("Barrier Level: $");
        double T = getInput("Time to Maturity (years): ");
        double r = getInput("Risk-free Rate (decimal): ");
        double sigma = getInput("Volatility (decimal): ");
        bool isCall = getBoolInput("Call option?");
        
        std::cout << "\nBarrier Types:" << std::endl;
        std::cout << "1. Up-and-Out  2. Up-and-In  3. Down-and-Out  4. Down-and-In" << std::endl;
        int barrierChoice;
        std::cout << "Choice: ";
        std::cin >> barrierChoice;
        
        BarrierOption::BarrierType barrierType = static_cast<BarrierOption::BarrierType>(barrierChoice - 1);
        double position = getInput("Position size (negative for short): ");
        
        auto option = std::make_unique<BarrierOption>(S0, K, T, r, sigma, barrier, barrierType, isCall);
        std::cout << "\nOption Price: $" << std::fixed << std::setprecision(4) << option->price() << std::endl;
        
        if (getBoolInput("Add to portfolio?")) {
            portfolio.addOption(std::move(option), position);
            std::cout << "Added to portfolio!" << std::endl;
        }
    }
    
    void createLookbackOption() {
        std::cout << "\n--- Lookback Option Parameters ---" << std::endl;
        double S0 = getInput("Current Stock Price: $");
        double K = getInput("Strike Price: $");
        double T = getInput("Time to Maturity (years): ");
        double r = getInput("Risk-free Rate (decimal): ");
        double sigma = getInput("Volatility (decimal): ");
        bool isCall = getBoolInput("Call option?");
        bool isFixed = getBoolInput("Fixed strike? (n for floating)");
        double position = getInput("Position size (negative for short): ");
        
        auto option = std::make_unique<LookbackOption>(S0, K, T, r, sigma, isCall, isFixed);
        std::cout << "\nOption Price: $" << std::fixed << std::setprecision(4) << option->price() << std::endl;
        
        if (getBoolInput("Add to portfolio?")) {
            portfolio.addOption(std::move(option), position);
            std::cout << "Added to portfolio!" << std::endl;
        }
    }
    
    void createDigitalOption() {
        std::cout << "\n--- Digital Option Parameters ---" << std::endl;
        double S0 = getInput("Current Stock Price: $");
        double K = getInput("Strike Price: $");
        double T = getInput("Time to Maturity (years): ");
        double r = getInput("Risk-free Rate (decimal): ");
        double sigma = getInput("Volatility (decimal): ");
        double cashAmount = getInput("Cash Amount (payout): $");
        bool isCall = getBoolInput("Call option?");
        double position = getInput("Position size (negative for short): ");
        
        auto option = std::make_unique<DigitalOption>(S0, K, T, r, sigma, isCall, cashAmount);
        std::cout << "\nOption Price: $" << std::fixed << std::setprecision(4) << option->price() << std::endl;
        
        if (getBoolInput("Add to portfolio?")) {
            portfolio.addOption(std::move(option), position);
            std::cout << "Added to portfolio!" << std::endl;
        }
    }
    
    void createRainbowOption() {
        std::cout << "\n--- Rainbow Option Parameters ---" << std::endl;
        double S1 = getInput("Asset 1 Price: $");
        double S2 = getInput("Asset 2 Price: $");
        double K = getInput("Strike Price: $");
        double T = getInput("Time to Maturity (years): ");
        double r = getInput("Risk-free Rate (decimal): ");
        double sigma1 = getInput("Asset 1 Volatility (decimal): ");
        double sigma2 = getInput("Asset 2 Volatility (decimal): ");
        double rho = getInput("Correlation (-1 to 1): ");
        bool isCall = getBoolInput("Call option?");
        bool isBest = getBoolInput("Best of two? (n for worst)");
        double position = getInput("Position size (negative for short): ");
        
        auto option = std::make_unique<RainbowOption>(S1, S2, K, T, r, sigma1, sigma2, rho, isCall, isBest);
        std::cout << "\nOption Price: $" << std::fixed << std::setprecision(4) << option->price() << std::endl;
        
        if (getBoolInput("Add to portfolio?")) {
            portfolio.addOption(std::move(option), position);
            std::cout << "Added to portfolio!" << std::endl;
        }
    }
};

int main() {
    ExoticOptionsCalculator calculator;
    calculator.run();
    return 0;
}