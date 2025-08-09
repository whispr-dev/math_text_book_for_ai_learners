#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>
#include <memory>
#include <limits>
#include <map>

namespace CalculusFunctions {

    const double PI = 3.141592653589793;
    const double E = 2.718281828459045;
    const double EPSILON = 1e-12;
    
    // ===== ENHANCED NUMERICAL INTEGRATION =====
    
    struct IntegrationResult {
        double value;
        double error_estimate;
        int function_evaluations;
    };
    
    // Gaussian Quadrature (5-point Gauss-Legendre)
    double gaussianQuadrature(const std::function<double(double)>& f, double a, double b) {
        // 5-point Gauss-Legendre weights and nodes
        static const std::vector<double> weights = {
            0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
            0.2369268850561891, 0.2369268850561891
        };
        static const std::vector<double> nodes = {
            0.0, -0.5384693101056831, 0.5384693101056831,
            -0.9061798459386640, 0.9061798459386640
        };
        
        double sum = 0.0;
        double half_length = (b - a) / 2.0;
        double mid_point = (a + b) / 2.0;
        
        for (size_t i = 0; i < weights.size(); ++i) {
            double x = mid_point + half_length * nodes[i];
            sum += weights[i] * f(x);
        }
        
        return half_length * sum;
    }
    
    // Adaptive Simpson's Rule with error control
    IntegrationResult adaptiveSimpson(const std::function<double(double)>& f, 
                                     double a, double b, double tolerance = 1e-10, int max_depth = 20) {
        int evaluations = 0;
        
        std::function<double(double, double, double, double, double, int)> adaptive;
        adaptive = [&](double a, double b, double fa, double fb, double fc, int depth) -> double {
            double c = (a + b) / 2.0;
            double h = b - a;
            double d = (a + c) / 2.0;
            double e = (c + b) / 2.0;
            double fd = f(d);
            double fe = f(e);
            evaluations += 2;
            
            double S1 = h * (fa + 4*fc + fb) / 6.0;  // Simpson on [a,b]
            double S2 = h * (fa + 4*fd + fc) / 12.0 + h * (fc + 4*fe + fb) / 12.0;  // Simpson on [a,c] + [c,b]
            
            if (depth >= max_depth || std::abs(S2 - S1) <= 15 * tolerance) {
                return S2 + (S2 - S1) / 15.0;  // Richardson extrapolation
            } else {
                return adaptive(a, c, fa, fd, fc, depth + 1) + 
                       adaptive(c, b, fc, fe, fb, depth + 1);
            }
        };
        
        double fa = f(a);
        double fb = f(b);
        double fc = f((a + b) / 2.0);
        evaluations += 3;
        
        double result = adaptive(a, b, fa, fb, fc, 0);
        double error = tolerance * std::abs(result);
        
        return {result, error, evaluations};
    }
    
    // Romberg Integration
    IntegrationResult rombergIntegration(const std::function<double(double)>& f, 
                                       double a, double b, int max_iterations = 10) {
        std::vector<std::vector<double>> R(max_iterations, std::vector<double>(max_iterations));
        double h = b - a;
        int evaluations = 2;
        
        // Initial trapezoidal rule
        R[0][0] = h * (f(a) + f(b)) / 2.0;
        
        for (int i = 1; i < max_iterations; ++i) {
            // Improve trapezoidal rule
            h /= 2.0;
            double sum = 0.0;
            for (int k = 1; k <= std::pow(2, i-1); ++k) {
                sum += f(a + (2*k - 1) * h);
                evaluations++;
            }
            R[i][0] = R[i-1][0] / 2.0 + h * sum;
            
            // Richardson extrapolation
            for (int j = 1; j <= i; ++j) {
                double coeff = std::pow(4, j);
                R[i][j] = (coeff * R[i][j-1] - R[i-1][j-1]) / (coeff - 1);
            }
            
            // Check convergence
            if (i > 1 && std::abs(R[i][i] - R[i-1][i-1]) < EPSILON) {
                return {R[i][i], std::abs(R[i][i] - R[i-1][i-1]), evaluations};
            }
        }
        
        return {R[max_iterations-1][max_iterations-1], EPSILON, evaluations};
    }
    
    // ===== SYMBOLIC INTEGRATION =====
    
    class SymbolicExpression {
    public:
        virtual ~SymbolicExpression() = default;
        virtual std::unique_ptr<SymbolicExpression> clone() const = 0;
        virtual std::unique_ptr<SymbolicExpression> differentiate(const std::string& var) const = 0;
        virtual std::unique_ptr<SymbolicExpression> integrate(const std::string& var) const = 0;
        virtual std::string toString() const = 0;
        virtual double evaluate(const std::map<std::string, double>& variables) const = 0;
    };
    
    class Constant : public SymbolicExpression {
        double value;
    public:
        Constant(double v) : value(v) {}
        
        std::unique_ptr<SymbolicExpression> clone() const override {
            return std::make_unique<Constant>(value);
        }
        
        std::unique_ptr<SymbolicExpression> differentiate(const std::string&) const override {
            return std::make_unique<Constant>(0);
        }
        
        std::unique_ptr<SymbolicExpression> integrate(const std::string& var) const override {
            return std::make_unique<Product>(clone(), std::make_unique<Variable>(var));
        }
        
        std::string toString() const override {
            return std::to_string(value);
        }
        
        double evaluate(const std::map<std::string, double>&) const override {
            return value;
        }
        
        double getValue() const { return value; }
    };
    
    class Variable : public SymbolicExpression {
        std::string name;
    public:
        Variable(const std::string& n) : name(n) {}
        
        std::unique_ptr<SymbolicExpression> clone() const override {
            return std::make_unique<Variable>(name);
        }
        
        std::unique_ptr<SymbolicExpression> differentiate(const std::string& var) const override {
            return std::make_unique<Constant>(name == var ? 1.0 : 0.0);
        }
        
        std::unique_ptr<SymbolicExpression> integrate(const std::string& var) const override {
            if (name == var) {
                // ∫x dx = x²/2
                return std::make_unique<Product>(
                    std::make_unique<Constant>(0.5),
                    std::make_unique<Power>(clone(), std::make_unique<Constant>(2))
                );
            } else {
                // ∫a dx = ax (where a is constant w.r.t. x)
                return std::make_unique<Product>(clone(), std::make_unique<Variable>(var));
            }
        }
        
        std::string toString() const override {
            return name;
        }
        
        double evaluate(const std::map<std::string, double>& variables) const override {
            auto it = variables.find(name);
            return it != variables.end() ? it->second : 0.0;
        }
        
        const std::string& getName() const { return name; }
    };
    
    class Sum : public SymbolicExpression {
        std::unique_ptr<SymbolicExpression> left, right;
    public:
        Sum(std::unique_ptr<SymbolicExpression> l, std::unique_ptr<SymbolicExpression> r)
            : left(std::move(l)), right(std::move(r)) {}
        
        std::unique_ptr<SymbolicExpression> clone() const override {
            return std::make_unique<Sum>(left->clone(), right->clone());
        }
        
        std::unique_ptr<SymbolicExpression> differentiate(const std::string& var) const override {
            return std::make_unique<Sum>(left->differentiate(var), right->differentiate(var));
        }
        
        std::unique_ptr<SymbolicExpression> integrate(const std::string& var) const override {
            return std::make_unique<Sum>(left->integrate(var), right->integrate(var));
        }
        
        std::string toString() const override {
            return "(" + left->toString() + " + " + right->toString() + ")";
        }
        
        double evaluate(const std::map<std::string, double>& variables) const override {
            return left->evaluate(variables) + right->evaluate(variables);
        }
    };
    
    class Product : public SymbolicExpression {
        std::unique_ptr<SymbolicExpression> left, right;
    public:
        Product(std::unique_ptr<SymbolicExpression> l, std::unique_ptr<SymbolicExpression> r)
            : left(std::move(l)), right(std::move(r)) {}
        
        std::unique_ptr<SymbolicExpression> clone() const override {
            return std::make_unique<Product>(left->clone(), right->clone());
        }
        
        std::unique_ptr<SymbolicExpression> differentiate(const std::string& var) const override {
            // Product rule: (fg)' = f'g + fg'
            return std::make_unique<Sum>(
                std::make_unique<Product>(left->differentiate(var), right->clone()),
                std::make_unique<Product>(left->clone(), right->differentiate(var))
            );
        }
        
        std::unique_ptr<SymbolicExpression> integrate(const std::string& var) const override {
            // Try simple cases first
            auto leftConst = dynamic_cast<const Constant*>(left.get());
            auto rightConst = dynamic_cast<const Constant*>(right.get());
            
            if (leftConst) {
                // ∫c*f(x) dx = c*∫f(x) dx
                return std::make_unique<Product>(left->clone(), right->integrate(var));
            } else if (rightConst) {
                // ∫f(x)*c dx = c*∫f(x) dx  
                return std::make_unique<Product>(right->clone(), left->integrate(var));
            }
            
            // For more complex products, we'd need integration by parts
            throw std::runtime_error("Complex product integration not implemented");
        }
        
        std::string toString() const override {
            return "(" + left->toString() + " * " + right->toString() + ")";
        }
        
        double evaluate(const std::map<std::string, double>& variables) const override {
            return left->evaluate(variables) * right->evaluate(variables);
        }
    };
    
    class Power : public SymbolicExpression {
        std::unique_ptr<SymbolicExpression> base;
        std::unique_ptr<SymbolicExpression> exponent;
    public:
        Power(std::unique_ptr<SymbolicExpression> b, std::unique_ptr<SymbolicExpression> e)
            : base(std::move(b)), exponent(std::move(e)) {}
        
        std::unique_ptr<SymbolicExpression> clone() const override {
            return std::make_unique<Power>(base->clone(), exponent->clone());
        }
        
        std::unique_ptr<SymbolicExpression> differentiate(const std::string& var) const override {
            auto expConst = dynamic_cast<const Constant*>(exponent.get());
            if (expConst) {
                // Power rule: (x^n)' = n*x^(n-1)*x'
                double n = expConst->getValue();
                if (n == 0) return std::make_unique<Constant>(0);
                
                return std::make_unique<Product>(
                    std::make_unique<Product>(
                        std::make_unique<Constant>(n),
                        std::make_unique<Power>(base->clone(), std::make_unique<Constant>(n-1))
                    ),
                    base->differentiate(var)
                );
            }
            throw std::runtime_error("General power differentiation not implemented");
        }
        
        std::unique_ptr<SymbolicExpression> integrate(const std::string& var) const override {
            auto expConst = dynamic_cast<const Constant*>(exponent.get());
            auto baseVar = dynamic_cast<const Variable*>(base.get());
            
            if (expConst && baseVar && baseVar->getName() == var) {
                double n = expConst->getValue();
                if (std::abs(n + 1) < EPSILON) {
                    throw std::runtime_error("∫x^(-1) dx = ln(x) not implemented");
                }
                // ∫x^n dx = x^(n+1)/(n+1)
                return std::make_unique<Product>(
                    std::make_unique<Constant>(1.0 / (n + 1)),
                    std::make_unique<Power>(base->clone(), std::make_unique<Constant>(n + 1))
                );
            }
            throw std::runtime_error("General power integration not implemented");
        }
        
        std::string toString() const override {
            return "(" + base->toString() + "^" + exponent->toString() + ")";
        }
        
        double evaluate(const std::map<std::string, double>& variables) const override {
            return std::pow(base->evaluate(variables), exponent->evaluate(variables));
        }
    };
    
    // ===== HIGHER-ORDER DERIVATIVES =====
    
    std::vector<double> higherOrderDerivative(const std::vector<double>& x, 
                                            const std::vector<double>& y, 
                                            int order, int point_index) {
        if (order <= 0 || point_index < 0 || point_index >= static_cast<int>(x.size())) {
            throw std::invalid_argument("Invalid parameters for higher-order derivative");
        }
        
        std::vector<double> current_y = y;
        
        for (int ord = 0; ord < order; ++ord) {
            std::vector<double> next_y(current_y.size() - 1);
            
            for (size_t i = 0; i < next_y.size(); ++i) {
                double h = (i + 1 < x.size()) ? x[i + 1] - x[i] : x[i] - x[i - 1];
                next_y[i] = (current_y[i + 1] - current_y[i]) / h;
            }
            
            current_y = next_y;
        }
        
        return current_y;
    }
    
    // ===== TAYLOR SERIES =====
    
    struct TaylorSeries {
        std::vector<double> coefficients;
        double center;
        int terms;
        
        double evaluate(double x) const {
            double result = 0.0;
            double dx = x - center;
            double power = 1.0;
            
            for (int i = 0; i < terms && i < static_cast<int>(coefficients.size()); ++i) {
                result += coefficients[i] * power;
                power *= dx;
            }
            return result;
        }
    };
    
    TaylorSeries taylorExpansion(const std::function<double(double)>& f, 
                               double center, int terms, double h = 1e-6) {
        TaylorSeries series;
        series.center = center;
        series.terms = terms;
        series.coefficients.resize(terms);
        
        // Calculate derivatives using finite differences
        std::vector<double> factorial(terms, 1.0);
        for (int i = 1; i < terms; ++i) {
            factorial[i] = factorial[i-1] * i;
        }
        
        for (int n = 0; n < terms; ++n) {
            if (n == 0) {
                series.coefficients[0] = f(center);
            } else {
                // Use central difference for higher derivatives
                double derivative = 0.0;
                double sign = 1.0;
                
                for (int k = 0; k <= n; ++k) {
                    double binomial = 1.0;
                    for (int j = 1; j <= k; ++j) {
                        binomial = binomial * (n - j + 1) / j;
                    }
                    derivative += sign * binomial * f(center + (n/2.0 - k) * h);
                    sign *= -1.0;
                }
                
                derivative /= std::pow(h, n);
                series.coefficients[n] = derivative / factorial[n];
            }
        }
        
        return series;
    }
    
    // ===== ROOT FINDING & OPTIMIZATION =====
    
    struct OptimizationResult {
        double value;
        double function_value;
        int iterations;
        bool converged;
    };
    
    // Newton-Raphson root finding
    OptimizationResult newtonRaphson(const std::function<double(double)>& f,
                                   const std::function<double(double)>& df,
                                   double initial_guess, double tolerance = 1e-10,
                                   int max_iterations = 100) {
        double x = initial_guess;
        
        for (int i = 0; i < max_iterations; ++i) {
            double fx = f(x);
            double dfx = df(x);
            
            if (std::abs(dfx) < tolerance) {
                return {x, fx, i, false}; // Derivative too small
            }
            
            double x_new = x - fx / dfx;
            
            if (std::abs(x_new - x) < tolerance) {
                return {x_new, f(x_new), i + 1, true};
            }
            
            x = x_new;
        }
        
        return {x, f(x), max_iterations, false};
    }
    
    // Bisection method (more robust)
    OptimizationResult bisection(const std::function<double(double)>& f,
                               double a, double b, double tolerance = 1e-10,
                               int max_iterations = 100) {
        if (f(a) * f(b) >= 0) {
            throw std::invalid_argument("Function must have opposite signs at endpoints");
        }
        
        for (int i = 0; i < max_iterations; ++i) {
            double c = (a + b) / 2.0;
            double fc = f(c);
            
            if (std::abs(fc) < tolerance || std::abs(b - a) < tolerance) {
                return {c, fc, i + 1, true};
            }
            
            if (f(a) * fc < 0) {
                b = c;
            } else {
                a = c;
            }
        }
        
        double c = (a + b) / 2.0;
        return {c, f(c), max_iterations, false};
    }
    
    // Golden section search for optimization
    OptimizationResult goldenSectionMinimize(const std::function<double(double)>& f,
                                           double a, double b, double tolerance = 1e-10,
                                           int max_iterations = 100) {
        const double phi = (1.0 + std::sqrt(5.0)) / 2.0; // Golden ratio
        const double resphi = 2.0 - phi;
        
        double x1 = a + resphi * (b - a);
        double x2 = b - resphi * (b - a);
        double f1 = f(x1);
        double f2 = f(x2);
        
        for (int i = 0; i < max_iterations; ++i) {
            if (std::abs(b - a) < tolerance) {
                double x_min = (a + b) / 2.0;
                return {x_min, f(x_min), i + 1, true};
            }
            
            if (f1 < f2) {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = a + resphi * (b - a);
                f1 = f(x1);
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = b - resphi * (b - a);
                f2 = f(x2);
            }
        }
        
        double x_min = (a + b) / 2.0;
        return {x_min, f(x_min), max_iterations, false};
    }
    
} // namespace CalculusFunctions

// ===== DEMONSTRATION =====

int main() {
    using namespace CalculusFunctions;
    
    std::cout << "=== ENHANCED CALCULUS FUNCTIONS DEMO ===" << std::endl << std::endl;
    
    // Test function: f(x) = x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)
    auto test_func = [](double x) { return x*x*x - 6*x*x + 11*x - 6; };
    auto test_func_derivative = [](double x) { return 3*x*x - 12*x + 11; };
    
    // === ADVANCED INTEGRATION ===
    std::cout << "1. ADVANCED NUMERICAL INTEGRATION" << std::endl;
    
    // Compare integration methods
    double a = 0.0, b = 4.0;
    
    // Gaussian quadrature
    double gauss_result = gaussianQuadrature(test_func, a, b);
    std::cout << "Gaussian Quadrature: " << gauss_result << std::endl;
    
    // Adaptive Simpson
    auto adaptive_result = adaptiveSimpson(test_func, a, b, 1e-10);
    std::cout << "Adaptive Simpson: " << adaptive_result.value 
              << " (error: " << adaptive_result.error_estimate 
              << ", evals: " << adaptive_result.function_evaluations << ")" << std::endl;
    
    // Romberg integration
    auto romberg_result = rombergIntegration(test_func, a, b, 8);
    std::cout << "Romberg: " << romberg_result.value 
              << " (error: " << romberg_result.error_estimate 
              << ", evals: " << romberg_result.function_evaluations << ")" << std::endl;
    std::cout << std::endl;
    
    // === SYMBOLIC INTEGRATION ===
    std::cout << "2. SYMBOLIC INTEGRATION" << std::endl;
    
    // Create symbolic expressions: x^3 + 2x^2 + x + 1
    auto x = std::make_unique<Variable>("x");
    auto expr = std::make_unique<Sum>(
        std::make_unique<Sum>(
            std::make_unique<Sum>(
                std::make_unique<Power>(x->clone(), std::make_unique<Constant>(3)),
                std::make_unique<Product>(std::make_unique<Constant>(2), 
                                        std::make_unique<Power>(x->clone(), std::make_unique<Constant>(2)))
            ),
            x->clone()
        ),
        std::make_unique<Constant>(1)
    );
    
    std::cout << "Expression: " << expr->toString() << std::endl;
    
    auto derivative = expr->differentiate("x");
    std::cout << "Derivative: " << derivative->toString() << std::endl;
    
    auto integral = expr->integrate("x");  
    std::cout << "Integral: " << integral->toString() << std::endl;
    
    // Evaluate at x = 2
    std::map<std::string, double> vars = {{"x", 2.0}};
    std::cout << "f(2) = " << expr->evaluate(vars) << std::endl;
    std::cout << "f'(2) = " << derivative->evaluate(vars) << std::endl;
    std::cout << std::endl;
    
    // === TAYLOR SERIES ===
    std::cout << "3. TAYLOR SERIES EXPANSION" << std::endl;
    
    auto sin_func = [](double x) { return std::sin(x); };
    auto taylor_sin = taylorExpansion(sin_func, 0.0, 6); // sin(x) around x=0
    
    std::cout << "Taylor series for sin(x) around x=0:" << std::endl;
    for (int i = 0; i < taylor_sin.terms; ++i) {
        std::cout << "  a" << i << " = " << taylor_sin.coefficients[i] << std::endl;
    }
    
    double test_x = PI / 6.0;
    std::cout << "sin(" << test_x << ") = " << std::sin(test_x) << " (exact)" << std::endl;
    std::cout << "Taylor approximation = " << taylor_sin.evaluate(test_x) << std::endl;
    std::cout << std::endl;
    
    // === ROOT FINDING ===
    std::cout << "4. ROOT FINDING & OPTIMIZATION" << std::endl;
    
    // Find roots of test function
    auto newton_result = newtonRaphson(test_func, test_func_derivative, 0.5);
    std::cout << "Newton-Raphson root near 0.5: " << newton_result.value 
              << " (f=" << newton_result.function_value 
              << ", converged=" << newton_result.converged << ")" << std::endl;
    
    auto bisection_result = bisection(test_func, 0.5, 1.5);
    std::cout << "Bisection root in [0.5, 1.5]: " << bisection_result.value 
              << " (f=" << bisection_result.function_value 
              << ", converged=" << bisection_result.converged << ")" << std::endl;
    
    // Find minimum
    auto min_result = goldenSectionMinimize(test_func, 0.0, 4.0);
    std::cout << "Minimum at x = " << min_result.value 
              << " (f=" << min_result.function_value 
              << ", converged=" << min_result.converged << ")" << std::endl;
    std::cout << std::endl;
    
    // === PERFORMANCE COMPARISON ===
    std::cout << "5. PERFORMANCE COMPARISON" << std::endl;
    std::cout << "Function evaluations for ∫₀⁴ f(x)dx:" << std::endl;
    std::cout << "  Adaptive Simpson: " << adaptive_result.function_evaluations << std::endl;
    std::cout << "  Romberg: " << romberg_result.function_evaluations << std::endl;
    std::cout << "  Gaussian (5-point): 5 (fixed)" << std::endl;
    
    std::cout << std::endl << "=== CALCULUS MODULE ENHANCED ===" << std::endl;
    
    return 0;
}