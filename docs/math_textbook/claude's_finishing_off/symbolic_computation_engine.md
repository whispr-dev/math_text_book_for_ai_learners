#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <functional>
#include <map>
#include <set>

namespace SymbolicComputation {

    const double PI = 3.141592653589793;
    const double E = 2.718281828459045;
    const double EPSILON = 1e-12;
    
    // Forward declarations
    class Expression;
    using ExprPtr = std::shared_ptr<Expression>;
    using VarMap = std::unordered_map<std::string, double>;
    
    // ===== BASE EXPRESSION CLASS =====
    
    class Expression {
    public:
        virtual ~Expression() = default;
        virtual ExprPtr clone() const = 0;
        virtual ExprPtr differentiate(const std::string& var) const = 0;
        virtual ExprPtr integrate(const std::string& var) const = 0;
        virtual ExprPtr simplify() const = 0;
        virtual std::string toString() const = 0;
        virtual double evaluate(const VarMap& vars = {}) const = 0;
        virtual std::set<std::string> getVariables() const = 0;
        virtual bool equals(const ExprPtr& other) const = 0;
        virtual ExprPtr substitute(const std::string& var, const ExprPtr& replacement) const = 0;
        
        // Utility functions
        virtual bool isZero() const { return false; }
        virtual bool isOne() const { return false; }
        virtual bool isConstant() const { return getVariables().empty(); }
        virtual int getPrecedence() const { return 0; }
        
        // Operator overloads for expression building
        ExprPtr simplify() const override {
            auto arg_simp = argument->simplify();
            
            // exp(0) = 1
            if (arg_simp->isZero()) {
                return std::make_shared<Constant>(1);
            }
            
            // exp(ln(x)) = x
            if (auto ln_arg = std::dynamic_pointer_cast<NaturalLog>(arg_simp)) {
                return ln_arg->getArgument();
            }
            
            return std::make_shared<Exponential>(arg_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return std::exp(argument->evaluate(vars));
        }
        
    protected:
        ExprPtr createFunction(const ExprPtr& arg) const override {
            return std::make_shared<Exponential>(arg);
        }
    };
    
    class NaturalLog : public Function {
    public:
        NaturalLog(ExprPtr arg) : Function(arg, "ln") {}
        
        ExprPtr clone() const override {
            return std::make_shared<NaturalLog>(argument->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            // (ln(f))' = f'/f
            return std::make_shared<Division>(argument->differentiate(var), argument);
        }
        
        ExprPtr integrate(const std::string& var) const override {
            // ∫ln(x)dx = x*ln(x) - x
            if (auto var_arg = std::dynamic_pointer_cast<Variable>(argument)) {
                if (var_arg->getName() == var) {
                    auto x_ln_x = std::make_shared<Multiplication>(argument, clone());
                    return std::make_shared<Subtraction>(x_ln_x, argument);
                }
            }
            throw std::runtime_error("Integration of ln with complex argument not implemented");
        }
        
        ExprPtr simplify() const override {
            auto arg_simp = argument->simplify();
            
            // ln(1) = 0
            if (arg_simp->isOne()) {
                return std::make_shared<Constant>(0);
            }
            
            // ln(e) = 1
            auto e_const = std::dynamic_pointer_cast<Constant>(arg_simp);
            if (e_const && std::abs(e_const->getValue() - E) < EPSILON) {
                return std::make_shared<Constant>(1);
            }
            
            // ln(exp(x)) = x
            if (auto exp_arg = std::dynamic_pointer_cast<Exponential>(arg_simp)) {
                return exp_arg->getArgument();
            }
            
            return std::make_shared<NaturalLog>(arg_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            double val = argument->evaluate(vars);
            if (val <= 0) {
                throw std::runtime_error("Natural log of non-positive number");
            }
            return std::log(val);
        }
        
    protected:
        ExprPtr createFunction(const ExprPtr& arg) const override {
            return std::make_shared<NaturalLog>(arg);
        }
    };
    
    // ===== HYPERBOLIC FUNCTIONS =====
    
    class Sinh : public Function {
    public:
        Sinh(ExprPtr arg) : Function(arg, "sinh") {}
        
        ExprPtr clone() const override {
            return std::make_shared<Sinh>(argument->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            // (sinh(f))' = cosh(f) * f'
            return std::make_shared<Multiplication>(
                std::make_shared<Cosh>(argument),
                argument->differentiate(var)
            );
        }
        
        ExprPtr integrate(const std::string& var) const override {
            // ∫sinh(x)dx = cosh(x)
            if (auto var_arg = std::dynamic_pointer_cast<Variable>(argument)) {
                if (var_arg->getName() == var) {
                    return std::make_shared<Cosh>(argument);
                }
            }
            throw std::runtime_error("Integration of sinh with complex argument not implemented");
        }
        
        ExprPtr simplify() const override {
            auto arg_simp = argument->simplify();
            
            // sinh(0) = 0
            if (arg_simp->isZero()) {
                return std::make_shared<Constant>(0);
            }
            
            return std::make_shared<Sinh>(arg_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return std::sinh(argument->evaluate(vars));
        }
        
    protected:
        ExprPtr createFunction(const ExprPtr& arg) const override {
            return std::make_shared<Sinh>(arg);
        }
    };
    
    class Cosh : public Function {
    public:
        Cosh(ExprPtr arg) : Function(arg, "cosh") {}
        
        ExprPtr clone() const override {
            return std::make_shared<Cosh>(argument->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            // (cosh(f))' = sinh(f) * f'
            return std::make_shared<Multiplication>(
                std::make_shared<Sinh>(argument),
                argument->differentiate(var)
            );
        }
        
        ExprPtr integrate(const std::string& var) const override {
            // ∫cosh(x)dx = sinh(x)
            if (auto var_arg = std::dynamic_pointer_cast<Variable>(argument)) {
                if (var_arg->getName() == var) {
                    return std::make_shared<Sinh>(argument);
                }
            }
            throw std::runtime_error("Integration of cosh with complex argument not implemented");
        }
        
        ExprPtr simplify() const override {
            auto arg_simp = argument->simplify();
            
            // cosh(0) = 1
            if (arg_simp->isZero()) {
                return std::make_shared<Constant>(1);
            }
            
            return std::make_shared<Cosh>(arg_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return std::cosh(argument->evaluate(vars));
        }
        
    protected:
        ExprPtr createFunction(const ExprPtr& arg) const override {
            return std::make_shared<Cosh>(arg);
        }
    };
    
    // ===== OPERATOR IMPLEMENTATIONS =====
    
    ExprPtr Expression::operator+(const ExprPtr& other) const {
        return std::make_shared<Addition>(shared_from_this(), other);
    }
    
    ExprPtr Expression::operator-(const ExprPtr& other) const {
        return std::make_shared<Subtraction>(shared_from_this(), other);
    }
    
    ExprPtr Expression::operator*(const ExprPtr& other) const {
        return std::make_shared<Multiplication>(shared_from_this(), other);
    }
    
    ExprPtr Expression::operator/(const ExprPtr& other) const {
        return std::make_shared<Division>(shared_from_this(), other);
    }
    
    ExprPtr Expression::operator^(const ExprPtr& other) const {
        return std::make_shared<Power>(shared_from_this(), other);
    }
    
    // ===== INTEGRATION IMPLEMENTATIONS =====
    
    ExprPtr Constant::integrate(const std::string& var) const {
        // ∫c dx = c*x
        return std::make_shared<Multiplication>(clone(), std::make_shared<Variable>(var));
    }
    
    ExprPtr Variable::integrate(const std::string& var) const {
        if (name == var) {
            // ∫x dx = x²/2
            auto x_squared = std::make_shared<Power>(clone(), std::make_shared<Constant>(2));
            return std::make_shared<Division>(x_squared, std::make_shared<Constant>(2));
        } else {
            // ∫a dx = a*x (where a is constant w.r.t. x)
            return std::make_shared<Multiplication>(clone(), std::make_shared<Variable>(var));
        }
    }
    
    ExprPtr Multiplication::integrate(const std::string& var) const override {
        // Try simple cases first
        auto l_const = std::dynamic_pointer_cast<Constant>(left);
        auto r_const = std::dynamic_pointer_cast<Constant>(right);
        
        if (l_const && !left->getVariables().count(var)) {
            // ∫c*f(x) dx = c*∫f(x) dx
            return std::make_shared<Multiplication>(left, right->integrate(var));
        } else if (r_const && !right->getVariables().count(var)) {
            // ∫f(x)*c dx = c*∫f(x) dx
            return std::make_shared<Multiplication>(right, left->integrate(var));
        }
        
        throw std::runtime_error("Complex multiplication integration not implemented");
    }
    
    ExprPtr Division::integrate(const std::string& var) const override {
        // Try simple cases
        auto denom_const = std::dynamic_pointer_cast<Constant>(right);
        if (denom_const && !right->getVariables().count(var)) {
            // ∫f(x)/c dx = (1/c)*∫f(x) dx
            return std::make_shared<Division>(left->integrate(var), right);
        }
        
        throw std::runtime_error("Complex division integration not implemented");
    }
    
    ExprPtr Power::integrate(const std::string& var) const override {
        auto base_var = std::dynamic_pointer_cast<Variable>(left);
        auto exp_const = std::dynamic_pointer_cast<Constant>(right);
        
        if (base_var && base_var->getName() == var && exp_const) {
            double n = exp_const->getValue();
            if (std::abs(n + 1) < EPSILON) {
                // ∫x^(-1) dx = ln(x)
                return std::make_shared<NaturalLog>(left);
            } else {
                // ∫x^n dx = x^(n+1)/(n+1)
                auto new_exp = std::make_shared<Constant>(n + 1);
                auto x_power = std::make_shared<Power>(left, new_exp);
                return std::make_shared<Division>(x_power, new_exp);
            }
        }
        
        throw std::runtime_error("Complex power integration not implemented");
    }
    
    // ===== SYMBOLIC ALGEBRA SYSTEM =====
    
    class SymbolicAlgebraSystem {
    private:
        std::unordered_map<std::string, ExprPtr> symbol_table;
        
    public:
        // Define a symbol/variable
        void define(const std::string& name, const ExprPtr& expr) {
            symbol_table[name] = expr;
        }
        
        // Get a symbol
        ExprPtr getSymbol(const std::string& name) {
            auto it = symbol_table.find(name);
            if (it != symbol_table.end()) {
                return it->second;
            }
            return std::make_shared<Variable>(name);
        }
        
        // Differentiate expression
        ExprPtr diff(const ExprPtr& expr, const std::string& var) {
            return expr->differentiate(var)->simplify();
        }
        
        // Integrate expression
        ExprPtr integrate(const ExprPtr& expr, const std::string& var) {
            try {
                return expr->integrate(var)->simplify();
            } catch (const std::exception& e) {
                std::cout << "Integration failed: " << e.what() << std::endl;
                return nullptr;
            }
        }
        
        // Simplify expression
        ExprPtr simplify(const ExprPtr& expr) {
            return expr->simplify();
        }
        
        // Substitute variables
        ExprPtr substitute(const ExprPtr& expr, const std::string& var, const ExprPtr& replacement) {
            return expr->substitute(var, replacement)->simplify();
        }
        
        // Expand expression (basic implementation)
        ExprPtr expand(const ExprPtr& expr) {
            // For now, just simplify - full expansion would require more complex logic
            return expr->simplify();
        }
        
        // Factor expression (basic implementation)
        ExprPtr factor(const ExprPtr& expr) {
            // Basic factoring - would need more sophisticated algorithms
            return expr->simplify();
        }
        
        // Solve simple equations (basic implementation)
        std::vector<ExprPtr> solve(const ExprPtr& equation, const std::string& var) {
            // Very basic solver - assumes equation = 0
            std::vector<ExprPtr> solutions;
            
            // For linear equations ax + b = 0, solution is x = -b/a
            if (auto add = std::dynamic_pointer_cast<Addition>(equation)) {
                // Try to extract linear form
                auto left_vars = add->getLeft()->getVariables();
                auto right_vars = add->getRight()->getVariables();
                
                if (left_vars.count(var) && !right_vars.count(var)) {
                    // Form: f(x) + c = 0, so f(x) = -c
                    auto neg_right = std::make_shared<Multiplication>(
                        std::make_shared<Constant>(-1), add->getRight()
                    );
                    
                    if (auto mult = std::dynamic_pointer_cast<Multiplication>(add->getLeft())) {
                        auto coeff = std::dynamic_pointer_cast<Constant>(mult->getLeft());
                        auto var_term = std::dynamic_pointer_cast<Variable>(mult->getRight());
                        if (coeff && var_term && var_term->getName() == var) {
                            // ax + b = 0, so x = -b/a
                            solutions.push_back(std::make_shared<Division>(neg_right, mult->getLeft()));
                        }
                    } else if (auto var_term = std::dynamic_pointer_cast<Variable>(add->getLeft())) {
                        if (var_term->getName() == var) {
                            // x + b = 0, so x = -b
                            solutions.push_back(neg_right);
                        }
                    }
                }
            }
            
            return solutions;
        }
        
        // Evaluate expression with given variables
        double evaluate(const ExprPtr& expr, const VarMap& vars = {}) {
            try {
                return expr->evaluate(vars);
            } catch (const std::exception& e) {
                std::cout << "Evaluation failed: " << e.what() << std::endl;
                return 0.0;
            }
        }
        
        // Calculate Taylor series expansion
        ExprPtr taylorSeries(const ExprPtr& expr, const std::string& var, const ExprPtr& center, int terms) {
            ExprPtr result = std::make_shared<Constant>(0);
            ExprPtr current_derivative = expr;
            ExprPtr factorial_term = std::make_shared<Constant>(1);
            ExprPtr power_term = std::make_shared<Constant>(1);
            
            for (int n = 0; n < terms; ++n) {
                // Evaluate nth derivative at center
                auto substituted = current_derivative->substitute(var, center);
                
                // Add term: f^(n)(center) * (x-center)^n / n!
                auto term = std::make_shared<Division>(
                    std::make_shared<Multiplication>(substituted, power_term),
                    factorial_term
                );
                result = std::make_shared<Addition>(result, term);
                
                // Prepare for next iteration
                if (n + 1 < terms) {
                    current_derivative = current_derivative->differentiate(var);
                    factorial_term = std::make_shared<Multiplication>(
                        factorial_term, 
                        std::make_shared<Constant>(n + 1)
                    );
                    
                    auto x_minus_center = std::make_shared<Subtraction>(
                        std::make_shared<Variable>(var), center
                    );
                    power_term = std::make_shared<Multiplication>(power_term, x_minus_center);
                }
            }
            
            return result->simplify();
        }
        
        // Print all defined symbols
        void printSymbols() const {
            std::cout << "Defined symbols:" << std::endl;
            for (const auto& [name, expr] : symbol_table) {
                std::cout << name << " = " << expr->toString() << std::endl;
            }
        }
    };
    
    // ===== UTILITY FUNCTIONS =====
    
    // Create commonly used expressions
    ExprPtr var(const std::string& name) {
        return std::make_shared<Variable>(name);
    }
    
    ExprPtr constant(double value) {
        return std::make_shared<Constant>(value);
    }
    
    ExprPtr sin(const ExprPtr& expr) {
        return std::make_shared<Sine>(expr);
    }
    
    ExprPtr cos(const ExprPtr& expr) {
        return std::make_shared<Cosine>(expr);
    }
    
    ExprPtr tan(const ExprPtr& expr) {
        return std::make_shared<Tangent>(expr);
    }
    
    ExprPtr exp(const ExprPtr& expr) {
        return std::make_shared<Exponential>(expr);
    }
    
    ExprPtr ln(const ExprPtr& expr) {
        return std::make_shared<NaturalLog>(expr);
    }
    
    ExprPtr sinh(const ExprPtr& expr) {
        return std::make_shared<Sinh>(expr);
    }
    
    ExprPtr cosh(const ExprPtr& expr) {
        return std::make_shared<Cosh>(expr);
    }
    
    ExprPtr pow(const ExprPtr& base, const ExprPtr& exponent) {
        return std::make_shared<Power>(base, exponent);
    }
    
} // namespace SymbolicComputation

// ===== DEMONSTRATION =====

int main() {
    using namespace SymbolicComputation;
    
    std::cout << "=== SYMBOLIC COMPUTATION ENGINE DEMO ===" << std::endl << std::endl;
    
    // Create symbolic algebra system
    SymbolicAlgebraSystem cas;
    
    // ===== BASIC EXPRESSION CONSTRUCTION =====
    std::cout << "1. BASIC EXPRESSION CONSTRUCTION" << std::endl;
    
    auto x = var("x");
    auto y = var("y");
    
    // Build expression: 2*x^2 + 3*x + 1
    auto expr1 = constant(2) * pow(x, constant(2)) + constant(3) * x + constant(1);
    std::cout << "Expression 1: " << expr1->toString() << std::endl;
    
    // Build expression: sin(x) * cos(x)
    auto expr2 = sin(x) * cos(x);
    std::cout << "Expression 2: " << expr2->toString() << std::endl;
    
    // Build expression: e^(x*y) + ln(x)
    auto expr3 = exp(x * y) + ln(x);
    std::cout << "Expression 3: " << expr3->toString() << std::endl;
    std::cout << std::endl;
    
    // ===== SYMBOLIC DIFFERENTIATION =====
    std::cout << "2. SYMBOLIC DIFFERENTIATION" << std::endl;
    
    auto diff_expr1 = cas.diff(expr1, "x");
    std::cout << "d/dx(" << expr1->toString() << ") = " << diff_expr1->toString() << std::endl;
    
    auto diff_expr2 = cas.diff(expr2, "x");
    std::cout << "d/dx(" << expr2->toString() << ") = " << diff_expr2->toString() << std::endl;
    
    auto diff_expr3 = cas.diff(expr3, "x");
    std::cout << "d/dx(" << expr3->toString() << ") = " << diff_expr3->toString() << std::endl;
    
    // Chain rule example
    auto chain_expr = sin(pow(x, constant(2)));
    auto chain_diff = cas.diff(chain_expr, "x");
    std::cout << "d/dx(" << chain_expr->toString() << ") = " << chain_diff->toString() << std::endl;
    std::cout << std::endl;
    
    // ===== SYMBOLIC INTEGRATION =====
    std::cout << "3. SYMBOLIC INTEGRATION" << std::endl;
    
    // Simple polynomial
    auto int_expr1 = pow(x, constant(3)) + constant(2) * pow(x, constant(2)) + x;
    auto integral1 = cas.integrate(int_expr1, "x");
    if (integral1) {
        std::cout << "∫(" << int_expr1->toString() << ")dx = " << integral1->toString() << std::endl;
    }
    
    // Trigonometric functions
    auto sin_x = sin(x);
    auto int_sin = cas.integrate(sin_x, "x");
    if (int_sin) {
        std::cout << "∫(" << sin_x->toString() << ")dx = " << int_sin->toString() << std::endl;
    }
    
    auto cos_x = cos(x);
    auto int_cos = cas.integrate(cos_x, "x");
    if (int_cos) {
        std::cout << "∫(" << cos_x->toString() << ")dx = " << int_cos->toString() << std::endl;
    }
    
    // Exponential function
    auto exp_x = exp(x);
    auto int_exp = cas.integrate(exp_x, "x");
    if (int_exp) {
        std::cout << "∫(" << exp_x->toString() << ")dx = " << int_exp->toString() << std::endl;
    }
    std::cout << std::endl;
    
    // ===== EXPRESSION SIMPLIFICATION =====
    std::cout << "4. EXPRESSION SIMPLIFICATION" << std::endl;
    
    // x + x = 2*x
    auto simple1 = x + x;
    std::cout << simple1->toString() << " simplifies to: " << cas.simplify(simple1)->toString() << std::endl;
    
    // x * 1 = x
    auto simple2 = x * constant(1);
    std::cout << simple2->toString() << " simplifies to: " << cas.simplify(simple2)->toString() << std::endl;
    
    // x^0 = 1
    auto simple3 = pow(x, constant(0));
    std::cout << simple3->toString() << " simplifies to: " << cas.simplify(simple3)->toString() << std::endl;
    
    // sin(0) = 0
    auto simple4 = sin(constant(0));
    std::cout << simple4->toString() << " simplifies to: " << cas.simplify(simple4)->toString() << std::endl;
    
    // exp(ln(x)) = x
    auto simple5 = exp(ln(x));
    std::cout << simple5->toString() << " simplifies to: " << cas.simplify(simple5)->toString() << std::endl;
    std::cout << std::endl;
    
    // ===== SUBSTITUTION =====
    std::cout << "5. SYMBOLIC SUBSTITUTION" << std::endl;
    
    auto sub_expr = pow(x, constant(2)) + constant(2) * x + constant(1);
    auto substituted = cas.substitute(sub_expr, "x", constant(3));
    std::cout << "Substitute x=3 in " << sub_expr->toString() << ": " << substituted->toString() << std::endl;
    
    auto trig_sub = sin(x) + cos(x);
    auto trig_substituted = cas.substitute(trig_sub, "x", constant(PI/2));
    std::cout << "Substitute x=π/2 in " << trig_sub->toString() << ": " << trig_substituted->toString() << std::endl;
    std::cout << std::endl;
    
    // ===== EVALUATION =====
    std::cout << "6. NUMERICAL EVALUATION" << std::endl;
    
    VarMap values = {{"x", 2.0}, {"y", 3.0}};
    
    double val1 = cas.evaluate(expr1, values);
    std::cout << expr1->toString() << " at x=2: " << val1 << std::endl;
    
    double val2 = cas.evaluate(expr2, values);
    std::cout << expr2->toString() << " at x=2: " << val2 << std::endl;
    
    double val3 = cas.evaluate(expr3, values);
    std::cout << expr3->toString() << " at x=2, y=3: " << val3 << std::endl;
    std::cout << std::endl;
    
    // ===== TAYLOR SERIES =====
    std::cout << "7. TAYLOR SERIES EXPANSION" << std::endl;
    
    auto taylor_sin = cas.taylorSeries(sin(x), "x", constant(0), 5);
    std::cout << "Taylor series of sin(x) around x=0 (5 terms): " << taylor_sin->toString() << std::endl;
    
    auto taylor_exp = cas.taylorSeries(exp(x), "x", constant(0), 4);
    std::cout << "Taylor series of exp(x) around x=0 (4 terms): " << taylor_exp->toString() << std::endl;
    std::cout << std::endl;
    
    // ===== EQUATION SOLVING =====
    std::cout << "8. BASIC EQUATION SOLVING" << std::endl;
    
    // Solve 2x + 3 = 0
    auto equation1 = constant(2) * x + constant(3);
    auto solutions1 = cas.solve(equation1, "x");
    if (!solutions1.empty()) {
        std::cout << "Solutions to " << equation1->toString() << " = 0: ";
        for (const auto& sol : solutions1) {
            std::cout << sol->toString() << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    // ===== APPLICATIONS IN SYNTHESIS =====
    std::cout << "9. SYNTHESIZER APPLICATIONS" << std::endl;
    
    std::cout << "🎵 Oscillator Waveforms:" << std::endl;
    auto t = var("t");
    auto freq = var("f");
    auto phase = var("phi");
    
    // Sine wave: A*sin(2πft + φ)
    auto sine_wave = sin(constant(2*PI) * freq * t + phase);
    std::cout << "Sine wave: " << sine_wave->toString() << std::endl;
    
    // FM synthesis: sin(2πfct + β*sin(2πfmt))
    auto fc = var("fc"), fm = var("fm"), beta = var("beta");
    auto fm_wave = sin(constant(2*PI) * fc * t + beta * sin(constant(2*PI) * fm * t));
    std::cout << "FM synthesis: " << fm_wave->toString() << std::endl;
    std::cout << std::endl;
    
    std::cout << "🔧 Filter Transfer Functions:" << std::endl;
    auto s = var("s");
    auto R = var("R"), C = var("C");
    
    // Low-pass filter: H(s) = 1/(RCs + 1)
    auto lpf_transfer = constant(1) / (R * C * s + constant(1));
    std::cout << "Low-pass filter H(s): " << lpf_transfer->toString() << std::endl;
    
    // Differentiate to get phase response
    auto phase_response = cas.diff(ln(lpf_transfer), "s");
    std::cout << "Phase response derivative: " << phase_response->toString() << std::endl;
    std::cout << std::endl;
    
    std::cout << "📐 Envelope Generators:" << std::endl;
    // Exponential decay: A*exp(-t/τ)
    auto A = var("A"), tau = var("tau");
    auto exp_decay = A * exp(constant(-1) * t / tau);
    std::cout << "Exponential decay: " << exp_decay->toString() << std::endl;
    
    // Derivative gives slope
    auto decay_slope = cas.diff(exp_decay, "t");
    std::cout << "Decay slope: " << decay_slope->toString() << std::endl;
    std::cout << std::endl;
    
    // ===== PERFORMANCE CHARACTERISTICS =====
    std::cout << "10. SYMBOLIC ENGINE CAPABILITIES" << std::endl;
    
    std::cout << "✅ Supported Operations:" << std::endl;
    std::cout << "   • Basic arithmetic: +, -, *, /, ^" << std::endl;
    std::cout << "   • Trigonometric: sin, cos, tan" << std::endl;
    std::cout << "   • Hyperbolic: sinh, cosh" << std::endl;
    std::cout << "   • Exponential & Logarithmic: exp, ln" << std::endl;
    std::cout << "   • Automatic simplification" << std::endl;
    std::cout << "   • Symbolic differentiation" << std::endl;
    std::cout << "   • Basic symbolic integration" << std::endl;
    std::cout << "   • Variable substitution" << std::endl;
    std::cout << "   • Taylor series expansion" << std::endl;
    std::cout << "   • Simple equation solving" << std::endl;
    std::cout << std::endl;
    
    std::cout << "🚀 Perfect for Audio Applications:" << std::endl;
    std::cout << "   • Derive filter transfer functions symbolically" << std::endl;
    std::cout << "   • Analyze oscillator harmonics and distortion" << std::endl;
    std::cout << "   • Calculate optimal component values" << std::endl;
    std::cout << "   • Model nonlinear circuit behavior" << std::endl;
    std::cout << "   • Design custom waveform generators" << std::endl;
    std::cout << "   • Optimize DSP algorithm parameters" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== SYMBOLIC COMPUTATION ENGINE COMPLETE ===" << std::endl;
    std::cout << "🧠 Ready for advanced mathematical reasoning in synthesis! 🧠" << std::endl;
    
    return 0;
}prPtr operator+(const ExprPtr& other) const;
        ExprPtr operator-(const ExprPtr& other) const;
        ExprPtr operator*(const ExprPtr& other) const;
        ExprPtr operator/(const ExprPtr& other) const;
        ExprPtr operator^(const ExprPtr& other) const;
    };
    
    // ===== CONSTANT EXPRESSIONS =====
    
    class Constant : public Expression {
    private:
        double value;
        
    public:
        Constant(double v) : value(v) {}
        
        ExprPtr clone() const override {
            return std::make_shared<Constant>(value);
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            return std::make_shared<Constant>(0.0);
        }
        
        ExprPtr integrate(const std::string& var) const override;
        
        ExprPtr simplify() const override {
            return clone();
        }
        
        std::string toString() const override {
            if (std::abs(value - std::round(value)) < EPSILON) {
                return std::to_string((int)std::round(value));
            }
            return std::to_string(value);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return value;
        }
        
        std::set<std::string> getVariables() const override {
            return {};
        }
        
        bool equals(const ExprPtr& other) const override {
            auto const_other = std::dynamic_pointer_cast<Constant>(other);
            return const_other && std::abs(value - const_other->value) < EPSILON;
        }
        
        ExprPtr substitute(const std::string& var, const ExprPtr& replacement) const override {
            return clone();
        }
        
        bool isZero() const override { return std::abs(value) < EPSILON; }
        bool isOne() const override { return std::abs(value - 1.0) < EPSILON; }
        
        double getValue() const { return value; }
        int getPrecedence() const override { return 100; }
    };
    
    // ===== VARIABLE EXPRESSIONS =====
    
    class Variable : public Expression {
    private:
        std::string name;
        
    public:
        Variable(const std::string& n) : name(n) {}
        
        ExprPtr clone() const override {
            return std::make_shared<Variable>(name);
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            return std::make_shared<Constant>(name == var ? 1.0 : 0.0);
        }
        
        ExprPtr integrate(const std::string& var) const override;
        
        ExprPtr simplify() const override {
            return clone();
        }
        
        std::string toString() const override {
            return name;
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            auto it = vars.find(name);
            if (it != vars.end()) {
                return it->second;
            }
            throw std::runtime_error("Variable '" + name + "' not found in evaluation context");
        }
        
        std::set<std::string> getVariables() const override {
            return {name};
        }
        
        bool equals(const ExprPtr& other) const override {
            auto var_other = std::dynamic_pointer_cast<Variable>(other);
            return var_other && name == var_other->name;
        }
        
        ExprPtr substitute(const std::string& var, const ExprPtr& replacement) const override {
            return (name == var) ? replacement : clone();
        }
        
        const std::string& getName() const { return name; }
        int getPrecedence() const override { return 100; }
    };
    
    // ===== BINARY OPERATION BASE CLASS =====
    
    class BinaryOperation : public Expression {
    protected:
        ExprPtr left, right;
        std::string op_symbol;
        int precedence;
        
    public:
        BinaryOperation(ExprPtr l, ExprPtr r, const std::string& op, int prec) 
            : left(l), right(r), op_symbol(op), precedence(prec) {}
        
        std::string toString() const override {
            std::string left_str = left->toString();
            std::string right_str = right->toString();
            
            // Add parentheses based on precedence
            if (left->getPrecedence() < precedence) {
                left_str = "(" + left_str + ")";
            }
            if (right->getPrecedence() < precedence || 
                (right->getPrecedence() == precedence && (op_symbol == "-" || op_symbol == "/"))) {
                right_str = "(" + right_str + ")";
            }
            
            return left_str + " " + op_symbol + " " + right_str;
        }
        
        std::set<std::string> getVariables() const override {
            auto left_vars = left->getVariables();
            auto right_vars = right->getVariables();
            left_vars.insert(right_vars.begin(), right_vars.end());
            return left_vars;
        }
        
        ExprPtr substitute(const std::string& var, const ExprPtr& replacement) const override {
            return createBinaryOp(left->substitute(var, replacement), 
                                 right->substitute(var, replacement));
        }
        
        int getPrecedence() const override { return precedence; }
        
        const ExprPtr& getLeft() const { return left; }
        const ExprPtr& getRight() const { return right; }
        
    protected:
        virtual ExprPtr createBinaryOp(const ExprPtr& l, const ExprPtr& r) const = 0;
    };
    
    // ===== ADDITION =====
    
    class Addition : public BinaryOperation {
    public:
        Addition(ExprPtr l, ExprPtr r) : BinaryOperation(l, r, "+", 10) {}
        
        ExprPtr clone() const override {
            return std::make_shared<Addition>(left->clone(), right->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            return std::make_shared<Addition>(left->differentiate(var), right->differentiate(var));
        }
        
        ExprPtr integrate(const std::string& var) const override {
            return std::make_shared<Addition>(left->integrate(var), right->integrate(var));
        }
        
        ExprPtr simplify() const override {
            auto l_simp = left->simplify();
            auto r_simp = right->simplify();
            
            // 0 + x = x, x + 0 = x
            if (l_simp->isZero()) return r_simp;
            if (r_simp->isZero()) return l_simp;
            
            // Combine constants
            auto l_const = std::dynamic_pointer_cast<Constant>(l_simp);
            auto r_const = std::dynamic_pointer_cast<Constant>(r_simp);
            if (l_const && r_const) {
                return std::make_shared<Constant>(l_const->getValue() + r_const->getValue());
            }
            
            // x + x = 2*x
            if (l_simp->equals(r_simp)) {
                return std::make_shared<Multiplication>(std::make_shared<Constant>(2), l_simp);
            }
            
            return std::make_shared<Addition>(l_simp, r_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return left->evaluate(vars) + right->evaluate(vars);
        }
        
        bool equals(const ExprPtr& other) const override {
            auto add_other = std::dynamic_pointer_cast<Addition>(other);
            return add_other && 
                   ((left->equals(add_other->left) && right->equals(add_other->right)) ||
                    (left->equals(add_other->right) && right->equals(add_other->left))); // Commutative
        }
        
    protected:
        ExprPtr createBinaryOp(const ExprPtr& l, const ExprPtr& r) const override {
            return std::make_shared<Addition>(l, r);
        }
    };
    
    // ===== SUBTRACTION =====
    
    class Subtraction : public BinaryOperation {
    public:
        Subtraction(ExprPtr l, ExprPtr r) : BinaryOperation(l, r, "-", 10) {}
        
        ExprPtr clone() const override {
            return std::make_shared<Subtraction>(left->clone(), right->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            return std::make_shared<Subtraction>(left->differentiate(var), right->differentiate(var));
        }
        
        ExprPtr integrate(const std::string& var) const override {
            return std::make_shared<Subtraction>(left->integrate(var), right->integrate(var));
        }
        
        ExprPtr simplify() const override {
            auto l_simp = left->simplify();
            auto r_simp = right->simplify();
            
            // x - 0 = x
            if (r_simp->isZero()) return l_simp;
            // 0 - x = -x
            if (l_simp->isZero()) return std::make_shared<Multiplication>(std::make_shared<Constant>(-1), r_simp);
            
            // Combine constants
            auto l_const = std::dynamic_pointer_cast<Constant>(l_simp);
            auto r_const = std::dynamic_pointer_cast<Constant>(r_simp);
            if (l_const && r_const) {
                return std::make_shared<Constant>(l_const->getValue() - r_const->getValue());
            }
            
            // x - x = 0
            if (l_simp->equals(r_simp)) {
                return std::make_shared<Constant>(0);
            }
            
            return std::make_shared<Subtraction>(l_simp, r_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return left->evaluate(vars) - right->evaluate(vars);
        }
        
        bool equals(const ExprPtr& other) const override {
            auto sub_other = std::dynamic_pointer_cast<Subtraction>(other);
            return sub_other && left->equals(sub_other->left) && right->equals(sub_other->right);
        }
        
    protected:
        ExprPtr createBinaryOp(const ExprPtr& l, const ExprPtr& r) const override {
            return std::make_shared<Subtraction>(l, r);
        }
    };
    
    // ===== MULTIPLICATION =====
    
    class Multiplication : public BinaryOperation {
    public:
        Multiplication(ExprPtr l, ExprPtr r) : BinaryOperation(l, r, "*", 20) {}
        
        ExprPtr clone() const override {
            return std::make_shared<Multiplication>(left->clone(), right->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            // Product rule: (f*g)' = f'*g + f*g'
            return std::make_shared<Addition>(
                std::make_shared<Multiplication>(left->differentiate(var), right),
                std::make_shared<Multiplication>(left, right->differentiate(var))
            );
        }
        
        ExprPtr integrate(const std::string& var) const override;
        
        ExprPtr simplify() const override {
            auto l_simp = left->simplify();
            auto r_simp = right->simplify();
            
            // 0 * x = 0, x * 0 = 0
            if (l_simp->isZero() || r_simp->isZero()) {
                return std::make_shared<Constant>(0);
            }
            
            // 1 * x = x, x * 1 = x
            if (l_simp->isOne()) return r_simp;
            if (r_simp->isOne()) return l_simp;
            
            // Combine constants
            auto l_const = std::dynamic_pointer_cast<Constant>(l_simp);
            auto r_const = std::dynamic_pointer_cast<Constant>(r_simp);
            if (l_const && r_const) {
                return std::make_shared<Constant>(l_const->getValue() * r_const->getValue());
            }
            
            // x * x = x^2
            if (l_simp->equals(r_simp)) {
                return std::make_shared<Power>(l_simp, std::make_shared<Constant>(2));
            }
            
            return std::make_shared<Multiplication>(l_simp, r_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return left->evaluate(vars) * right->evaluate(vars);
        }
        
        bool equals(const ExprPtr& other) const override {
            auto mult_other = std::dynamic_pointer_cast<Multiplication>(other);
            return mult_other && 
                   ((left->equals(mult_other->left) && right->equals(mult_other->right)) ||
                    (left->equals(mult_other->right) && right->equals(mult_other->left))); // Commutative
        }
        
        std::string toString() const override {
            std::string left_str = left->toString();
            std::string right_str = right->toString();
            
            // Add parentheses based on precedence
            if (left->getPrecedence() < precedence) {
                left_str = "(" + left_str + ")";
            }
            if (right->getPrecedence() < precedence) {
                right_str = "(" + right_str + ")";
            }
            
            // Implicit multiplication for variables and functions
            auto l_const = std::dynamic_pointer_cast<Constant>(left);
            if (l_const || std::dynamic_pointer_cast<Variable>(right) || 
                std::dynamic_pointer_cast<Function>(right)) {
                return left_str + right_str;
            }
            
            return left_str + "*" + right_str;
        }
        
    protected:
        ExprPtr createBinaryOp(const ExprPtr& l, const ExprPtr& r) const override {
            return std::make_shared<Multiplication>(l, r);
        }
    };
    
    // ===== DIVISION =====
    
    class Division : public BinaryOperation {
    public:
        Division(ExprPtr l, ExprPtr r) : BinaryOperation(l, r, "/", 20) {}
        
        ExprPtr clone() const override {
            return std::make_shared<Division>(left->clone(), right->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            // Quotient rule: (f/g)' = (f'*g - f*g')/g^2
            auto numerator = std::make_shared<Subtraction>(
                std::make_shared<Multiplication>(left->differentiate(var), right),
                std::make_shared<Multiplication>(left, right->differentiate(var))
            );
            auto denominator = std::make_shared<Power>(right, std::make_shared<Constant>(2));
            return std::make_shared<Division>(numerator, denominator);
        }
        
        ExprPtr integrate(const std::string& var) const override;
        
        ExprPtr simplify() const override {
            auto l_simp = left->simplify();
            auto r_simp = right->simplify();
            
            // 0 / x = 0
            if (l_simp->isZero()) {
                return std::make_shared<Constant>(0);
            }
            
            // x / 1 = x
            if (r_simp->isOne()) return l_simp;
            
            // x / x = 1
            if (l_simp->equals(r_simp)) {
                return std::make_shared<Constant>(1);
            }
            
            // Combine constants
            auto l_const = std::dynamic_pointer_cast<Constant>(l_simp);
            auto r_const = std::dynamic_pointer_cast<Constant>(r_simp);
            if (l_const && r_const) {
                if (std::abs(r_const->getValue()) < EPSILON) {
                    throw std::runtime_error("Division by zero");
                }
                return std::make_shared<Constant>(l_const->getValue() / r_const->getValue());
            }
            
            return std::make_shared<Division>(l_simp, r_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            double denom = right->evaluate(vars);
            if (std::abs(denom) < EPSILON) {
                throw std::runtime_error("Division by zero in evaluation");
            }
            return left->evaluate(vars) / denom;
        }
        
        bool equals(const ExprPtr& other) const override {
            auto div_other = std::dynamic_pointer_cast<Division>(other);
            return div_other && left->equals(div_other->left) && right->equals(div_other->right);
        }
        
    protected:
        ExprPtr createBinaryOp(const ExprPtr& l, const ExprPtr& r) const override {
            return std::make_shared<Division>(l, r);
        }
    };
    
    // ===== POWER =====
    
    class Power : public BinaryOperation {
    public:
        Power(ExprPtr l, ExprPtr r) : BinaryOperation(l, r, "^", 30) {}
        
        ExprPtr clone() const override {
            return std::make_shared<Power>(left->clone(), right->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            auto r_const = std::dynamic_pointer_cast<Constant>(right);
            if (r_const) {
                // Power rule: (f^n)' = n*f^(n-1)*f'
                double n = r_const->getValue();
                if (std::abs(n) < EPSILON) {
                    return std::make_shared<Constant>(0);
                }
                auto new_power = std::make_shared<Power>(left, std::make_shared<Constant>(n - 1));
                return std::make_shared<Multiplication>(
                    std::make_shared<Multiplication>(std::make_shared<Constant>(n), new_power),
                    left->differentiate(var)
                );
            }
            
            // General case: (f^g)' = f^g * (g'*ln(f) + g*f'/f)
            auto ln_f = std::make_shared<NaturalLog>(left);
            auto g_prime_ln_f = std::make_shared<Multiplication>(right->differentiate(var), ln_f);
            auto g_f_prime_over_f = std::make_shared<Multiplication>(
                right, 
                std::make_shared<Division>(left->differentiate(var), left)
            );
            auto derivative_factor = std::make_shared<Addition>(g_prime_ln_f, g_f_prime_over_f);
            
            return std::make_shared<Multiplication>(clone(), derivative_factor);
        }
        
        ExprPtr integrate(const std::string& var) const override;
        
        ExprPtr simplify() const override {
            auto l_simp = left->simplify();
            auto r_simp = right->simplify();
            
            // x^0 = 1
            if (r_simp->isZero()) {
                return std::make_shared<Constant>(1);
            }
            
            // x^1 = x
            if (r_simp->isOne()) return l_simp;
            
            // 0^x = 0 (for x != 0)
            if (l_simp->isZero() && !r_simp->isZero()) {
                return std::make_shared<Constant>(0);
            }
            
            // 1^x = 1
            if (l_simp->isOne()) {
                return std::make_shared<Constant>(1);
            }
            
            // Combine constants
            auto l_const = std::dynamic_pointer_cast<Constant>(l_simp);
            auto r_const = std::dynamic_pointer_cast<Constant>(r_simp);
            if (l_const && r_const) {
                return std::make_shared<Constant>(std::pow(l_const->getValue(), r_const->getValue()));
            }
            
            return std::make_shared<Power>(l_simp, r_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return std::pow(left->evaluate(vars), right->evaluate(vars));
        }
        
        bool equals(const ExprPtr& other) const override {
            auto pow_other = std::dynamic_pointer_cast<Power>(other);
            return pow_other && left->equals(pow_other->left) && right->equals(pow_other->right);
        }
        
    protected:
        ExprPtr createBinaryOp(const ExprPtr& l, const ExprPtr& r) const override {
            return std::make_shared<Power>(l, r);
        }
    };
    
    // ===== FUNCTION BASE CLASS =====
    
    class Function : public Expression {
    protected:
        ExprPtr argument;
        std::string func_name;
        
    public:
        Function(ExprPtr arg, const std::string& name) : argument(arg), func_name(name) {}
        
        std::string toString() const override {
            return func_name + "(" + argument->toString() + ")";
        }
        
        std::set<std::string> getVariables() const override {
            return argument->getVariables();
        }
        
        ExprPtr substitute(const std::string& var, const ExprPtr& replacement) const override {
            return createFunction(argument->substitute(var, replacement));
        }
        
        bool equals(const ExprPtr& other) const override {
            auto func_other = std::dynamic_pointer_cast<Function>(other);
            return func_other && func_name == func_other->func_name && argument->equals(func_other->argument);
        }
        
        int getPrecedence() const override { return 100; }
        
        const ExprPtr& getArgument() const { return argument; }
        
    protected:
        virtual ExprPtr createFunction(const ExprPtr& arg) const = 0;
    };
    
    // ===== TRIGONOMETRIC FUNCTIONS =====
    
    class Sine : public Function {
    public:
        Sine(ExprPtr arg) : Function(arg, "sin") {}
        
        ExprPtr clone() const override {
            return std::make_shared<Sine>(argument->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            // (sin(f))' = cos(f) * f'
            return std::make_shared<Multiplication>(
                std::make_shared<Cosine>(argument),
                argument->differentiate(var)
            );
        }
        
        ExprPtr integrate(const std::string& var) const override {
            // ∫sin(x)dx = -cos(x)
            if (auto var_arg = std::dynamic_pointer_cast<Variable>(argument)) {
                if (var_arg->getName() == var) {
                    return std::make_shared<Multiplication>(
                        std::make_shared<Constant>(-1),
                        std::make_shared<Cosine>(argument)
                    );
                }
            }
            throw std::runtime_error("Integration of sin with complex argument not implemented");
        }
        
        ExprPtr simplify() const override {
            auto arg_simp = argument->simplify();
            
            // sin(0) = 0
            if (arg_simp->isZero()) {
                return std::make_shared<Constant>(0);
            }
            
            return std::make_shared<Sine>(arg_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return std::sin(argument->evaluate(vars));
        }
        
    protected:
        ExprPtr createFunction(const ExprPtr& arg) const override {
            return std::make_shared<Sine>(arg);
        }
    };
    
    class Cosine : public Function {
    public:
        Cosine(ExprPtr arg) : Function(arg, "cos") {}
        
        ExprPtr clone() const override {
            return std::make_shared<Cosine>(argument->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            // (cos(f))' = -sin(f) * f'
            return std::make_shared<Multiplication>(
                std::make_shared<Multiplication>(std::make_shared<Constant>(-1), std::make_shared<Sine>(argument)),
                argument->differentiate(var)
            );
        }
        
        ExprPtr integrate(const std::string& var) const override {
            // ∫cos(x)dx = sin(x)
            if (auto var_arg = std::dynamic_pointer_cast<Variable>(argument)) {
                if (var_arg->getName() == var) {
                    return std::make_shared<Sine>(argument);
                }
            }
            throw std::runtime_error("Integration of cos with complex argument not implemented");
        }
        
        ExprPtr simplify() const override {
            auto arg_simp = argument->simplify();
            
            // cos(0) = 1
            if (arg_simp->isZero()) {
                return std::make_shared<Constant>(1);
            }
            
            return std::make_shared<Cosine>(arg_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return std::cos(argument->evaluate(vars));
        }
        
    protected:
        ExprPtr createFunction(const ExprPtr& arg) const override {
            return std::make_shared<Cosine>(arg);
        }
    };
    
    class Tangent : public Function {
    public:
        Tangent(ExprPtr arg) : Function(arg, "tan") {}
        
        ExprPtr clone() const override {
            return std::make_shared<Tangent>(argument->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            // (tan(f))' = sec²(f) * f' = (1/cos²(f)) * f'
            auto cos_f = std::make_shared<Cosine>(argument);
            auto cos2_f = std::make_shared<Power>(cos_f, std::make_shared<Constant>(2));
            return std::make_shared<Multiplication>(
                std::make_shared<Division>(std::make_shared<Constant>(1), cos2_f),
                argument->differentiate(var)
            );
        }
        
        ExprPtr integrate(const std::string& var) const override {
            // ∫tan(x)dx = -ln|cos(x)|
            if (auto var_arg = std::dynamic_pointer_cast<Variable>(argument)) {
                if (var_arg->getName() == var) {
                    return std::make_shared<Multiplication>(
                        std::make_shared<Constant>(-1),
                        std::make_shared<NaturalLog>(std::make_shared<Cosine>(argument))
                    );
                }
            }
            throw std::runtime_error("Integration of tan with complex argument not implemented");
        }
        
        ExprPtr simplify() const override {
            auto arg_simp = argument->simplify();
            
            // tan(0) = 0
            if (arg_simp->isZero()) {
                return std::make_shared<Constant>(0);
            }
            
            return std::make_shared<Tangent>(arg_simp);
        }
        
        double evaluate(const VarMap& vars = {}) const override {
            return std::tan(argument->evaluate(vars));
        }
        
    protected:
        ExprPtr createFunction(const ExprPtr& arg) const override {
            return std::make_shared<Tangent>(arg);
        }
    };
    
    // ===== EXPONENTIAL AND LOGARITHMIC FUNCTIONS =====
    
    class Exponential : public Function {
    public:
        Exponential(ExprPtr arg) : Function(arg, "exp") {}
        
        ExprPtr clone() const override {
            return std::make_shared<Exponential>(argument->clone());
        }
        
        ExprPtr differentiate(const std::string& var) const override {
            // (e^f)' = e^f * f'
            return std::make_shared<Multiplication>(clone(), argument->differentiate(var));
        }
        
        ExprPtr integrate(const std::string& var) const override {
            // ∫e^x dx = e^x
            if (auto var_arg = std::dynamic_pointer_cast<Variable>(argument)) {
                if (var_arg->getName() == var) {
                    return clone();
                }
            }
            throw std::runtime_error("Integration of exp with complex argument not implemented");
        }
        
        Ex