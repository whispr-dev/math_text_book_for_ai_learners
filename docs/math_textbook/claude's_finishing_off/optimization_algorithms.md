#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>
#include <random>
#include <limits>
#include <memory>
#include <map>
#include <queue>

namespace OptimizationAlgorithms {

    const double EPSILON = 1e-12;
    const double INF = std::numeric_limits<double>::infinity();
    
    // ===== OPTIMIZATION RESULT STRUCTURES =====
    
    struct OptimizationResult {
        std::vector<double> solution;
        double objective_value;
        int iterations;
        bool converged;
        std::string method;
        double gradient_norm;
        
        void print() const {
            std::cout << "Method: " << method << std::endl;
            std::cout << "Converged: " << (converged ? "Yes" : "No") << std::endl;
            std::cout << "Iterations: " << iterations << std::endl;
            std::cout << "Objective: " << objective_value << std::endl;
            std::cout << "Solution: [";
            for (size_t i = 0; i < solution.size(); ++i) {
                std::cout << solution[i];
                if (i < solution.size() - 1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
            if (gradient_norm >= 0) {
                std::cout << "Gradient norm: " << gradient_norm << std::endl;
            }
            std::cout << std::endl;
        }
    };
    
    struct LinearProgramResult {
        std::vector<double> solution;
        double objective_value;
        bool optimal;
        bool unbounded;
        bool infeasible;
        int iterations;
        
        void print() const {
            std::cout << "Linear Program Result:" << std::endl;
            if (optimal) {
                std::cout << "Status: OPTIMAL" << std::endl;
                std::cout << "Objective: " << objective_value << std::endl;
                std::cout << "Solution: [";
                for (size_t i = 0; i < solution.size(); ++i) {
                    std::cout << solution[i];
                    if (i < solution.size() - 1) std::cout << ", ";
                }
                std::cout << "]" << std::endl;
            } else if (unbounded) {
                std::cout << "Status: UNBOUNDED" << std::endl;
            } else if (infeasible) {
                std::cout << "Status: INFEASIBLE" << std::endl;
            }
            std::cout << "Iterations: " << iterations << std::endl << std::endl;
        }
    };
    
    // ===== UTILITY FUNCTIONS =====
    
    // Vector operations
    std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b) {
        std::vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }
    
    std::vector<double> vectorScale(const std::vector<double>& v, double scalar) {
        std::vector<double> result(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            result[i] = v[i] * scalar;
        }
        return result;
    }
    
    double vectorNorm(const std::vector<double>& v) {
        double sum = 0.0;
        for (double x : v) {
            sum += x * x;
        }
        return std::sqrt(sum);
    }
    
    double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            result += a[i] * b[i];
        }
        return result;
    }
    
    // Numerical gradient computation
    std::vector<double> computeGradient(const std::function<double(const std::vector<double>&)>& f,
                                       const std::vector<double>& x, double h = 1e-8) {
        std::vector<double> gradient(x.size());
        
        for (size_t i = 0; i < x.size(); ++i) {
            std::vector<double> x_plus = x, x_minus = x;
            x_plus[i] += h;
            x_minus[i] -= h;
            gradient[i] = (f(x_plus) - f(x_minus)) / (2.0 * h);
        }
        
        return gradient;
    }
    
    // ===== LINEAR PROGRAMMING - SIMPLEX METHOD =====
    
    class SimplexSolver {
    private:
        std::vector<std::vector<double>> tableau;
        std::vector<int> basic_vars;
        int m, n; // m constraints, n variables
        
        void pivot(int pivot_row, int pivot_col) {
            // Normalize pivot row
            double pivot_element = tableau[pivot_row][pivot_col];
            for (int j = 0; j <= n; ++j) {
                tableau[pivot_row][j] /= pivot_element;
            }
            
            // Eliminate other rows
            for (int i = 0; i <= m; ++i) {
                if (i != pivot_row && std::abs(tableau[i][pivot_col]) > EPSILON) {
                    double multiplier = tableau[i][pivot_col];
                    for (int j = 0; j <= n; ++j) {
                        tableau[i][j] -= multiplier * tableau[pivot_row][j];
                    }
                }
            }
            
            basic_vars[pivot_row] = pivot_col;
        }
        
        int findEnteringVariable() {
            int entering = -1;
            double min_ratio = 0.0;
            
            for (int j = 0; j < n; ++j) {
                if (tableau[0][j] < min_ratio) {
                    min_ratio = tableau[0][j];
                    entering = j;
                }
            }
            
            return entering;
        }
        
        int findLeavingVariable(int entering_var) {
            int leaving = -1;
            double min_ratio = INF;
            
            for (int i = 1; i <= m; ++i) {
                if (tableau[i][entering_var] > EPSILON) {
                    double ratio = tableau[i][n] / tableau[i][entering_var];
                    if (ratio < min_ratio) {
                        min_ratio = ratio;
                        leaving = i;
                    }
                }
            }
            
            return leaving;
        }
        
        bool isOptimal() {
            for (int j = 0; j < n; ++j) {
                if (tableau[0][j] < -EPSILON) {
                    return false;
                }
            }
            return true;
        }
        
        bool isUnbounded(int entering_var) {
            for (int i = 1; i <= m; ++i) {
                if (tableau[i][entering_var] > EPSILON) {
                    return false;
                }
            }
            return true;
        }
        
    public:
        // Solve: minimize c^T x subject to Ax <= b, x >= 0
        LinearProgramResult solve(const std::vector<double>& c,
                                const std::vector<std::vector<double>>& A,
                                const std::vector<double>& b) {
            m = A.size();    // number of constraints
            n = c.size();    // number of variables
            
            // Initialize tableau with slack variables
            tableau.assign(m + 1, std::vector<double>(n + m + 1, 0.0));
            basic_vars.assign(m, 0);
            
            // Set up objective function (row 0)
            for (int j = 0; j < n; ++j) {
                tableau[0][j] = c[j]; // For minimization
            }
            
            // Set up constraints (rows 1 to m)
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    tableau[i + 1][j] = A[i][j];
                }
                tableau[i + 1][n + i] = 1.0; // Slack variable
                tableau[i + 1][n + m] = b[i]; // RHS
                basic_vars[i] = n + i; // Slack variables are initially basic
            }
            
            LinearProgramResult result;
            result.iterations = 0;
            result.optimal = false;
            result.unbounded = false;
            result.infeasible = false;
            
            const int MAX_ITERATIONS = 1000;
            
            while (!isOptimal() && result.iterations < MAX_ITERATIONS) {
                int entering_var = findEnteringVariable();
                if (entering_var == -1) break; // Should not happen if not optimal
                
                if (isUnbounded(entering_var)) {
                    result.unbounded = true;
                    return result;
                }
                
                int leaving_var = findLeavingVariable(entering_var);
                if (leaving_var == -1) {
                    result.unbounded = true;
                    return result;
                }
                
                pivot(leaving_var, entering_var);
                result.iterations++;
            }
            
            if (isOptimal()) {
                result.optimal = true;
                result.objective_value = tableau[0][n + m];
                
                // Extract solution
                result.solution.assign(n, 0.0);
                for (int i = 0; i < m; ++i) {
                    if (basic_vars[i] < n) {
                        result.solution[basic_vars[i]] = tableau[i + 1][n + m];
                    }
                }
            } else {
                result.infeasible = true;
            }
            
            return result;
        }
    };
    
    // ===== GRADIENT DESCENT =====
    
    OptimizationResult gradientDescent(const std::function<double(const std::vector<double>&)>& objective,
                                     const std::vector<double>& initial_point,
                                     double learning_rate = 0.01,
                                     double tolerance = 1e-6,
                                     int max_iterations = 1000) {
        std::vector<double> x = initial_point;
        OptimizationResult result;
        result.method = "Gradient Descent";
        
        for (int iter = 0; iter < max_iterations; ++iter) {
            std::vector<double> gradient = computeGradient(objective, x);
            double grad_norm = vectorNorm(gradient);
            
            if (grad_norm < tolerance) {
                result.converged = true;
                result.gradient_norm = grad_norm;
                break;
            }
            
            // Update: x = x - α * ∇f(x)
            std::vector<double> step = vectorScale(gradient, -learning_rate);
            x = vectorAdd(x, step);
        }
        
        result.solution = x;
        result.objective_value = objective(x);
        result.iterations = max_iterations;
        if (!result.converged) {
            result.gradient_norm = vectorNorm(computeGradient(objective, x));
        }
        
        return result;
    }
    
    // ===== NEWTON'S METHOD =====
    
    // Approximate Hessian using finite differences
    std::vector<std::vector<double>> computeHessian(const std::function<double(const std::vector<double>&)>& f,
                                                   const std::vector<double>& x, double h = 1e-6) {
        int n = x.size();
        std::vector<std::vector<double>> hessian(n, std::vector<double>(n));
        
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::vector<double> x_pp = x, x_pm = x, x_mp = x, x_mm = x;
                x_pp[i] += h; x_pp[j] += h;
                x_pm[i] += h; x_pm[j] -= h;
                x_mp[i] -= h; x_mp[j] += h;
                x_mm[i] -= h; x_mm[j] -= h;
                
                hessian[i][j] = (f(x_pp) - f(x_pm) - f(x_mp) + f(x_mm)) / (4.0 * h * h);
            }
        }
        
        return hessian;
    }
    
    // Simple matrix inversion using Gaussian elimination
    std::vector<std::vector<double>> invertMatrix(std::vector<std::vector<double>> A) {
        int n = A.size();
        std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));
        
        // Initialize identity matrix
        for (int i = 0; i < n; ++i) {
            inv[i][i] = 1.0;
        }
        
        // Gaussian elimination
        for (int i = 0; i < n; ++i) {
            // Find pivot
            int max_row = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                    max_row = k;
                }
            }
            
            if (std::abs(A[max_row][i]) < EPSILON) {
                throw std::runtime_error("Matrix is singular");
            }
            
            if (max_row != i) {
                std::swap(A[i], A[max_row]);
                std::swap(inv[i], inv[max_row]);
            }
            
            // Scale pivot row
            double pivot = A[i][i];
            for (int j = 0; j < n; ++j) {
                A[i][j] /= pivot;
                inv[i][j] /= pivot;
            }
            
            // Eliminate column
            for (int k = 0; k < n; ++k) {
                if (k != i && std::abs(A[k][i]) > EPSILON) {
                    double factor = A[k][i];
                    for (int j = 0; j < n; ++j) {
                        A[k][j] -= factor * A[i][j];
                        inv[k][j] -= factor * inv[i][j];
                    }
                }
            }
        }
        
        return inv;
    }
    
    OptimizationResult newtonMethod(const std::function<double(const std::vector<double>&)>& objective,
                                  const std::vector<double>& initial_point,
                                  double tolerance = 1e-6,
                                  int max_iterations = 100) {
        std::vector<double> x = initial_point;
        OptimizationResult result;
        result.method = "Newton's Method";
        
        for (int iter = 0; iter < max_iterations; ++iter) {
            std::vector<double> gradient = computeGradient(objective, x);
            double grad_norm = vectorNorm(gradient);
            
            if (grad_norm < tolerance) {
                result.converged = true;
                result.gradient_norm = grad_norm;
                break;
            }
            
            try {
                auto hessian = computeHessian(objective, x);
                auto hessian_inv = invertMatrix(hessian);
                
                // Newton step: x = x - H^(-1) * ∇f(x)
                std::vector<double> newton_step(x.size(), 0.0);
                for (size_t i = 0; i < x.size(); ++i) {
                    for (size_t j = 0; j < x.size(); ++j) {
                        newton_step[i] -= hessian_inv[i][j] * gradient[j];
                    }
                }
                
                x = vectorAdd(x, newton_step);
                
            } catch (const std::exception&) {
                // Fallback to gradient descent if Hessian is singular
                std::vector<double> step = vectorScale(gradient, -0.01);
                x = vectorAdd(x, step);
            }
        }
        
        result.solution = x;
        result.objective_value = objective(x);
        result.iterations = max_iterations;
        if (!result.converged) {
            result.gradient_norm = vectorNorm(computeGradient(objective, x));
        }
        
        return result;
    }
    
    // ===== SIMULATED ANNEALING =====
    
    OptimizationResult simulatedAnnealing(const std::function<double(const std::vector<double>&)>& objective,
                                        const std::vector<double>& initial_point,
                                        const std::vector<std::pair<double, double>>& bounds,
                                        double initial_temperature = 100.0,
                                        double cooling_rate = 0.95,
                                        int iterations_per_temp = 100,
                                        int max_iterations = 10000) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        
        std::vector<double> current = initial_point;
        std::vector<double> best = current;
        double current_value = objective(current);
        double best_value = current_value;
        double temperature = initial_temperature;
        
        OptimizationResult result;
        result.method = "Simulated Annealing";
        result.converged = false;
        
        for (int iter = 0; iter < max_iterations; ++iter) {
            // Generate neighbor
            std::vector<double> neighbor = current;
            for (size_t i = 0; i < neighbor.size(); ++i) {
                std::normal_distribution<double> normal(0.0, temperature * 0.01);
                neighbor[i] += normal(gen);
                
                // Apply bounds
                if (i < bounds.size()) {
                    neighbor[i] = std::max(bounds[i].first, 
                                         std::min(bounds[i].second, neighbor[i]));
                }
            }
            
            double neighbor_value = objective(neighbor);
            double delta = neighbor_value - current_value;
            
            // Accept or reject
            if (delta < 0 || uniform(gen) < std::exp(-delta / temperature)) {
                current = neighbor;
                current_value = neighbor_value;
                
                if (current_value < best_value) {
                    best = current;
                    best_value = current_value;
                }
            }
            
            // Cool down
            if (iter % iterations_per_temp == 0) {
                temperature *= cooling_rate;
            }
            
            // Check for convergence
            if (temperature < 1e-6) {
                result.converged = true;
                break;
            }
        }
        
        result.solution = best;
        result.objective_value = best_value;
        result.iterations = max_iterations;
        result.gradient_norm = -1; // Not applicable
        
        return result;
    }
    
    // ===== CONSTRAINT OPTIMIZATION - PENALTY METHOD =====
    
    class PenaltyMethod {
    private:
        std::function<double(const std::vector<double>&)> objective;
        std::vector<std::function<double(const std::vector<double>&)>> equality_constraints;
        std::vector<std::function<double(const std::vector<double>&)>> inequality_constraints;
        
        double penaltyFunction(const std::vector<double>& x, double penalty_param) {
            double f = objective(x);
            
            // Equality constraints: penalty = μ * Σ(g_i(x)²)
            for (const auto& constraint : equality_constraints) {
                double g = constraint(x);
                f += penalty_param * g * g;
            }
            
            // Inequality constraints: penalty = μ * Σ(max(0, h_i(x))²)
            for (const auto& constraint : inequality_constraints) {
                double h = constraint(x);
                if (h > 0) {
                    f += penalty_param * h * h;
                }
            }
            
            return f;
        }
        
    public:
        PenaltyMethod(const std::function<double(const std::vector<double>&)>& obj) 
            : objective(obj) {}
        
        void addEqualityConstraint(const std::function<double(const std::vector<double>&)>& constraint) {
            equality_constraints.push_back(constraint);
        }
        
        void addInequalityConstraint(const std::function<double(const std::vector<double>&)>& constraint) {
            inequality_constraints.push_back(constraint);
        }
        
        OptimizationResult solve(const std::vector<double>& initial_point,
                               double initial_penalty = 1.0,
                               double penalty_growth = 10.0,
                               int max_outer_iterations = 20) {
            std::vector<double> x = initial_point;
            double penalty_param = initial_penalty;
            
            OptimizationResult result;
            result.method = "Penalty Method";
            result.iterations = 0;
            
            for (int outer_iter = 0; outer_iter < max_outer_iterations; ++outer_iter) {
                // Create penalized objective
                auto penalized_obj = [this, penalty_param](const std::vector<double>& vars) {
                    return penaltyFunction(vars, penalty_param);
                };
                
                // Solve unconstrained problem
                auto inner_result = gradientDescent(penalized_obj, x, 0.001, 1e-6, 500);
                x = inner_result.solution;
                result.iterations += inner_result.iterations;
                
                // Check constraint satisfaction
                bool constraints_satisfied = true;
                for (const auto& constraint : equality_constraints) {
                    if (std::abs(constraint(x)) > 1e-6) {
                        constraints_satisfied = false;
                        break;
                    }
                }
                
                if (constraints_satisfied) {
                    for (const auto& constraint : inequality_constraints) {
                        if (constraint(x) > 1e-6) {
                            constraints_satisfied = false;
                            break;
                        }
                    }
                }
                
                if (constraints_satisfied) {
                    result.converged = true;
                    break;
                }
                
                penalty_param *= penalty_growth;
            }
            
            result.solution = x;
            result.objective_value = objective(x);
            result.gradient_norm = vectorNorm(computeGradient(objective, x));
            
            return result;
        }
    };
    
    // ===== GENETIC ALGORITHM =====
    
    OptimizationResult geneticAlgorithm(const std::function<double(const std::vector<double>&)>& objective,
                                       const std::vector<std::pair<double, double>>& bounds,
                                       int population_size = 50,
                                       int generations = 200,
                                       double mutation_rate = 0.1,
                                       double crossover_rate = 0.7) {
        std::random_device rd;
        std::mt19937 gen(rd());
        
        int dimensions = bounds.size();
        
        // Initialize population
        std::vector<std::vector<double>> population(population_size, std::vector<double>(dimensions));
        for (int i = 0; i < population_size; ++i) {
            for (int j = 0; j < dimensions; ++j) {
                std::uniform_real_distribution<double> dist(bounds[j].first, bounds[j].second);
                population[i][j] = dist(gen);
            }
        }
        
        OptimizationResult result;
        result.method = "Genetic Algorithm";
        result.iterations = generations;
        
        std::vector<double> best_individual;
        double best_fitness = INF;
        
        for (int gen = 0; gen < generations; ++gen) {
            // Evaluate fitness
            std::vector<double> fitness(population_size);
            for (int i = 0; i < population_size; ++i) {
                fitness[i] = objective(population[i]);
                if (fitness[i] < best_fitness) {
                    best_fitness = fitness[i];
                    best_individual = population[i];
                }
            }
            
            // Selection (tournament selection)
            std::vector<std::vector<double>> new_population;
            for (int i = 0; i < population_size; ++i) {
                std::uniform_int_distribution<int> idx_dist(0, population_size - 1);
                int idx1 = idx_dist(gen);
                int idx2 = idx_dist(gen);
                
                if (fitness[idx1] < fitness[idx2]) {
                    new_population.push_back(population[idx1]);
                } else {
                    new_population.push_back(population[idx2]);
                }
            }
            
            // Crossover and mutation
            for (int i = 0; i < population_size - 1; i += 2) {
                std::uniform_real_distribution<double> prob(0.0, 1.0);
                
                // Crossover
                if (prob(gen) < crossover_rate) {
                    std::uniform_int_distribution<int> point_dist(1, dimensions - 1);
                    int crossover_point = point_dist(gen);
                    
                    for (int j = crossover_point; j < dimensions; ++j) {
                        std::swap(new_population[i][j], new_population[i + 1][j]);
                    }
                }
                
                // Mutation
                for (int j = 0; j < dimensions; ++j) {
                    if (prob(gen) < mutation_rate) {
                        std::uniform_real_distribution<double> mut_dist(bounds[j].first, bounds[j].second);
                        new_population[i][j] = mut_dist(gen);
                    }
                    if (prob(gen) < mutation_rate) {
                        std::uniform_real_distribution<double> mut_dist(bounds[j].first, bounds[j].second);
                        new_population[i + 1][j] = mut_dist(gen);
                    }
                }
            }
            
            population = new_population;
        }
        
        result.solution = best_individual;
        result.objective_value = best_fitness;
        result.converged = true; // GA always "converges" after set generations
        result.gradient_norm = -1; // Not applicable
        
        return result;
    }
    
} // namespace OptimizationAlgorithms

// ===== DEMONSTRATION =====

int main() {
    using namespace OptimizationAlgorithms;
    
    std::cout << "=== OPTIMIZATION ALGORITHMS MODULE DEMO ===" << std::endl << std::endl;
    
    // ===== LINEAR PROGRAMMING =====
    std::cout << "1. LINEAR PROGRAMMING (Simplex Method)" << std::endl;
    std::cout << "Problem: minimize 3x + 2y subject to:" << std::endl;
    std::cout << "         x + y >= 1" << std::endl;
    std::cout << "         2x + y >= 2" << std::endl;
    std::cout << "         x, y >= 0" << std::endl << std::endl;
    
    // Convert to standard form: minimize c^T x subject to Ax <= b
    // We need to negate constraints: -x - y <= -1, -2x - y <= -2
    std::vector<double> c = {3.0, 2.0};
    std::vector<std::vector<double>> A = {{-1.0, -1.0}, {-2.0, -1.0}};
    std::vector<double> b = {-1.0, -2.0};
    
    SimplexSolver simplex;
    auto lp_result = simplex.solve(c, A, b);
    lp_result.print();
    
    // ===== NONLINEAR OPTIMIZATION =====
    std::cout << "2. NONLINEAR OPTIMIZATION" << std::endl;
    std::cout << "Problem: minimize f(x,y) = (x-3)² + (y-2)²" << std::endl;
    std::cout << "Optimal solution should be (3, 2) with f=0" << std::endl << std::endl;
    
    auto quadratic_objective = [](const std::vector<double>& x) {
        return (x[0] - 3.0) * (x[0] - 3.0) + (x[1] - 2.0) * (x[1] - 2.0);
    };
    
    std::vector<double> initial_point = {0.0, 0.0};
    
    // Gradient Descent
    auto gd_result = gradientDescent(quadratic_objective, initial_point, 0.1, 1e-8, 1000);
    gd_result.print();
    
    // Newton's Method
    auto newton_result = newtonMethod(quadratic_objective, initial_point, 1e-8, 100);
    newton_result.print();
    
    // ===== SIMULATED ANNEALING =====
    std::cout << "3. SIMULATED ANNEALING" << std::endl;
    std::cout << "Problem: minimize Rosenbrock function f(x,y) = 100(y-x²)² + (1-x)²" << std::endl;
    std::cout << "Optimal solution should be (1, 1) with f=0" << std::endl << std::endl;
    
    auto rosenbrock = [](const std::vector<double>& x) {
        double dx = x[1] - x[0] * x[0];
        double dy = 1.0 - x[0];
        return 100.0 * dx * dx + dy * dy;
    };
    
    std::vector<std::pair<double, double>> bounds = {{-5.0, 5.0}, {-5.0, 5.0}};
    auto sa_result = simulatedAnnealing(rosenbrock, {-1.0, -1.0}, bounds, 10.0, 0.95, 50, 5000);
    sa_result.print();
    
    // ===== CONSTRAINED OPTIMIZATION =====
    std::cout << "4. CONSTRAINED OPTIMIZATION (Penalty Method)" << std::endl;
    std::cout << "Problem: minimize f(x,y) = x² + y²" << std::endl;
    std::cout << "Subject to: x + y = 2 (equality constraint)" << std::endl;
    std::cout << "Optimal solution should be (1, 1) with f=2" << std::endl << std::endl;
    
    auto circle_objective = [](const std::vector<double>& x) {
        return x[0] * x[0] + x[1] * x[1];
    };
    
    PenaltyMethod penalty_solver(circle_objective);
    
    // Add equality constraint: g(x,y) = x + y - 2 = 0
    penalty_solver.addEqualityConstraint([](const std::vector<double>& x) {
        return x[0] + x[1] - 2.0;
    });
    
    auto penalty_result = penalty_solver.solve({0.0, 0.0}, 1.0, 10.0, 15);
    penalty_result.print();
    
    // ===== GENETIC ALGORITHM =====
    std::cout << "5. GENETIC ALGORITHM" << std::endl;
    std::cout << "Problem: minimize Rastrigin function (highly multimodal)" << std::endl;
    std::cout << "f(x,y) = 20 + x² + y² - 10cos(2πx) - 10cos(2πy)" << std::endl;
    std::cout << "Optimal solution should be (0, 0) with f=0" << std::endl << std::endl;
    
    auto rastrigin = [](const std::vector<double>& x) {
        const double PI = 3.14159265359;
        double sum = 20.0;
        for (double xi : x) {
            sum += xi * xi - 10.0 * std::cos(2.0 * PI * xi);
        }
        return sum;
    };
    
    std::vector<std::pair<double, double>> rastrigin_bounds = {{-5.12, 5.12}, {-5.12, 5.12}};
    auto ga_result = geneticAlgorithm(rastrigin, rastrigin_bounds, 100, 300, 0.1, 0.7);
    ga_result.print();
    
    // ===== COMPARISON SUMMARY =====
    std::cout << "6. ALGORITHM COMPARISON SUMMARY" << std::endl;
    std::cout << "Algorithm                 | Converged | Iterations | Final Objective" << std::endl;
    std::cout << "========================= | ========= | ========== | ===============" << std::endl;
    std::cout << "Gradient Descent         | " << (gd_result.converged ? "Yes" : "No") 
              << "       | " << std::setw(10) << gd_result.iterations 
              << " | " << std::scientific << gd_result.objective_value << std::endl;
    std::cout << "Newton's Method          | " << (newton_result.converged ? "Yes" : "No") 
              << "       | " << std::setw(10) << newton_result.iterations 
              << " | " << std::scientific << newton_result.objective_value << std::endl;
    std::cout << "Simulated Annealing      | " << (sa_result.converged ? "Yes" : "No") 
              << "       | " << std::setw(10) << sa_result.iterations 
              << " | " << std::scientific << sa_result.objective_value << std::endl;
    std::cout << "Penalty Method           | " << (penalty_result.converged ? "Yes" : "No") 
              << "       | " << std::setw(10) << penalty_result.iterations 
              << " | " << std::scientific << penalty_result.objective_value << std::endl;
    std::cout << "Genetic Algorithm        | " << (ga_result.converged ? "Yes" : "No") 
              << "       | " << std::setw(10) << ga_result.iterations 
              << " | " << std::scientific << ga_result.objective_value << std::endl;
    std::cout << std::fixed << std::endl;
    
    // ===== ALGORITHM SELECTION GUIDE =====
    std::cout << "7. ALGORITHM SELECTION GUIDE" << std::endl;
    std::cout << "• Linear Programming (Simplex): Use for linear objectives with linear constraints" << std::endl;
    std::cout << "• Gradient Descent: Fast, simple, good for smooth convex functions" << std::endl;
    std::cout << "• Newton's Method: Very fast convergence, requires Hessian, good for smooth functions" << std::endl;
    std::cout << "• Simulated Annealing: Global optimization, handles non-convex/discontinuous functions" << std::endl;
    std::cout << "• Penalty Method: Handles equality/inequality constraints systematically" << std::endl;
    std::cout << "• Genetic Algorithm: Global optimization, no gradient needed, handles discrete variables" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== OPTIMIZATION ALGORITHMS MODULE COMPLETE ===" << std::endl;
    
    return 0;
}












