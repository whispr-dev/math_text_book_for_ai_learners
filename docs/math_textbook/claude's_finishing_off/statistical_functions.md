#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <map>

namespace StatisticalFunctions {

    // ===== REGRESSION ANALYSIS =====
    
    struct LinearRegressionResult {
        double slope;
        double intercept;
        double correlation_coefficient;
        double r_squared;
    };
    
    // Linear Regression: y = mx + b
    LinearRegressionResult linearRegression(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size() || x.empty()) {
            throw std::invalid_argument("X and Y vectors must be the same size and non-empty");
        }
        
        size_t n = x.size();
        double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
        double sum_y = std::accumulate(y.begin(), y.end(), 0.0);
        double sum_xy = 0.0, sum_x2 = 0.0, sum_y2 = 0.0;
        
        for (size_t i = 0; i < n; ++i) {
            sum_xy += x[i] * y[i];
            sum_x2 += x[i] * x[i];
            sum_y2 += y[i] * y[i];
        }
        
        double mean_x = sum_x / n;
        double mean_y = sum_y / n;
        
        // Calculate slope and intercept
        double numerator = n * sum_xy - sum_x * sum_y;
        double denominator = n * sum_x2 - sum_x * sum_x;
        
        if (std::abs(denominator) < 1e-10) {
            throw std::runtime_error("Cannot compute regression: vertical line");
        }
        
        double slope = numerator / denominator;
        double intercept = mean_y - slope * mean_x;
        
        // Calculate correlation coefficient
        double r_num = numerator;
        double r_den = std::sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y));
        double correlation = r_den > 1e-10 ? r_num / r_den : 0.0;
        
        return {slope, intercept, correlation, correlation * correlation};
    }
    
    // Polynomial Regression using least squares
    std::vector<double> polynomialRegression(const std::vector<double>& x, const std::vector<double>& y, int degree) {
        if (x.size() != y.size() || x.empty() || degree < 0) {
            throw std::invalid_argument("Invalid input for polynomial regression");
        }
        
        size_t n = x.size();
        int m = degree + 1;
        
        // Create Vandermonde matrix
        std::vector<std::vector<double>> A(n, std::vector<double>(m));
        for (size_t i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                A[i][j] = std::pow(x[i], j);
            }
        }
        
        // Solve normal equations: A^T * A * coeff = A^T * y
        std::vector<std::vector<double>> AtA(m, std::vector<double>(m, 0.0));
        std::vector<double> Aty(m, 0.0);
        
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                for (size_t k = 0; k < n; ++k) {
                    AtA[i][j] += A[k][i] * A[k][j];
                }
            }
            for (size_t k = 0; k < n; ++k) {
                Aty[i] += A[k][i] * y[k];
            }
        }
        
        // Gaussian elimination to solve AtA * coeff = Aty
        std::vector<double> coefficients(m);
        for (int i = 0; i < m; ++i) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k < m; ++k) {
                if (std::abs(AtA[k][i]) > std::abs(AtA[maxRow][i])) {
                    maxRow = k;
                }
            }
            std::swap(AtA[i], AtA[maxRow]);
            std::swap(Aty[i], Aty[maxRow]);
            
            // Forward elimination
            for (int k = i + 1; k < m; ++k) {
                double factor = AtA[k][i] / AtA[i][i];
                for (int j = i; j < m; ++j) {
                    AtA[k][j] -= factor * AtA[i][j];
                }
                Aty[k] -= factor * Aty[i];
            }
        }
        
        // Back substitution
        for (int i = m - 1; i >= 0; --i) {
            coefficients[i] = Aty[i];
            for (int j = i + 1; j < m; ++j) {
                coefficients[i] -= AtA[i][j] * coefficients[j];
            }
            coefficients[i] /= AtA[i][i];
        }
        
        return coefficients;
    }
    
    // ===== CORRELATION AND COVARIANCE =====
    
    // Pearson correlation coefficient
    double correlation(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size() || x.empty()) {
            throw std::invalid_argument("Vectors must be same size and non-empty");
        }
        
        size_t n = x.size();
        double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
        double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
        
        double numerator = 0.0, sum_sq_x = 0.0, sum_sq_y = 0.0;
        
        for (size_t i = 0; i < n; ++i) {
            double dx = x[i] - mean_x;
            double dy = y[i] - mean_y;
            numerator += dx * dy;
            sum_sq_x += dx * dx;
            sum_sq_y += dy * dy;
        }
        
        double denominator = std::sqrt(sum_sq_x * sum_sq_y);
        return denominator > 1e-10 ? numerator / denominator : 0.0;
    }
    
    // Covariance
    double covariance(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size() || x.empty()) {
            throw std::invalid_argument("Vectors must be same size and non-empty");
        }
        
        size_t n = x.size();
        double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
        double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
        
        double sum = 0.0;
        for (size_t i = 0; i < n; ++i) {
            sum += (x[i] - mean_x) * (y[i] - mean_y);
        }
        
        return sum / (n - 1); // Sample covariance
    }
    
    // Covariance matrix for multiple variables
    std::vector<std::vector<double>> covarianceMatrix(const std::vector<std::vector<double>>& data) {
        if (data.empty() || data[0].empty()) {
            throw std::invalid_argument("Data matrix cannot be empty");
        }
        
        size_t vars = data.size();
        size_t samples = data[0].size();
        
        // Check all variables have same number of samples
        for (const auto& var : data) {
            if (var.size() != samples) {
                throw std::invalid_argument("All variables must have same number of samples");
            }
        }
        
        std::vector<std::vector<double>> cov_matrix(vars, std::vector<double>(vars));
        
        for (size_t i = 0; i < vars; ++i) {
            for (size_t j = 0; j < vars; ++j) {
                cov_matrix[i][j] = covariance(data[i], data[j]);
            }
        }
        
        return cov_matrix;
    }
    
    // ===== RANDOM NUMBER GENERATION =====
    
    class RandomGenerator {
    private:
        std::mt19937 gen;
        
    public:
        RandomGenerator() : gen(std::random_device{}()) {}
        explicit RandomGenerator(unsigned seed) : gen(seed) {}
        
        // Uniform distribution [min, max]
        std::vector<double> uniform(size_t count, double min = 0.0, double max = 1.0) {
            std::uniform_real_distribution<double> dist(min, max);
            std::vector<double> result(count);
            for (size_t i = 0; i < count; ++i) {
                result[i] = dist(gen);
            }
            return result;
        }
        
        // Normal (Gaussian) distribution
        std::vector<double> normal(size_t count, double mean = 0.0, double stddev = 1.0) {
            std::normal_distribution<double> dist(mean, stddev);
            std::vector<double> result(count);
            for (size_t i = 0; i < count; ++i) {
                result[i] = dist(gen);
            }
            return result;
        }
        
        // Binomial distribution
        std::vector<int> binomial(size_t count, int trials, double probability) {
            std::binomial_distribution<int> dist(trials, probability);
            std::vector<int> result(count);
            for (size_t i = 0; i < count; ++i) {
                result[i] = dist(gen);
            }
            return result;
        }
        
        // Exponential distribution
        std::vector<double> exponential(size_t count, double lambda = 1.0) {
            std::exponential_distribution<double> dist(lambda);
            std::vector<double> result(count);
            for (size_t i = 0; i < count; ++i) {
                result[i] = dist(gen);
            }
            return result;
        }
        
        // Poisson distribution
        std::vector<int> poisson(size_t count, double mean) {
            std::poisson_distribution<int> dist(mean);
            std::vector<int> result(count);
            for (size_t i = 0; i < count; ++i) {
                result[i] = dist(gen);
            }
            return result;
        }
    };
    
    // ===== HYPOTHESIS TESTING =====
    
    struct TTestResult {
        double t_statistic;
        double degrees_of_freedom;
        double p_value_approx; // Simplified p-value approximation
    };
    
    // One-sample t-test
    TTestResult oneSampleTTest(const std::vector<double>& sample, double population_mean) {
        if (sample.empty()) {
            throw std::invalid_argument("Sample cannot be empty");
        }
        
        size_t n = sample.size();
        double sample_mean = std::accumulate(sample.begin(), sample.end(), 0.0) / n;
        
        // Calculate sample standard deviation
        double sum_sq_diff = 0.0;
        for (double x : sample) {
            double diff = x - sample_mean;
            sum_sq_diff += diff * diff;
        }
        double sample_std = std::sqrt(sum_sq_diff / (n - 1));
        
        // Calculate t-statistic
        double t_stat = (sample_mean - population_mean) / (sample_std / std::sqrt(n));
        double df = n - 1;
        
        // Simplified p-value approximation (for demonstration)
        double p_value = 2.0 * (1.0 - std::abs(t_stat) / (std::abs(t_stat) + df));
        
        return {t_stat, df, p_value};
    }
    
    // Two-sample t-test (assuming equal variances)
    TTestResult twoSampleTTest(const std::vector<double>& sample1, const std::vector<double>& sample2) {
        if (sample1.empty() || sample2.empty()) {
            throw std::invalid_argument("Samples cannot be empty");
        }
        
        size_t n1 = sample1.size(), n2 = sample2.size();
        
        double mean1 = std::accumulate(sample1.begin(), sample1.end(), 0.0) / n1;
        double mean2 = std::accumulate(sample2.begin(), sample2.end(), 0.0) / n2;
        
        // Calculate sample variances
        double var1 = 0.0, var2 = 0.0;
        for (double x : sample1) {
            double diff = x - mean1;
            var1 += diff * diff;
        }
        var1 /= (n1 - 1);
        
        for (double x : sample2) {
            double diff = x - mean2;
            var2 += diff * diff;
        }
        var2 /= (n2 - 1);
        
        // Pooled standard error
        double pooled_se = std::sqrt(var1/n1 + var2/n2);
        double t_stat = (mean1 - mean2) / pooled_se;
        double df = n1 + n2 - 2;
        
        // Simplified p-value approximation
        double p_value = 2.0 * (1.0 - std::abs(t_stat) / (std::abs(t_stat) + df));
        
        return {t_stat, df, p_value};
    }
    
    struct ChiSquareResult {
        double chi_square_statistic;
        double degrees_of_freedom;
        double p_value_approx;
    };
    
    // Chi-square goodness of fit test
    ChiSquareResult chiSquareTest(const std::vector<int>& observed, const std::vector<double>& expected) {
        if (observed.size() != expected.size() || observed.empty()) {
            throw std::invalid_argument("Observed and expected must be same size and non-empty");
        }
        
        double chi_square = 0.0;
        for (size_t i = 0; i < observed.size(); ++i) {
            if (expected[i] <= 0) {
                throw std::invalid_argument("Expected frequencies must be positive");
            }
            double diff = observed[i] - expected[i];
            chi_square += (diff * diff) / expected[i];
        }
        
        double df = observed.size() - 1;
        
        // Very simplified p-value approximation
        double p_value = std::exp(-chi_square / 2.0);
        
        return {chi_square, df, p_value};
    }
    
} // namespace StatisticalFunctions

// ===== DEMONSTRATION AND USAGE =====

int main() {
    using namespace StatisticalFunctions;
    
    std::cout << "=== STATISTICAL FUNCTIONS MODULE DEMO ===" << std::endl << std::endl;
    
    // Sample data
    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> y = {2.1, 3.9, 6.2, 7.8, 10.1, 12.2, 13.9, 16.1, 18.0, 20.2};
    
    // === LINEAR REGRESSION ===
    std::cout << "1. LINEAR REGRESSION" << std::endl;
    try {
        auto lr_result = linearRegression(x, y);
        std::cout << "Slope: " << lr_result.slope << std::endl;
        std::cout << "Intercept: " << lr_result.intercept << std::endl;
        std::cout << "Correlation: " << lr_result.correlation_coefficient << std::endl;
        std::cout << "R-squared: " << lr_result.r_squared << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    std::cout << std::endl;
    
    // === POLYNOMIAL REGRESSION ===
    std::cout << "2. POLYNOMIAL REGRESSION (degree 2)" << std::endl;
    try {
        auto poly_coeff = polynomialRegression(x, y, 2);
        std::cout << "Coefficients: ";
        for (size_t i = 0; i < poly_coeff.size(); ++i) {
            std::cout << "a" << i << "=" << poly_coeff[i] << " ";
        }
        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    std::cout << std::endl;
    
    // === CORRELATION AND COVARIANCE ===
    std::cout << "3. CORRELATION AND COVARIANCE" << std::endl;
    std::cout << "Correlation: " << correlation(x, y) << std::endl;
    std::cout << "Covariance: " << covariance(x, y) << std::endl;
    std::cout << std::endl;
    
    // === RANDOM NUMBER GENERATION ===
    std::cout << "4. RANDOM NUMBER GENERATION" << std::endl;
    RandomGenerator rng(42); // Fixed seed for reproducibility
    
    auto uniform_nums = rng.uniform(5, 0.0, 10.0);
    std::cout << "Uniform [0,10]: ";
    for (double val : uniform_nums) std::cout << val << " ";
    std::cout << std::endl;
    
    auto normal_nums = rng.normal(5, 0.0, 1.0);
    std::cout << "Normal (0,1): ";
    for (double val : normal_nums) std::cout << val << " ";
    std::cout << std::endl;
    
    auto binomial_nums = rng.binomial(5, 10, 0.5);
    std::cout << "Binomial (10,0.5): ";
    for (int val : binomial_nums) std::cout << val << " ";
    std::cout << std::endl;
    std::cout << std::endl;
    
    // === HYPOTHESIS TESTING ===
    std::cout << "5. HYPOTHESIS TESTING" << std::endl;
    
    // One-sample t-test
    std::vector<double> sample = {2.3, 2.1, 2.5, 2.0, 2.4, 2.2, 2.6, 2.1, 2.3, 2.4};
    auto t_result = oneSampleTTest(sample, 2.0);
    std::cout << "One-sample t-test (H0: μ = 2.0):" << std::endl;
    std::cout << "  t-statistic: " << t_result.t_statistic << std::endl;
    std::cout << "  df: " << t_result.degrees_of_freedom << std::endl;
    std::cout << "  p-value (approx): " << t_result.p_value_approx << std::endl;
    
    // Chi-square test
    std::vector<int> observed = {20, 25, 15, 10};
    std::vector<double> expected = {17.5, 17.5, 17.5, 17.5};
    auto chi_result = chiSquareTest(observed, expected);
    std::cout << "Chi-square goodness of fit:" << std::endl;
    std::cout << "  χ² statistic: " << chi_result.chi_square_statistic << std::endl;
    std::cout << "  df: " << chi_result.degrees_of_freedom << std::endl;
    std::cout << "  p-value (approx): " << chi_result.p_value_approx << std::endl;
    
    return 0;
}