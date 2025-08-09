#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>

namespace StatisticalFunctions {
    
    // Helper function to calculate mean
    double mean(const std::vector<double>& data) {
        return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    }
    
    // Helper function to calculate variance
    double variance(const std::vector<double>& data) {
        double m = mean(data);
        double sum = 0.0;
        for (double x : data) {
            sum += (x - m) * (x - m);
        }
        return sum / (data.size() - 1);
    }
    
    // Helper function to calculate standard deviation
    double standardDeviation(const std::vector<double>& data) {
        return std::sqrt(variance(data));
    }
    
    // Linear Regression - returns pair of (slope, intercept)
    std::pair<double, double> linearRegression(const std::vector<double>& x, const std::vector<double>& y) {
        size_t n = x.size();
        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
        
        for (size_t i = 0; i < n; i++) {
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += x[i] * y[i];
            sum_x2 += x[i] * x[i];
        }
        
        double slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
        double intercept = (sum_y - slope * sum_x) / n;
        
        return {slope, intercept};
    }
    
    // Polynomial Regression - returns coefficients
    std::vector<double> polynomialRegression(const std::vector<double>& x, const std::vector<double>& y, int degree) {
        int n = x.size();
        int np = degree + 1;  // number of parameters
        
        // Build the Vandermonde matrix
        std::vector<std::vector<double>> A(n, std::vector<double>(np));
        for (int i = 0; i < n; i++) {
            double xi = x[i];
            A[i][0] = 1.0;
            for (int j = 1; j < np; j++) {
                A[i][j] = A[i][j-1] * xi;
            }
        }
        
        // Solve using normal equations (A^T * A) * coeffs = A^T * y
        std::vector<std::vector<double>> ATA(np, std::vector<double>(np, 0));
        std::vector<double> ATy(np, 0);
        
        for (int i = 0; i < np; i++) {
            for (int j = 0; j < np; j++) {
                for (int k = 0; k < n; k++) {
                    ATA[i][j] += A[k][i] * A[k][j];
                }
            }
            for (int k = 0; k < n; k++) {
                ATy[i] += A[k][i] * y[k];
            }
        }
        
        // Gaussian elimination to solve the system
        for (int i = 0; i < np; i++) {
            // Partial pivoting
            int maxRow = i;
            for (int k = i + 1; k < np; k++) {
                if (std::abs(ATA[k][i]) > std::abs(ATA[maxRow][i])) {
                    maxRow = k;
                }
            }
            std::swap(ATA[maxRow], ATA[i]);
            std::swap(ATy[maxRow], ATy[i]);
            
            // Forward elimination
            for (int k = i + 1; k < np; k++) {
                double factor = ATA[k][i] / ATA[i][i];
                for (int j = i; j < np; j++) {
                    ATA[k][j] -= factor * ATA[i][j];
                }
                ATy[k] -= factor * ATy[i];
            }
        }
        
        // Back substitution
        std::vector<double> coeffs(np);
        for (int i = np - 1; i >= 0; i--) {
            coeffs[i] = ATy[i];
            for (int j = i + 1; j < np; j++) {
                coeffs[i] -= ATA[i][j] * coeffs[j];
            }
            coeffs[i] /= ATA[i][i];
        }
        
        return coeffs;
    }
    
    // Correlation coefficient (Pearson)
    double correlation(const std::vector<double>& x, const std::vector<double>& y) {
        size_t n = x.size();
        double mean_x = mean(x);
        double mean_y = mean(y);
        
        double sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
        for (size_t i = 0; i < n; i++) {
            double dx = x[i] - mean_x;
            double dy = y[i] - mean_y;
            sum_xy += dx * dy;
            sum_x2 += dx * dx;
            sum_y2 += dy * dy;
        }
        
        return sum_xy / (std::sqrt(sum_x2) * std::sqrt(sum_y2));
    }
    
    // Covariance matrix
    std::vector<std::vector<double>> covarianceMatrix(const std::vector<std::vector<double>>& data) {
        size_t n = data.size();     // number of observations
        size_t m = data[0].size();  // number of variables
        
        // Calculate means for each variable
        std::vector<double> means(m, 0);
        for (size_t j = 0; j < m; j++) {
            for (size_t i = 0; i < n; i++) {
                means[j] += data[i][j];
            }
            means[j] /= n;
        }
        
        // Calculate covariance matrix
        std::vector<std::vector<double>> cov(m, std::vector<double>(m, 0));
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < m; j++) {
                double sum = 0;
                for (size_t k = 0; k < n; k++) {
                    sum += (data[k][i] - means[i]) * (data[k][j] - means[j]);
                }
                cov[i][j] = sum / (n - 1);
            }
        }
        
        return cov;
    }
    
    // T-test (two-sample, returns t-statistic and p-value approximation)
    std::pair<double, double> tTest(const std::vector<double>& sample1, const std::vector<double>& sample2) {
        double mean1 = mean(sample1);
        double mean2 = mean(sample2);
        double var1 = variance(sample1);
        double var2 = variance(sample2);
        size_t n1 = sample1.size();
        size_t n2 = sample2.size();
        
        // Calculate t-statistic
        double pooledStdError = std::sqrt(var1/n1 + var2/n2);
        double tStat = (mean1 - mean2) / pooledStdError;
        
        // Degrees of freedom (Welch's approximation)
        double df = std::pow(var1/n1 + var2/n2, 2) / 
                   (std::pow(var1/n1, 2)/(n1-1) + std::pow(var2/n2, 2)/(n2-1));
        
        // Approximate p-value using normal distribution (simplified)
        // For production, use a proper t-distribution CDF
        double pValue = 2 * (1 - 0.5 * std::erfc(-std::abs(tStat) / std::sqrt(2)));
        
        return {tStat, pValue};
    }
    
    // Chi-square test for independence
    double chiSquareTest(const std::vector<std::vector<int>>& observed) {
        size_t rows = observed.size();
        size_t cols = observed[0].size();
        
        // Calculate row and column totals
        std::vector<int> rowTotals(rows, 0);
        std::vector<int> colTotals(cols, 0);
        int total = 0;
        
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                rowTotals[i] += observed[i][j];
                colTotals[j] += observed[i][j];
                total += observed[i][j];
            }
        }
        
        // Calculate chi-square statistic
        double chiSquare = 0;
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                double expected = static_cast<double>(rowTotals[i] * colTotals[j]) / total;
                double diff = observed[i][j] - expected;
                chiSquare += (diff * diff) / expected;
            }
        }
        
        return chiSquare;
    }
    
    // Random Number Generators for different distributions
    class RandomGenerator {
    private:
        std::mt19937 gen;
        
    public:
        RandomGenerator() : gen(std::random_device{}()) {}
        
        // Uniform distribution
        std::vector<double> uniform(size_t n, double min, double max) {
            std::uniform_real_distribution<> dist(min, max);
            std::vector<double> result(n);
            for (size_t i = 0; i < n; i++) {
                result[i] = dist(gen);
            }
            return result;
        }
        
        // Normal distribution
        std::vector<double> normal(size_t n, double mean, double stddev) {
            std::normal_distribution<> dist(mean, stddev);
            std::vector<double> result(n);
            for (size_t i = 0; i < n; i++) {
                result[i] = dist(gen);
            }
            return result;
        }
        
        // Binomial distribution
        std::vector<int> binomial(size_t n, int trials, double p) {
            std::binomial_distribution<> dist(trials, p);
            std::vector<int> result(n);
            for (size_t i = 0; i < n; i++) {
                result[i] = dist(gen);
            }
            return result;
        }
        
        // Poisson distribution
        std::vector<int> poisson(size_t n, double lambda) {
            std::poisson_distribution<> dist(lambda);
            std::vector<int> result(n);
            for (size_t i = 0; i < n; i++) {
                result[i] = dist(gen);
            }
            return result;
        }
        
        // Exponential distribution
        std::vector<double> exponential(size_t n, double lambda) {
            std::exponential_distribution<> dist(lambda);
            std::vector<double> result(n);
            for (size_t i = 0; i < n; i++) {
                result[i] = dist(gen);
            }
            return result;
        }
    };
    
    // General function to demonstrate statistical functions
    void calc() {
        // Example data
        std::vector<double> x = {1, 2, 3, 4, 5};
        std::vector<double> y = {2.1, 3.9, 6.1, 8.0, 10.2};
        
        // Linear regression
        auto [slope, intercept] = linearRegression(x, y);
        std::cout << "Linear Regression: y = " << slope << "x + " << intercept << std::endl;
        
        // Polynomial regression (degree 2)
        auto polyCoeffs = polynomialRegression(x, y, 2);
        std::cout << "Polynomial Regression (degree 2): y = ";
        for (int i = polyCoeffs.size() - 1; i >= 0; i--) {
            std::cout << polyCoeffs[i];
            if (i > 0) std::cout << "x^" << i << " + ";
        }
        std::cout << std::endl;
        
        // Correlation
        double corr = correlation(x, y);
        std::cout << "Correlation coefficient: " << corr << std::endl;
        
        // Covariance matrix
        std::vector<std::vector<double>> data = {{1, 2}, {2, 3}, {3, 5}, {4, 7}, {5, 8}};
        auto cov = covarianceMatrix(data);
        std::cout << "Covariance Matrix:" << std::endl;
        for (const auto& row : cov) {
            for (double val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
        
        // Random number generation
        RandomGenerator rng;
        auto normalSample = rng.normal(5, 0, 1);
        std::cout << "Normal distribution sample (mean=0, stddev=1): ";
        for (double val : normalSample) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        
        // T-test
        std::vector<double> sample1 = {1.2, 2.3, 3.1, 4.0, 5.2};
        std::vector<double> sample2 = {2.1, 3.2, 4.0, 5.1, 6.0};
        auto [tStat, pValue] = tTest(sample1, sample2);
        std::cout << "T-test: t-statistic = " << tStat << ", p-value ≈ " << pValue << std::endl;
    }
}

int main() {
    StatisticalFunctions::calc();
    return 0;
}