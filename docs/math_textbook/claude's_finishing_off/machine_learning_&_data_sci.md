#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <limits>
#include <map>
#include <memory>

namespace MachineLearningUtilities {
    
    // Helper function to compute mean of a vector
    double mean(const std::vector<double>& v) {
        double sum = 0;
        for (double x : v) sum += x;
        return sum / v.size();
    }
    
    // Helper function for matrix multiplication
    std::vector<std::vector<double>> matrixMultiply(
        const std::vector<std::vector<double>>& A,
        const std::vector<std::vector<double>>& B) {
        
        size_t m = A.size();
        size_t n = B[0].size();
        size_t p = B.size();
        
        std::vector<std::vector<double>> C(m, std::vector<double>(n, 0));
        
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                for (size_t k = 0; k < p; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }
    
    // Helper function for matrix transpose
    std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& A) {
        size_t m = A.size();
        size_t n = A[0].size();
        std::vector<std::vector<double>> AT(n, std::vector<double>(m));
        
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                AT[j][i] = A[i][j];
            }
        }
        return AT;
    }
    
    // Principal Component Analysis (PCA)
    class PCA {
    private:
        std::vector<std::vector<double>> components;
        std::vector<double> eigenvalues;
        std::vector<double> means;
        int n_components;
        
        // Power iteration method to find dominant eigenvector
        std::pair<std::vector<double>, double> powerIteration(
            const std::vector<std::vector<double>>& matrix, int maxIter = 100) {
            
            size_t n = matrix.size();
            std::vector<double> v(n);
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0, 1);
            
            // Initialize with random vector
            for (size_t i = 0; i < n; i++) {
                v[i] = dis(gen);
            }
            
            double eigenvalue = 0;
            for (int iter = 0; iter < maxIter; iter++) {
                // Multiply matrix by vector
                std::vector<double> Av(n, 0);
                for (size_t i = 0; i < n; i++) {
                    for (size_t j = 0; j < n; j++) {
                        Av[i] += matrix[i][j] * v[j];
                    }
                }
                
                // Calculate eigenvalue (Rayleigh quotient)
                double vAv = 0, vv = 0;
                for (size_t i = 0; i < n; i++) {
                    vAv += v[i] * Av[i];
                    vv += v[i] * v[i];
                }
                eigenvalue = vAv / vv;
                
                // Normalize
                double norm = 0;
                for (double x : Av) norm += x * x;
                norm = std::sqrt(norm);
                
                for (size_t i = 0; i < n; i++) {
                    v[i] = Av[i] / norm;
                }
            }
            
            return {v, eigenvalue};
        }
        
        // Compute covariance matrix
        std::vector<std::vector<double>> computeCovariance(
            const std::vector<std::vector<double>>& data) {
            
            size_t n = data.size();
            size_t m = data[0].size();
            
            // Compute means
            means.resize(m);
            for (size_t j = 0; j < m; j++) {
                double sum = 0;
                for (size_t i = 0; i < n; i++) {
                    sum += data[i][j];
                }
                means[j] = sum / n;
            }
            
            // Compute covariance
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
        
    public:
        PCA(int n_comp = 2) : n_components(n_comp) {}
        
        void fit(const std::vector<std::vector<double>>& data) {
            auto cov = computeCovariance(data);
            components.clear();
            eigenvalues.clear();
            
            // Extract principal components using power iteration
            for (int comp = 0; comp < n_components; comp++) {
                auto [eigenvector, eigenvalue] = powerIteration(cov);
                components.push_back(eigenvector);
                eigenvalues.push_back(eigenvalue);
                
                // Deflate matrix for next component
                for (size_t i = 0; i < cov.size(); i++) {
                    for (size_t j = 0; j < cov.size(); j++) {
                        cov[i][j] -= eigenvalue * eigenvector[i] * eigenvector[j];
                    }
                }
            }
        }
        
        std::vector<std::vector<double>> transform(const std::vector<std::vector<double>>& data) {
            size_t n = data.size();
            std::vector<std::vector<double>> transformed(n, std::vector<double>(n_components));
            
            for (size_t i = 0; i < n; i++) {
                for (int j = 0; j < n_components; j++) {
                    double sum = 0;
                    for (size_t k = 0; k < data[i].size(); k++) {
                        sum += (data[i][k] - means[k]) * components[j][k];
                    }
                    transformed[i][j] = sum;
                }
            }
            
            return transformed;
        }
        
        double explainedVarianceRatio(int component) {
            double totalVar = 0;
            for (double ev : eigenvalues) totalVar += ev;
            return eigenvalues[component] / totalVar;
        }
    };
    
    // K-Means Clustering
    class KMeans {
    private:
        int k;
        int maxIters;
        std::vector<std::vector<double>> centroids;
        std::vector<int> labels;
        
        double euclideanDistance(const std::vector<double>& a, const std::vector<double>& b) {
            double sum = 0;
            for (size_t i = 0; i < a.size(); i++) {
                double diff = a[i] - b[i];
                sum += diff * diff;
            }
            return std::sqrt(sum);
        }
        
    public:
        KMeans(int num_clusters, int max_iterations = 100) 
            : k(num_clusters), maxIters(max_iterations) {}
        
        void fit(const std::vector<std::vector<double>>& data) {
            size_t n = data.size();
            size_t m = data[0].size();
            labels.resize(n);
            
            // Initialize centroids randomly from data points
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(0, n - 1);
            
            centroids.clear();
            std::set<int> chosen;
            while (centroids.size() < k) {
                int idx = dis(gen);
                if (chosen.find(idx) == chosen.end()) {
                    centroids.push_back(data[idx]);
                    chosen.insert(idx);
                }
            }
            
            // K-means iterations
            for (int iter = 0; iter < maxIters; iter++) {
                // Assignment step
                bool changed = false;
                for (size_t i = 0; i < n; i++) {
                    double minDist = std::numeric_limits<double>::max();
                    int bestCluster = 0;
                    
                    for (int j = 0; j < k; j++) {
                        double dist = euclideanDistance(data[i], centroids[j]);
                        if (dist < minDist) {
                            minDist = dist;
                            bestCluster = j;
                        }
                    }
                    
                    if (labels[i] != bestCluster) {
                        changed = true;
                        labels[i] = bestCluster;
                    }
                }
                
                if (!changed) break;
                
                // Update step
                for (int j = 0; j < k; j++) {
                    std::vector<double> newCentroid(m, 0);
                    int count = 0;
                    
                    for (size_t i = 0; i < n; i++) {
                        if (labels[i] == j) {
                            for (size_t d = 0; d < m; d++) {
                                newCentroid[d] += data[i][d];
                            }
                            count++;
                        }
                    }
                    
                    if (count > 0) {
                        for (size_t d = 0; d < m; d++) {
                            centroids[j][d] = newCentroid[d] / count;
                        }
                    }
                }
            }
        }
        
        std::vector<int> predict(const std::vector<std::vector<double>>& data) {
            std::vector<int> predictions(data.size());
            
            for (size_t i = 0; i < data.size(); i++) {
                double minDist = std::numeric_limits<double>::max();
                int bestCluster = 0;
                
                for (int j = 0; j < k; j++) {
                    double dist = euclideanDistance(data[i], centroids[j]);
                    if (dist < minDist) {
                        minDist = dist;
                        bestCluster = j;
                    }
                }
                predictions[i] = bestCluster;
            }
            
            return predictions;
        }
        
        double inertia() {
            double totalInertia = 0;
            for (size_t i = 0; i < labels.size(); i++) {
                double dist = euclideanDistance(centroids[labels[i]], centroids[labels[i]]);
                totalInertia += dist * dist;
            }
            return totalInertia;
        }
        
        std::vector<std::vector<double>> getCentroids() { return centroids; }
        std::vector<int> getLabels() { return labels; }
    };
    
    // Decision Tree
    class DecisionTree {
    private:
        struct Node {
            bool isLeaf;
            int featureIndex;
            double threshold;
            double prediction;
            std::unique_ptr<Node> left;
            std::unique_ptr<Node> right;
            
            Node() : isLeaf(false), featureIndex(-1), threshold(0), prediction(0) {}
        };
        
        std::unique_ptr<Node> root;
        int maxDepth;
        int minSamplesSplit;
        
        double giniImpurity(const std::vector<int>& labels) {
            std::map<int, int> counts;
            for (int label : labels) {
                counts[label]++;
            }
            
            double impurity = 1.0;
            int total = labels.size();
            for (const auto& [label, count] : counts) {
                double prob = static_cast<double>(count) / total;
                impurity -= prob * prob;
            }
            return impurity;
        }
        
        std::pair<int, double> findBestSplit(
            const std::vector<std::vector<double>>& X,
            const std::vector<int>& y,
            const std::vector<int>& indices) {
            
            int bestFeature = -1;
            double bestThreshold = 0;
            double bestGain = 0;
            
            int n_features = X[0].size();
            int n_samples = indices.size();
            
            // Calculate parent impurity
            std::vector<int> parentLabels;
            for (int idx : indices) {
                parentLabels.push_back(y[idx]);
            }
            double parentImpurity = giniImpurity(parentLabels);
            
            // Try each feature
            for (int feat = 0; feat < n_features; feat++) {
                // Get unique values for this feature
                std::set<double> uniqueVals;
                for (int idx : indices) {
                    uniqueVals.insert(X[idx][feat]);
                }
                
                // Try each unique value as threshold
                for (double val : uniqueVals) {
                    std::vector<int> leftLabels, rightLabels;
                    
                    for (int idx : indices) {
                        if (X[idx][feat] <= val) {
                            leftLabels.push_back(y[idx]);
                        } else {
                            rightLabels.push_back(y[idx]);
                        }
                    }
                    
                    if (leftLabels.empty() || rightLabels.empty()) continue;
                    
                    // Calculate information gain
                    double leftWeight = static_cast<double>(leftLabels.size()) / n_samples;
                    double rightWeight = static_cast<double>(rightLabels.size()) / n_samples;
                    double gain = parentImpurity - 
                                 (leftWeight * giniImpurity(leftLabels) + 
                                  rightWeight * giniImpurity(rightLabels));
                    
                    if (gain > bestGain) {
                        bestGain = gain;
                        bestFeature = feat;
                        bestThreshold = val;
                    }
                }
            }
            
            return {bestFeature, bestThreshold};
        }
        
        std::unique_ptr<Node> buildTree(
            const std::vector<std::vector<double>>& X,
            const std::vector<int>& y,
            const std::vector<int>& indices,
            int depth) {
            
            auto node = std::make_unique<Node>();
            
            // Check stopping criteria
            if (depth >= maxDepth || indices.size() < minSamplesSplit) {
                node->isLeaf = true;
                // Calculate majority class
                std::map<int, int> counts;
                for (int idx : indices) {
                    counts[y[idx]]++;
                }
                int maxCount = 0;
                for (const auto& [label, count] : counts) {
                    if (count > maxCount) {
                        maxCount = count;
                        node->prediction = label;
                    }
                }
                return node;
            }
            
            // Find best split
            auto [bestFeature, bestThreshold] = findBestSplit(X, y, indices);
            
            if (bestFeature == -1) {
                // No valid split found
                node->isLeaf = true;
                std::map<int, int> counts;
                for (int idx : indices) {
                    counts[y[idx]]++;
                }
                int maxCount = 0;
                for (const auto& [label, count] : counts) {
                    if (count > maxCount) {
                        maxCount = count;
                        node->prediction = label;
                    }
                }
                return node;
            }
            
            // Split data
            std::vector<int> leftIndices, rightIndices;
            for (int idx : indices) {
                if (X[idx][bestFeature] <= bestThreshold) {
                    leftIndices.push_back(idx);
                } else {
                    rightIndices.push_back(idx);
                }
            }
            
            // Build subtrees
            node->featureIndex = bestFeature;
            node->threshold = bestThreshold;
            node->left = buildTree(X, y, leftIndices, depth + 1);
            node->right = buildTree(X, y, rightIndices, depth + 1);
            
            return node;
        }
        
        int predictSample(const std::vector<double>& x, Node* node) {
            if (node->isLeaf) {
                return static_cast<int>(node->prediction);
            }
            
            if (x[node->featureIndex] <= node->threshold) {
                return predictSample(x, node->left.get());
            } else {
                return predictSample(x, node->right.get());
            }
        }
        
    public:
        DecisionTree(int max_depth = 5, int min_samples = 2) 
            : maxDepth(max_depth), minSamplesSplit(min_samples) {}
        
        void fit(const std::vector<std::vector<double>>& X, const std::vector<int>& y) {
            std::vector<int> indices(X.size());
            for (size_t i = 0; i < X.size(); i++) {
                indices[i] = i;
            }
            root = buildTree(X, y, indices, 0);
        }
        
        std::vector<int> predict(const std::vector<std::vector<double>>& X) {
            std::vector<int> predictions(X.size());
            for (size_t i = 0; i < X.size(); i++) {
                predictions[i] = predictSample(X[i], root.get());
            }
            return predictions;
        }
        
        double accuracy(const std::vector<int>& yTrue, const std::vector<int>& yPred) {
            int correct = 0;
            for (size_t i = 0; i < yTrue.size(); i++) {
                if (yTrue[i] == yPred[i]) correct++;
            }
            return static_cast<double>(correct) / yTrue.size();
        }
    };
    
    // General demonstration function
    void calc() {
        std::cout << "=== PCA (Principal Component Analysis) ===" << std::endl;
        
        // Sample data
        std::vector<std::vector<double>> data = {
            {2.5, 2.4}, {0.5, 0.7}, {2.2, 2.9}, {1.9, 2.2}, {3.1, 3.0},
            {2.3, 2.7}, {2.0, 1.6}, {1.0, 1.1}, {1.5, 1.6}, {1.1, 0.9}
        };
        
        PCA pca(2);
        pca.fit(data);
        auto transformed = pca.transform(data);
        
        std::cout << "Original data -> PCA transformed (first 3 points):" << std::endl;
        for (int i = 0; i < 3; i++) {
            std::cout << "(" << data[i][0] << ", " << data[i][1] << ") -> ";
            std::cout << "(" << transformed[i][0] << ", " << transformed[i][1] << ")" << std::endl;
        }
        std::cout << "Explained variance ratio PC1: " << pca.explainedVarianceRatio(0) << std::endl;
        
        std::cout << "\n=== K-Means Clustering ===" << std::endl;
        
        KMeans kmeans(2);
        kmeans.fit(data);
        auto labels = kmeans.getLabels();
        
        std::cout << "Cluster assignments: ";
        for (int label : labels) {
            std::cout << label << " ";
        }
        std::cout << std::endl;
        
        auto centroids = kmeans.getCentroids();
        std::cout << "Centroids:" << std::endl;
        for (size_t i = 0; i < centroids.size(); i++) {
            std::cout << "Cluster " << i << ": (" << centroids[i][0] << ", " << centroids[i][1] << ")" << std::endl;
        }
        
        std::cout << "\n=== Decision Tree ===" << std::endl;
        
        // Classification data
        std::vector<std::vector<double>> X = {
            {2.5, 1.5}, {3.5, 1.0}, {1.0, 2.0}, {1.5, 2.5},
            {3.0, 3.0}, {2.0, 1.0}, {4.0, 2.0}, {0.5, 1.5}
        };
        std::vector<int> y = {0, 0, 1, 1, 0, 1, 0, 1};  // Binary classification
        
        DecisionTree dt(3);
        dt.fit(X, y);
        auto predictions = dt.predict(X);
        
        std::cout << "True labels:  ";
        for (int label : y) std::cout << label << " ";
        std::cout << "\nPredictions:  ";
        for (int pred : predictions) std::cout << pred << " ";
        std::cout << "\nAccuracy: " << dt.accuracy(y, predictions) << std::endl;
    }
}

int main() {
    MachineLearningUtilities::calc();
    return 0;
}