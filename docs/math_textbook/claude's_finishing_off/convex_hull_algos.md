#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>
#include <cmath>

namespace ConvexHullAlgorithms {
    
    // 2D Point structure
    struct Point2D {
        double x, y;
        int index;  // Original index in input
        
        Point2D(double x = 0, double y = 0, int idx = -1) : x(x), y(y), index(idx) {}
        
        bool operator<(const Point2D& p) const {
            return (x < p.x) || (x == p.x && y < p.y);
        }
        
        bool operator==(const Point2D& p) const {
            return x == p.x && y == p.y;
        }
    };
    
    // Cross product of vectors OA and OB where O is origin
    double crossProduct(const Point2D& O, const Point2D& A, const Point2D& B) {
        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
    }
    
    // Distance squared between two points
    double distanceSquared(const Point2D& A, const Point2D& B) {
        double dx = A.x - B.x;
        double dy = A.y - B.y;
        return dx * dx + dy * dy;
    }
    
    // Graham Scan Algorithm for Convex Hull
    class GrahamScan {
    private:
        Point2D pivot;
        
        // Compare points by polar angle with respect to pivot
        bool polarAngleCompare(const Point2D& A, const Point2D& B) {
            double cross = crossProduct(pivot, A, B);
            if (cross == 0) {
                // Collinear points - keep the one farther from pivot
                return distanceSquared(pivot, A) < distanceSquared(pivot, B);
            }
            return cross > 0;  // Counter-clockwise
        }
        
    public:
        std::vector<Point2D> computeHull(std::vector<Point2D> points) {
            int n = points.size();
            if (n < 3) return points;
            
            // Find the bottom-most point (and leftmost if tie)
            int minIdx = 0;
            for (int i = 1; i < n; i++) {
                if (points[i].y < points[minIdx].y || 
                    (points[i].y == points[minIdx].y && points[i].x < points[minIdx].x)) {
                    minIdx = i;
                }
            }
            
            std::swap(points[0], points[minIdx]);
            pivot = points[0];
            
            // Sort points by polar angle with respect to pivot
            std::sort(points.begin() + 1, points.end(), 
                     [this](const Point2D& a, const Point2D& b) {
                         return polarAngleCompare(a, b);
                     });
            
            // Remove points with same angle but keep farthest
            int m = 1;
            for (int i = 1; i < n; i++) {
                while (i < n - 1 && crossProduct(pivot, points[i], points[i + 1]) == 0) {
                    i++;
                }
                points[m++] = points[i];
            }
            
            if (m < 3) return points;
            
            // Build convex hull using stack
            std::stack<Point2D> hull;
            hull.push(points[0]);
            hull.push(points[1]);
            hull.push(points[2]);
            
            for (int i = 3; i < m; i++) {
                // Remove points that make clockwise turn
                Point2D top = hull.top();
                hull.pop();
                while (!hull.empty() && crossProduct(hull.top(), top, points[i]) <= 0) {
                    top = hull.top();
                    hull.pop();
                }
                hull.push(top);
                hull.push(points[i]);
            }
            
            // Convert stack to vector
            std::vector<Point2D> result;
            while (!hull.empty()) {
                result.push_back(hull.top());
                hull.pop();
            }
            std::reverse(result.begin(), result.end());
            
            return result;
        }
    };
    
    // Jarvis March (Gift Wrapping) Algorithm for Convex Hull
    class JarvisMarch {
    public:
        std::vector<Point2D> computeHull(const std::vector<Point2D>& points) {
            int n = points.size();
            if (n < 3) return points;
            
            std::vector<Point2D> hull;
            
            // Find leftmost point
            int leftmost = 0;
            for (int i = 1; i < n; i++) {
                if (points[i].x < points[leftmost].x) {
                    leftmost = i;
                }
            }
            
            int p = leftmost, q;
            do {
                hull.push_back(points[p]);
                q = (p + 1) % n;
                
                for (int i = 0; i < n; i++) {
                    // Find the most counter-clockwise point from points[p]
                    if (crossProduct(points[p], points[i], points[q]) > 0) {
                        q = i;
                    }
                }
                
                p = q;
            } while (p != leftmost);
            
            return hull;
        }
    };
    
    // QuickHull Algorithm for Convex Hull
    class QuickHull {
    private:
        // Find point farthest from line AB on the left side
        int findFarthestPoint(const std::vector<Point2D>& points,
                              const std::vector<int>& indices,
                              const Point2D& A, const Point2D& B) {
            double maxDist = 0;
            int maxIdx = -1;
            
            for (int idx : indices) {
                double dist = std::abs(crossProduct(A, B, points[idx]));
                if (dist > maxDist) {
                    maxDist = dist;
                    maxIdx = idx;
                }
            }
            
            return maxIdx;
        }
        
        // Recursive hull finding
        void quickHullRecursive(const std::vector<Point2D>& points,
                               const std::vector<int>& indices,
                               const Point2D& A, const Point2D& B,
                               std::vector<Point2D>& hull) {
            if (indices.empty()) return;
            
            int farthestIdx = findFarthestPoint(points, indices, A, B);
            if (farthestIdx == -1) return;
            
            Point2D C = points[farthestIdx];
            
            // Find points on left of AC and CB
            std::vector<int> leftAC, leftCB;
            for (int idx : indices) {
                if (idx == farthestIdx) continue;
                
                if (crossProduct(A, C, points[idx]) > 0) {
                    leftAC.push_back(idx);
                } else if (crossProduct(C, B, points[idx]) > 0) {
                    leftCB.push_back(idx);
                }
            }
            
            quickHullRecursive(points, leftAC, A, C, hull);
            hull.push_back(C);
            quickHullRecursive(points, leftCB, C, B, hull);
        }
        
    public:
        std::vector<Point2D> computeHull(const std::vector<Point2D>& points) {
            int n = points.size();
            if (n < 3) return points;
            
            // Find leftmost and rightmost points
            int minIdx = 0, maxIdx = 0;
            for (int i = 1; i < n; i++) {
                if (points[i].x < points[minIdx].x) minIdx = i;
                if (points[i].x > points[maxIdx].x) maxIdx = i;
            }
            
            Point2D A = points[minIdx];
            Point2D B = points[maxIdx];
            
            std::vector<Point2D> hull;
            hull.push_back(A);
            
            // Divide points into upper and lower sets
            std::vector<int> upper, lower;
            for (int i = 0; i < n; i++) {
                if (i == minIdx || i == maxIdx) continue;
                
                double cross = crossProduct(A, B, points[i]);
                if (cross > 0) {
                    upper.push_back(i);
                } else if (cross < 0) {
                    lower.push_back(i);
                }
            }
            
            // Find upper hull
            quickHullRecursive(points, upper, A, B, hull);
            hull.push_back(B);
            
            // Find lower hull
            quickHullRecursive(points, lower, B, A, hull);
            
            return hull;
        }
    };
    
    // 3D Convex Hull using QuickHull algorithm
    struct Point3D {
        double x, y, z;
        int index;
        
        Point3D(double x = 0, double y = 0, double z = 0, int idx = -1) 
            : x(x), y(y), z(z), index(idx) {}
        
        Point3D operator-(const Point3D& p) const {
            return Point3D(x - p.x, y - p.y, z - p.z);
        }
        
        Point3D cross(const Point3D& p) const {
            return Point3D(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
        }
        
        double dot(const Point3D& p) const {
            return x * p.x + y * p.y + z * p.z;
        }
    };
    
    struct Face3D {
        int a, b, c;  // Indices of vertices
        Point3D normal;
        
        Face3D(int a, int b, int c) : a(a), b(b), c(c) {}
    };
    
    // Incremental 3D Convex Hull
    class ConvexHull3D {
    public:
        std::vector<Face3D> computeHull(const std::vector<Point3D>& points) {
            int n = points.size();
            if (n < 4) return {};  // Need at least 4 points for 3D hull
            
            // Find initial tetrahedron (simplified - just use first 4 non-coplanar points)
            std::vector<Face3D> faces;
            faces.push_back(Face3D(0, 1, 2));
            faces.push_back(Face3D(0, 2, 3));
            faces.push_back(Face3D(0, 3, 1));
            faces.push_back(Face3D(1, 3, 2));
            
            // Compute normals for initial faces
            for (auto& face : faces) {
                Point3D v1 = points[face.b] - points[face.a];
                Point3D v2 = points[face.c] - points[face.a];
                face.normal = v1.cross(v2);
            }
            
            // Add remaining points incrementally (simplified version)
            // Full implementation would include horizon edge detection and face updates
            
            return faces;
        }
    };
    
    // General demonstration function
    void calc() {
        std::cout << "=== Convex Hull Algorithms ===" << std::endl;
        
        // Test points
        std::vector<Point2D> points = {
            Point2D(0, 0), Point2D(1, 1), Point2D(2, 0), Point2D(2, 2),
            Point2D(0, 2), Point2D(1, 0.5), Point2D(1.5, 1.5), Point2D(0.5, 1.8)
        };
        
        std::cout << "Input points: ";
        for (const auto& p : points) {
            std::cout << "(" << p.x << "," << p.y << ") ";
        }
        std::cout << std::endl;
        
        // Graham Scan
        std::cout << "\nGraham Scan Convex Hull:" << std::endl;
        GrahamScan graham;
        auto hullGraham = graham.computeHull(points);
        std::cout << "Hull vertices: ";
        for (const auto& p : hullGraham) {
            std::cout << "(" << p.x << "," << p.y << ") ";
        }
        std::cout << std::endl;
        
        // Jarvis March
        std::cout << "\nJarvis March Convex Hull:" << std::endl;
        JarvisMarch jarvis;
        auto hullJarvis = jarvis.computeHull(points);
        std::cout << "Hull vertices: ";
        for (const auto& p : hullJarvis) {
            std::cout << "(" << p.x << "," << p.y << ") ";
        }
        std::cout << std::endl;
        
        // QuickHull
        std::cout << "\nQuickHull Convex Hull:" << std::endl;
        QuickHull quick;
        auto hullQuick = quick.computeHull(points);
        std::cout << "Hull vertices: ";
        for (const auto& p : hullQuick) {
            std::cout << "(" << p.x << "," << p.y << ") ";
        }
        std::cout << std::endl;
        
        // 3D Convex Hull
        std::cout << "\n3D Convex Hull (simplified):" << std::endl;
        std::vector<Point3D> points3D = {
            Point3D(0, 0, 0), Point3D(1, 0, 0), Point3D(0, 1, 0),
            Point3D(0, 0, 1), Point3D(0.5, 0.5, 0.5)
        };
        
        ConvexHull3D hull3D;
        auto faces = hull3D.computeHull(points3D);
        std::cout << "Hull has " << faces.size() << " faces" << std::endl;
    }
}

int main() {
    ConvexHullAlgorithms::calc();
    return 0;
}