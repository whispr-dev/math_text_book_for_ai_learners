#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <array>
#include <functional>
#include <random>

namespace AdvancedGeometry {

    const double PI = 3.141592653589793;
    const double EPSILON = 1e-12;
    const double DEG_TO_RAD = PI / 180.0;
    const double RAD_TO_DEG = 180.0 / PI;
    
    // ===== FUNDAMENTAL GEOMETRIC STRUCTURES =====
    
    struct Point2D {
        double x, y;
        
        Point2D(double x = 0.0, double y = 0.0) : x(x), y(y) {}
        
        Point2D operator+(const Point2D& p) const { return Point2D(x + p.x, y + p.y); }
        Point2D operator-(const Point2D& p) const { return Point2D(x - p.x, y - p.y); }
        Point2D operator*(double s) const { return Point2D(x * s, y * s); }
        Point2D operator/(double s) const { return Point2D(x / s, y / s); }
        
        double dot(const Point2D& p) const { return x * p.x + y * p.y; }
        double cross(const Point2D& p) const { return x * p.y - y * p.x; }
        double magnitude() const { return std::sqrt(x * x + y * y); }
        double magnitudeSquared() const { return x * x + y * y; }
        Point2D normalize() const { double mag = magnitude(); return mag > 0 ? *this / mag : Point2D(0, 0); }
        
        void print() const { std::cout << "(" << x << ", " << y << ")"; }
    };
    
    struct Point3D {
        double x, y, z;
        
        Point3D(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}
        
        Point3D operator+(const Point3D& p) const { return Point3D(x + p.x, y + p.y, z + p.z); }
        Point3D operator-(const Point3D& p) const { return Point3D(x - p.x, y - p.y, z - p.z); }
        Point3D operator*(double s) const { return Point3D(x * s, y * s, z * s); }
        Point3D operator/(double s) const { return Point3D(x / s, y / s, z / s); }
        
        double dot(const Point3D& p) const { return x * p.x + y * p.y + z * p.z; }
        Point3D cross(const Point3D& p) const { 
            return Point3D(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x); 
        }
        double magnitude() const { return std::sqrt(x * x + y * y + z * z); }
        double magnitudeSquared() const { return x * x + y * y + z * z; }
        Point3D normalize() const { double mag = magnitude(); return mag > 0 ? *this / mag : Point3D(0, 0, 0); }
        
        // Convert to spherical coordinates (r, theta, phi)
        std::array<double, 3> toSpherical() const {
            double r = magnitude();
            double theta = std::atan2(y, x); // Azimuth
            double phi = r > 0 ? std::acos(z / r) : 0; // Elevation from z-axis
            return {r, theta, phi};
        }
        
        // Create from spherical coordinates
        static Point3D fromSpherical(double r, double theta, double phi) {
            return Point3D(
                r * std::sin(phi) * std::cos(theta),
                r * std::sin(phi) * std::sin(theta),
                r * std::cos(phi)
            );
        }
        
        void print() const { std::cout << "(" << x << ", " << y << ", " << z << ")"; }
    };
    
    struct Matrix3x3 {
        std::array<std::array<double, 3>, 3> m;
        
        Matrix3x3() {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    m[i][j] = (i == j) ? 1.0 : 0.0; // Identity matrix
                }
            }
        }
        
        Matrix3x3(const std::array<std::array<double, 3>, 3>& matrix) : m(matrix) {}
        
        Point3D multiply(const Point3D& p) const {
            return Point3D(
                m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z,
                m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z,
                m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z
            );
        }
        
        Matrix3x3 multiply(const Matrix3x3& other) const {
            Matrix3x3 result;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    result.m[i][j] = 0;
                    for (int k = 0; k < 3; ++k) {
                        result.m[i][j] += m[i][k] * other.m[k][j];
                    }
                }
            }
            return result;
        }
        
        double determinant() const {
            return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                   m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                   m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        }
        
        Matrix3x3 transpose() const {
            Matrix3x3 result;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    result.m[i][j] = m[j][i];
                }
            }
            return result;
        }
        
        void print() const {
            for (int i = 0; i < 3; ++i) {
                std::cout << "[";
                for (int j = 0; j < 3; ++j) {
                    std::cout << m[i][j];
                    if (j < 2) std::cout << ", ";
                }
                std::cout << "]" << std::endl;
            }
        }
    };
    
    // ===== SPHERICAL TRIGONOMETRY =====
    
    struct SphericalTriangle {
        double a, b, c; // Side lengths (in radians)
        double A, B, C; // Angles (in radians)
        
        SphericalTriangle(double a = 0, double b = 0, double c = 0) : a(a), b(b), c(c), A(0), B(0), C(0) {}
        
        void print() const {
            std::cout << "Spherical Triangle:" << std::endl;
            std::cout << "Sides (rad): a=" << a << ", b=" << b << ", c=" << c << std::endl;
            std::cout << "Angles (rad): A=" << A << ", B=" << B << ", C=" << C << std::endl;
            std::cout << "Sides (deg): a=" << a*RAD_TO_DEG << ", b=" << b*RAD_TO_DEG << ", c=" << c*RAD_TO_DEG << std::endl;
            std::cout << "Angles (deg): A=" << A*RAD_TO_DEG << ", B=" << B*RAD_TO_DEG << ", C=" << C*RAD_TO_DEG << std::endl;
        }
    };
    
    // Great circle distance between two points on a sphere (Haversine formula)
    double greatCircleDistance(const Point3D& p1, const Point3D& p2, double radius = 1.0) {
        auto [r1, theta1, phi1] = p1.toSpherical();
        auto [r2, theta2, phi2] = p2.toSpherical();
        
        double lat1 = PI/2 - phi1; // Convert to latitude
        double lat2 = PI/2 - phi2;
        double dlon = theta2 - theta1;
        
        double a = std::sin((lat2 - lat1) / 2) * std::sin((lat2 - lat1) / 2) +
                  std::cos(lat1) * std::cos(lat2) * 
                  std::sin(dlon / 2) * std::sin(dlon / 2);
        double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
        
        return radius * c;
    }
    
    // Solve spherical triangle using Law of Cosines
    SphericalTriangle solveSphericalTriangle(double a, double b, double c) {
        SphericalTriangle triangle(a, b, c);
        
        // Law of cosines for spherical trigonometry
        // cos(a) = cos(b)cos(c) + sin(b)sin(c)cos(A)
        triangle.A = std::acos((std::cos(a) - std::cos(b) * std::cos(c)) / 
                              (std::sin(b) * std::sin(c)));
        triangle.B = std::acos((std::cos(b) - std::cos(c) * std::cos(a)) / 
                              (std::sin(c) * std::sin(a)));
        triangle.C = std::acos((std::cos(c) - std::cos(a) * std::cos(b)) / 
                              (std::sin(a) * std::sin(b)));
        
        return triangle;
    }
    
    // Calculate spherical area (spherical excess)
    double sphericalTriangleArea(const SphericalTriangle& triangle, double radius = 1.0) {
        double excess = triangle.A + triangle.B + triangle.C - PI;
        return excess * radius * radius;
    }
    
    // ===== GEOMETRIC TRANSFORMATIONS =====
    
    // 2D Rotation matrix
    Matrix3x3 rotation2D(double angle) {
        double cos_a = std::cos(angle);
        double sin_a = std::sin(angle);
        return Matrix3x3({{
            {{cos_a, -sin_a, 0}},
            {{sin_a, cos_a, 0}},
            {{0, 0, 1}}
        }});
    }
    
    // 3D Rotation matrices
    Matrix3x3 rotationX(double angle) {
        double cos_a = std::cos(angle);
        double sin_a = std::sin(angle);
        return Matrix3x3({{
            {{1, 0, 0}},
            {{0, cos_a, -sin_a}},
            {{0, sin_a, cos_a}}
        }});
    }
    
    Matrix3x3 rotationY(double angle) {
        double cos_a = std::cos(angle);
        double sin_a = std::sin(angle);
        return Matrix3x3({{
            {{cos_a, 0, sin_a}},
            {{0, 1, 0}},
            {{-sin_a, 0, cos_a}}
        }});
    }
    
    Matrix3x3 rotationZ(double angle) {
        double cos_a = std::cos(angle);
        double sin_a = std::sin(angle);
        return Matrix3x3({{
            {{cos_a, -sin_a, 0}},
            {{sin_a, cos_a, 0}},
            {{0, 0, 1}}
        }});
    }
    
    // Euler angle rotation (ZYX convention)
    Matrix3x3 eulerRotation(double yaw, double pitch, double roll) {
        return rotationZ(yaw).multiply(rotationY(pitch)).multiply(rotationX(roll));
    }
    
    // Quaternion representation for rotations
    struct Quaternion {
        double w, x, y, z;
        
        Quaternion(double w = 1, double x = 0, double y = 0, double z = 0) : w(w), x(x), y(y), z(z) {}
        
        // Create from axis-angle representation
        static Quaternion fromAxisAngle(const Point3D& axis, double angle) {
            Point3D normalized_axis = axis.normalize();
            double half_angle = angle / 2.0;
            double sin_half = std::sin(half_angle);
            
            return Quaternion(
                std::cos(half_angle),
                normalized_axis.x * sin_half,
                normalized_axis.y * sin_half,
                normalized_axis.z * sin_half
            );
        }
        
        Quaternion multiply(const Quaternion& q) const {
            return Quaternion(
                w * q.w - x * q.x - y * q.y - z * q.z,
                w * q.x + x * q.w + y * q.z - z * q.y,
                w * q.y - x * q.z + y * q.w + z * q.x,
                w * q.z + x * q.y - y * q.x + z * q.w
            );
        }
        
        Quaternion conjugate() const {
            return Quaternion(w, -x, -y, -z);
        }
        
        double magnitude() const {
            return std::sqrt(w * w + x * x + y * y + z * z);
        }
        
        Quaternion normalize() const {
            double mag = magnitude();
            return mag > 0 ? Quaternion(w/mag, x/mag, y/mag, z/mag) : Quaternion();
        }
        
        Point3D rotatePoint(const Point3D& p) const {
            Quaternion point_quat(0, p.x, p.y, p.z);
            Quaternion result = multiply(point_quat).multiply(conjugate());
            return Point3D(result.x, result.y, result.z);
        }
        
        Matrix3x3 toMatrix() const {
            double xx = x * x, xy = x * y, xz = x * z, xw = x * w;
            double yy = y * y, yz = y * z, yw = y * w;
            double zz = z * z, zw = z * w;
            
            return Matrix3x3({{
                {{1 - 2*(yy + zz), 2*(xy - zw), 2*(xz + yw)}},
                {{2*(xy + zw), 1 - 2*(xx + zz), 2*(yz - xw)}},
                {{2*(xz - yw), 2*(yz + xw), 1 - 2*(xx + yy)}}
            }});
        }
        
        void print() const {
            std::cout << "Quaternion(w=" << w << ", x=" << x << ", y=" << y << ", z=" << z << ")";
        }
    };
    
    // ===== CONVEX HULL ALGORITHMS =====
    
    // Cross product for 2D points (z-component)
    double crossProduct2D(const Point2D& O, const Point2D& A, const Point2D& B) {
        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
    }
    
    // Andrew's monotone chain convex hull algorithm
    std::vector<Point2D> convexHull2D(std::vector<Point2D> points) {
        int n = points.size();
        if (n <= 1) return points;
        
        // Sort points lexicographically
        std::sort(points.begin(), points.end(), [](const Point2D& a, const Point2D& b) {
            return a.x < b.x || (a.x == b.x && a.y < b.y);
        });
        
        // Build lower hull
        std::vector<Point2D> hull;
        for (int i = 0; i < n; ++i) {
            while (hull.size() >= 2 && crossProduct2D(hull[hull.size()-2], hull[hull.size()-1], points[i]) <= 0) {
                hull.pop_back();
            }
            hull.push_back(points[i]);
        }
        
        // Build upper hull
        int t = hull.size() + 1;
        for (int i = n - 2; i >= 0; --i) {
            while (hull.size() >= t && crossProduct2D(hull[hull.size()-2], hull[hull.size()-1], points[i]) <= 0) {
                hull.pop_back();
            }
            hull.push_back(points[i]);
        }
        
        hull.pop_back(); // Remove last point as it's same as first
        return hull;
    }
    
    // Gift wrapping (Jarvis march) algorithm for 3D convex hull (simplified)
    std::vector<Point3D> convexHull3D_Simple(const std::vector<Point3D>& points) {
        // Simplified 3D convex hull - find extreme points
        if (points.size() < 4) return points;
        
        std::vector<Point3D> hull;
        
        // Find extreme points in each direction
        auto min_x = *std::min_element(points.begin(), points.end(), 
            [](const Point3D& a, const Point3D& b) { return a.x < b.x; });
        auto max_x = *std::max_element(points.begin(), points.end(), 
            [](const Point3D& a, const Point3D& b) { return a.x < b.x; });
        auto min_y = *std::min_element(points.begin(), points.end(), 
            [](const Point3D& a, const Point3D& b) { return a.y < b.y; });
        auto max_y = *std::max_element(points.begin(), points.end(), 
            [](const Point3D& a, const Point3D& b) { return a.y < b.y; });
        auto min_z = *std::min_element(points.begin(), points.end(), 
            [](const Point3D& a, const Point3D& b) { return a.z < b.z; });
        auto max_z = *std::max_element(points.begin(), points.end(), 
            [](const Point3D& a, const Point3D& b) { return a.z < b.z; });
        
        hull = {min_x, max_x, min_y, max_y, min_z, max_z};
        
        // Remove duplicates
        std::sort(hull.begin(), hull.end(), [](const Point3D& a, const Point3D& b) {
            return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z);
        });
        hull.erase(std::unique(hull.begin(), hull.end(), [](const Point3D& a, const Point3D& b) {
            return std::abs(a.x - b.x) < EPSILON && std::abs(a.y - b.y) < EPSILON && std::abs(a.z - b.z) < EPSILON;
        }), hull.end());
        
        return hull;
    }
    
    // ===== ADVANCED GEOMETRIC ALGORITHMS =====
    
    // Point in polygon test (ray casting)
    bool pointInPolygon(const Point2D& point, const std::vector<Point2D>& polygon) {
        int n = polygon.size();
        bool inside = false;
        
        for (int i = 0, j = n - 1; i < n; j = i++) {
            if (((polygon[i].y > point.y) != (polygon[j].y > point.y)) &&
                (point.x < (polygon[j].x - polygon[i].x) * (point.y - polygon[i].y) / 
                 (polygon[j].y - polygon[i].y) + polygon[i].x)) {
                inside = !inside;
            }
        }
        
        return inside;
    }
    
    // Closest pair of points using divide and conquer
    struct PointPair {
        Point2D p1, p2;
        double distance;
        
        void print() const {
            std::cout << "Closest pair: ";
            p1.print();
            std::cout << " and ";
            p2.print();
            std::cout << " (distance: " << distance << ")";
        }
    };
    
    PointPair closestPairOfPoints(std::vector<Point2D> points) {
        std::sort(points.begin(), points.end(), [](const Point2D& a, const Point2D& b) {
            return a.x < b.x;
        });
        
        std::function<PointPair(std::vector<Point2D>&, int, int)> closestPairRec = 
            [&](std::vector<Point2D>& px, int left, int right) -> PointPair {
            
            int n = right - left + 1;
            
            // Base case: brute force for small arrays
            if (n <= 3) {
                PointPair result;
                result.distance = std::numeric_limits<double>::infinity();
                
                for (int i = left; i <= right; ++i) {
                    for (int j = i + 1; j <= right; ++j) {
                        double dist = (px[i] - px[j]).magnitude();
                        if (dist < result.distance) {
                            result.distance = dist;
                            result.p1 = px[i];
                            result.p2 = px[j];
                        }
                    }
                }
                return result;
            }
            
            // Divide
            int mid = (left + right) / 2;
            auto left_result = closestPairRec(px, left, mid);
            auto right_result = closestPairRec(px, mid + 1, right);
            
            // Find minimum of the two
            PointPair min_result = (left_result.distance < right_result.distance) ? 
                                  left_result : right_result;
            
            // Find points close to the dividing line
            std::vector<Point2D> strip;
            for (int i = left; i <= right; ++i) {
                if (std::abs(px[i].x - px[mid].x) < min_result.distance) {
                    strip.push_back(px[i]);
                }
            }
            
            // Sort strip by y-coordinate
            std::sort(strip.begin(), strip.end(), [](const Point2D& a, const Point2D& b) {
                return a.y < b.y;
            });
            
            // Check points in strip
            for (size_t i = 0; i < strip.size(); ++i) {
                for (size_t j = i + 1; j < strip.size() && 
                     (strip[j].y - strip[i].y) < min_result.distance; ++j) {
                    double dist = (strip[i] - strip[j]).magnitude();
                    if (dist < min_result.distance) {
                        min_result.distance = dist;
                        min_result.p1 = strip[i];
                        min_result.p2 = strip[j];
                    }
                }
            }
            
            return min_result;
        };
        
        return closestPairRec(points, 0, points.size() - 1);
    }
    
    // Delaunay triangulation (simplified - returns triangle centers)
    std::vector<Point2D> delaunayTriangulation(const std::vector<Point2D>& points) {
        // Simplified implementation - returns circumcenters of triangles
        std::vector<Point2D> centers;
        
        // For simplicity, create triangles from every 3 consecutive points on convex hull
        auto hull = convexHull2D(const_cast<std::vector<Point2D>&>(points));
        
        for (size_t i = 0; i < hull.size(); ++i) {
            if (i + 2 < hull.size()) {
                Point2D a = hull[i];
                Point2D b = hull[i + 1];
                Point2D c = hull[i + 2];
                
                // Calculate circumcenter
                double d = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
                if (std::abs(d) > EPSILON) {
                    double ux = ((a.x*a.x + a.y*a.y) * (b.y - c.y) + 
                                (b.x*b.x + b.y*b.y) * (c.y - a.y) + 
                                (c.x*c.x + c.y*c.y) * (a.y - b.y)) / d;
                    double uy = ((a.x*a.x + a.y*a.y) * (c.x - b.x) + 
                                (b.x*b.x + b.y*b.y) * (a.x - c.x) + 
                                (c.x*c.x + c.y*c.y) * (b.x - a.x)) / d;
                    centers.push_back(Point2D(ux, uy));
                }
            }
        }
        
        return centers;
    }
    
    // ===== 3D AUDIO GEOMETRY =====
    
    struct AudioSource {
        Point3D position;
        double intensity;
        
        AudioSource(const Point3D& pos, double intens = 1.0) : position(pos), intensity(intens) {}
    };
    
    struct AudioListener {
        Point3D position;
        Point3D forward_direction;
        Point3D up_direction;
        
        AudioListener(const Point3D& pos, const Point3D& forward = Point3D(0, 0, -1), 
                     const Point3D& up = Point3D(0, 1, 0)) 
            : position(pos), forward_direction(forward.normalize()), up_direction(up.normalize()) {}
    };
    
    // Calculate 3D audio parameters
    struct AudioParams {
        double distance;
        double attenuation;
        double azimuth;    // Horizontal angle
        double elevation;  // Vertical angle
        double doppler_factor;
        
        void print() const {
            std::cout << "Audio Parameters:" << std::endl;
            std::cout << "Distance: " << distance << std::endl;
            std::cout << "Attenuation: " << attenuation << std::endl;
            std::cout << "Azimuth: " << azimuth * RAD_TO_DEG << "°" << std::endl;
            std::cout << "Elevation: " << elevation * RAD_TO_DEG << "°" << std::endl;
            std::cout << "Doppler factor: " << doppler_factor << std::endl;
        }
    };
    
    AudioParams calculate3DAudio(const AudioSource& source, const AudioListener& listener,
                                const Point3D& source_velocity = Point3D(0,0,0),
                                const Point3D& listener_velocity = Point3D(0,0,0),
                                double sound_speed = 343.0) { // m/s at 20°C
        
        Point3D source_to_listener = listener.position - source.position;
        double distance = source_to_listener.magnitude();
        
        // Distance-based attenuation (inverse square law)
        double attenuation = source.intensity / (distance * distance + 1.0);
        
        // Convert to listener's local coordinate system
        Point3D right = listener.forward_direction.cross(listener.up_direction);
        Point3D local_source = Point3D(
            source_to_listener.dot(right),
            source_to_listener.dot(listener.up_direction),
            -source_to_listener.dot(listener.forward_direction) // Forward is negative z
        ).normalize();
        
        // Calculate spherical coordinates in listener space
        auto [r, azimuth, elevation_from_z] = local_source.toSpherical();
        double elevation = PI/2 - elevation_from_z; // Convert to elevation from horizontal
        
        // Doppler effect calculation
        Point3D relative_velocity = source_velocity - listener_velocity;
        Point3D direction_to_listener = source_to_listener.normalize();
        double radial_velocity = relative_velocity.dot(direction_to_listener);
        double doppler_factor = (sound_speed - radial_velocity) / sound_speed;
        
        return {distance, attenuation, azimuth, elevation, doppler_factor};
    }
    
    // Room acoustics - calculate reflection points on walls
    std::vector<Point3D> calculateReflectionPoints(const AudioSource& source, 
                                                   const AudioListener& listener,
                                                   const std::vector<Point3D>& room_corners) {
        std::vector<Point3D> reflections;
        
        // Simplified room acoustics - find reflections on axis-aligned walls
        if (room_corners.size() >= 8) { // Assume rectangular room
            // Find room bounds
            double min_x = room_corners[0].x, max_x = room_corners[0].x;
            double min_y = room_corners[0].y, max_y = room_corners[0].y;
            double min_z = room_corners[0].z, max_z = room_corners[0].z;
            
            for (const auto& corner : room_corners) {
                min_x = std::min(min_x, corner.x); max_x = std::max(max_x, corner.x);
                min_y = std::min(min_y, corner.y); max_y = std::max(max_y, corner.y);
                min_z = std::min(min_z, corner.z); max_z = std::max(max_z, corner.z);
            }
            
            // Calculate first-order reflections
            std::vector<std::pair<Point3D, Point3D>> walls = {
                {Point3D(min_x, 0, 0), Point3D(1, 0, 0)}, // Left wall
                {Point3D(max_x, 0, 0), Point3D(-1, 0, 0)}, // Right wall
                {Point3D(0, min_y, 0), Point3D(0, 1, 0)}, // Floor
                {Point3D(0, max_y, 0), Point3D(0, -1, 0)}, // Ceiling
                {Point3D(0, 0, min_z), Point3D(0, 0, 1)}, // Back wall
                {Point3D(0, 0, max_z), Point3D(0, 0, -1)}  // Front wall
            };
            
            for (const auto& [wall_point, wall_normal] : walls) {
                // Mirror source across wall
                Point3D to_wall = wall_point - source.position;
                double distance_to_wall = to_wall.dot(wall_normal);
                Point3D mirrored_source = source.position + wall_normal * (2.0 * distance_to_wall);
                
                // Find intersection point on wall
                Point3D wall_to_listener = listener.position - wall_point;
                Point3D wall_to_mirrored = mirrored_source - wall_point;
                
                // Simple reflection point calculation
                Point3D reflection_point = wall_point + (wall_to_listener + wall_to_mirrored) * 0.5;
                reflections.push_back(reflection_point);
            }
        }
        
        return reflections;
    }
    
    // ===== WAVE GEOMETRY AND INTERFERENCE =====
    
    struct Wave {
        Point3D origin;
        double frequency;
        double amplitude;
        double phase;
        double wavelength;
        
        Wave(const Point3D& orig, double freq, double amp = 1.0, double ph = 0.0, double sound_speed = 343.0) 
            : origin(orig), frequency(freq), amplitude(amp), phase(ph), wavelength(sound_speed / freq) {}
        
        // Calculate wave amplitude at a point and time
        double amplitudeAt(const Point3D& point, double time) const {
            double distance = (point - origin).magnitude();
            double wave_phase = 2.0 * PI * (frequency * time - distance / wavelength) + phase;
            return amplitude * std::sin(wave_phase) / (distance + 1.0); // Distance attenuation
        }
    };
    
    // Calculate interference pattern from multiple wave sources
    std::vector<std::vector<double>> calculateInterferencePattern(
        const std::vector<Wave>& waves,
        const Point2D& area_min, const Point2D& area_max,
        int grid_width, int grid_height,
        double time = 0.0, double z_plane = 0.0) {
        
        std::vector<std::vector<double>> pattern(grid_height, std::vector<double>(grid_width, 0.0));
        
        double dx = (area_max.x - area_min.x) / (grid_width - 1);
        double dy = (area_max.y - area_min.y) / (grid_height - 1);
        
        for (int i = 0; i < grid_height; ++i) {
            for (int j = 0; j < grid_width; ++j) {
                Point3D point(area_min.x + j * dx, area_min.y + i * dy, z_plane);
                
                double total_amplitude = 0.0;
                for (const auto& wave : waves) {
                    total_amplitude += wave.amplitudeAt(point, time);
                }
                
                pattern[i][j] = total_amplitude;
            }
        }
        
        return pattern;
    }
    
    // ===== ADVANCED GEOMETRIC COMPUTATIONS =====
    
    // Möller–Trumbore ray-triangle intersection
    struct RayTriangleIntersection {
        bool hit;
        double t; // Parameter along ray
        Point3D intersection_point;
        Point2D barycentric_coords;
        
        void print() const {
            if (hit) {
                std::cout << "Ray-triangle intersection at t=" << t << ", point: ";
                intersection_point.print();
                std::cout << ", barycentric: (" << barycentric_coords.x << ", " << barycentric_coords.y << ")";
            } else {
                std::cout << "No intersection";
            }
        }
    };
    
    RayTriangleIntersection rayTriangleIntersect(const Point3D& ray_origin, const Point3D& ray_direction,
                                               const Point3D& v0, const Point3D& v1, const Point3D& v2) {
        RayTriangleIntersection result;
        result.hit = false;
        
        Point3D edge1 = v1 - v0;
        Point3D edge2 = v2 - v0;
        Point3D h = ray_direction.cross(edge2);
        double a = edge1.dot(h);
        
        if (std::abs(a) < EPSILON) return result; // Ray parallel to triangle
        
        double f = 1.0 / a;
        Point3D s = ray_origin - v0;
        double u = f * s.dot(h);
        
        if (u < 0.0 || u > 1.0) return result;
        
        Point3D q = s.cross(edge1);
        double v = f * ray_direction.dot(q);
        
        if (v < 0.0 || u + v > 1.0) return result;
        
        double t = f * edge2.dot(q);
        
        if (t > EPSILON) {
            result.hit = true;
            result.t = t;
            result.intersection_point = ray_origin + ray_direction * t;
            result.barycentric_coords = Point2D(u, v);
        }
        
        return result;
    }
    
    // Bounding box calculations
    struct BoundingBox3D {
        Point3D min_corner, max_corner;
        
        BoundingBox3D(const Point3D& min_pt, const Point3D& max_pt) : min_corner(min_pt), max_corner(max_pt) {}
        
        static BoundingBox3D fromPoints(const std::vector<Point3D>& points) {
            if (points.empty()) return BoundingBox3D(Point3D(0,0,0), Point3D(0,0,0));
            
            Point3D min_pt = points[0], max_pt = points[0];
            for (const auto& p : points) {
                min_pt.x = std::min(min_pt.x, p.x); max_pt.x = std::max(max_pt.x, p.x);
                min_pt.y = std::min(min_pt.y, p.y); max_pt.y = std::max(max_pt.y, p.y);
                min_pt.z = std::min(min_pt.z, p.z); max_pt.z = std::max(max_pt.z, p.z);
            }
            return BoundingBox3D(min_pt, max_pt);
        }
        
        Point3D center() const { return (min_corner + max_corner) * 0.5; }
        Point3D size() const { return max_corner - min_corner; }
        double volume() const { 
            Point3D s = size(); 
            return s.x * s.y * s.z; 
        }
        
        bool contains(const Point3D& point) const {
            return point.x >= min_corner.x && point.x <= max_corner.x &&
                   point.y >= min_corner.y && point.y <= max_corner.y &&
                   point.z >= min_corner.z && point.z <= max_corner.z;
        }
        
        bool intersects(const BoundingBox3D& other) const {
            return !(max_corner.x < other.min_corner.x || min_corner.x > other.max_corner.x ||
                     max_corner.y < other.min_corner.y || min_corner.y > other.max_corner.y ||
                     max_corner.z < other.min_corner.z || min_corner.z > other.max_corner.z);
        }
        
        void print() const {
            std::cout << "BoundingBox3D: min=";
            min_corner.print();
            std::cout << ", max=";
            max_corner.print();
            std::cout << ", volume=" << volume();
        }
    };
    
    // ===== COMPUTATIONAL GEOMETRY UTILITIES =====
    
    // Generate random points in various distributions
    class GeometryGenerator {
    private:
        std::mt19937 rng;
        
    public:
        GeometryGenerator(unsigned seed = std::random_device{}()) : rng(seed) {}
        
        // Random points in circle
        std::vector<Point2D> randomPointsInCircle(int count, const Point2D& center, double radius) {
            std::vector<Point2D> points;
            std::uniform_real_distribution<double> angle_dist(0, 2 * PI);
            std::uniform_real_distribution<double> radius_dist(0, 1);
            
            for (int i = 0; i < count; ++i) {
                double angle = angle_dist(rng);
                double r = radius * std::sqrt(radius_dist(rng)); // Uniform distribution in circle
                points.push_back(Point2D(center.x + r * std::cos(angle), center.y + r * std::sin(angle)));
            }
            return points;
        }
        
        // Random points in sphere
        std::vector<Point3D> randomPointsInSphere(int count, const Point3D& center, double radius) {
            std::vector<Point3D> points;
            std::uniform_real_distribution<double> uniform(0, 1);
            std::normal_distribution<double> normal(0, 1);
            
            for (int i = 0; i < count; ++i) {
                // Generate uniform random point in sphere using normal distribution
                double x = normal(rng), y = normal(rng), z = normal(rng);
                double mag = std::sqrt(x*x + y*y + z*z);
                double scale = radius * std::cbrt(uniform(rng)) / mag;
                
                points.push_back(Point3D(center.x + x * scale, center.y + y * scale, center.z + z * scale));
            }
            return points;
        }
        
        // Random points on sphere surface
        std::vector<Point3D> randomPointsOnSphere(int count, const Point3D& center, double radius) {
            std::vector<Point3D> points;
            std::normal_distribution<double> normal(0, 1);
            
            for (int i = 0; i < count; ++i) {
                double x = normal(rng), y = normal(rng), z = normal(rng);
                double mag = std::sqrt(x*x + y*y + z*z);
                points.push_back(Point3D(center.x + radius * x/mag, center.y + radius * y/mag, center.z + radius * z/mag));
            }
            return points;
        }
    };
    
} // namespace AdvancedGeometry

// ===== DEMONSTRATION =====

int main() {
    using namespace AdvancedGeometry;
    
    std::cout << "=== ADVANCED TRIGONOMETRY AND GEOMETRY MODULE DEMO ===" << std::endl << std::endl;
    
    // ===== SPHERICAL TRIGONOMETRY =====
    std::cout << "1. SPHERICAL TRIGONOMETRY" << std::endl;
    
    // Calculate great circle distance between two cities
    Point3D london_sphere = Point3D::fromSpherical(1.0, -0.1278 * DEG_TO_RAD, (90 - 51.5074) * DEG_TO_RAD);
    Point3D tokyo_sphere = Point3D::fromSpherical(1.0, 139.6917 * DEG_TO_RAD, (90 - 35.6762) * DEG_TO_RAD);
    
    double distance_km = greatCircleDistance(london_sphere, tokyo_sphere, 6371.0);
    std::cout << "Great circle distance London to Tokyo: " << distance_km << " km" << std::endl;
    
    // Solve spherical triangle
    auto triangle = solveSphericalTriangle(PI/3, PI/4, PI/2); // 60°, 45°, 90° sides
    triangle.print();
    std::cout << "Triangle area: " << sphericalTriangleArea(triangle, 1.0) << " steradians" << std::endl;
    std::cout << std::endl;
    
    // ===== 3D TRANSFORMATIONS =====
    std::cout << "2. 3D GEOMETRIC TRANSFORMATIONS" << std::endl;
    
    Point3D original_point(1, 0, 0);
    std::cout << "Original point: ";
    original_point.print();
    std::cout << std::endl;
    
    // Euler rotations
    auto rotation_matrix = eulerRotation(PI/4, PI/6, PI/3); // 45°, 30°, 60°
    Point3D rotated_euler = rotation_matrix.multiply(original_point);
    std::cout << "After Euler rotation: ";
    rotated_euler.print();
    std::cout << std::endl;
    
    // Quaternion rotation
    Quaternion quat = Quaternion::fromAxisAngle(Point3D(0, 0, 1), PI/4); // 45° around Z
    Point3D rotated_quat = quat.rotatePoint(original_point);
    std::cout << "After quaternion rotation (45° around Z): ";
    rotated_quat.print();
    std::cout << std::endl;
    std::cout << std::endl;
    
    // ===== CONVEX HULL =====
    std::cout << "3. CONVEX HULL ALGORITHMS" << std::endl;
    
    std::vector<Point2D> points_2d = {{0, 3}, {1, 1}, {2, 2}, {4, 4}, {0, 0}, {1, 2}, {3, 1}, {3, 3}};
    auto hull_2d = convexHull2D(points_2d);
    
    std::cout << "2D Convex Hull points: ";
    for (const auto& p : hull_2d) {
        p.print();
        std::cout << " ";
    }
    std::cout << std::endl;
    
    GeometryGenerator gen(42);
    auto points_3d = gen.randomPointsInSphere(20, Point3D(0,0,0), 5.0);
    auto hull_3d = convexHull3D_Simple(points_3d);
    
    std::cout << "3D Convex Hull has " << hull_3d.size() << " extreme points" << std::endl;
    std::cout << std::endl;
    
    // ===== 3D AUDIO GEOMETRY =====
    std::cout << "4. 3D AUDIO GEOMETRY" << std::endl;
    
    AudioSource source(Point3D(2, 1, 0), 1.0);
    AudioListener listener(Point3D(0, 0, 0), Point3D(0, 0, -1), Point3D(0, 1, 0));
    
    auto audio_params = calculate3DAudio(source, listener);
    audio_params.print();
    std::cout << std::endl;
    
    // Room acoustics
    std::vector<Point3D> room_corners = {
        {-5, -3, -2}, {5, -3, -2}, {5, 3, -2}, {-5, 3, -2}, // Floor
        {-5, -3, 2}, {5, -3, 2}, {5, 3, 2}, {-5, 3, 2}    // Ceiling
    };
    
    auto reflections = calculateReflectionPoints(source, listener, room_corners);
    std::cout << "First-order reflections: " << reflections.size() << " points" << std::endl;
    std::cout << std::endl;
    
    // ===== WAVE INTERFERENCE =====
    std::cout << "5. WAVE INTERFERENCE PATTERNS" << std::endl;
    
    std::vector<Wave> waves = {
        Wave(Point3D(-1, 0, 0), 440.0, 1.0, 0.0),     // 440 Hz source on left
        Wave(Point3D(1, 0, 0), 440.0, 1.0, PI/2)      // 440 Hz source on right (90° phase shift)
    };
    
    auto interference = calculateInterferencePattern(waves, Point2D(-3, -2), Point2D(3, 2), 50, 30, 0.0, 0.0);
    
    std::cout << "Interference pattern calculated for 50x30 grid" << std::endl;
    std::cout << "Sample values from center row:" << std::endl;
    int center_row = interference.size() / 2;
    for (int j = 0; j < std::min(10, (int)interference[center_row].size()); j += 2) {
        std::cout << interference[center_row][j] << " ";
    }
    std::cout << std::endl << std::endl;
    
    // ===== ADVANCED COMPUTATIONS =====
    std::cout << "6. ADVANCED GEOMETRIC COMPUTATIONS" << std::endl;
    
    // Ray-triangle intersection
    Point3D ray_origin(0, 0, 5);
    Point3D ray_direction(0, 0, -1);
    Point3D v0(0, 0, 0), v1(1, 0, 0), v2(0, 1, 0);
    
    auto intersection = rayTriangleIntersect(ray_origin, ray_direction, v0, v1, v2);
    std::cout << "Ray-triangle intersection: ";
    intersection.print();
    std::cout << std::endl;
    
    // Closest pair of points
    auto random_points = gen.randomPointsInCircle(20, Point2D(0, 0), 10.0);
    auto closest_pair = closestPairOfPoints(random_points);
    std::cout << "Closest pair of points: ";
    closest_pair.print();
    std::cout << std::endl;
    
    // Bounding box
    auto sphere_points = gen.randomPointsInSphere(100, Point3D(0, 0, 0), 5.0);
    auto bbox = BoundingBox3D::fromPoints(sphere_points);
    std::cout << "Bounding box of 100 random points: ";
    bbox.print();
    std::cout << std::endl << std::endl;
    
    // ===== SYNTHESIZER APPLICATIONS =====
    std::cout << "7. SYNTHESIZER AND AUDIO APPLICATIONS" << std::endl;
    
    std::cout << "🎵 3D Audio Processing:" << std::endl;
    std::cout << "   • Spatial audio positioning using spherical coordinates" << std::endl;
    std::cout << "   • HRTF calculation based on azimuth/elevation angles" << std::endl;
    std::cout << "   • Distance-based attenuation and Doppler effects" << std::endl;
    std::cout << "   • Room impulse response from reflection geometry" << std::endl;
    std::cout << std::endl;
    
    std::cout << "🎛️ Physical Modeling:" << std::endl;
    std::cout << "   • String vibration using geometric wave equations" << std::endl;
    std::cout << "   • Drum membrane modes with circular geometry" << std::endl;
    std::cout << "   • Brass instrument bore calculations" << std::endl;
    std::cout << "   • Reed dynamics in woodwind instruments" << std::endl;
    std::cout << std::endl;
    
    std::cout << "📐 Speaker Array Design:" << std::endl;
    std::cout << "   • Optimal speaker placement using convex hull algorithms" << std::endl;
    std::cout << "   • Interference pattern optimization" << std::endl;
    std::cout << "   • Sweet spot calculation in listening rooms" << std::endl;
    std::cout << "   • Line array beam forming geometry" << std::endl;
    std::cout << std::endl;
    
    std::cout << "🔧 Circuit Layout Optimization:" << std::endl;
    std::cout << "   • Component placement using geometric algorithms" << std::endl;
    std::cout << "   • Trace routing with collision detection" << std::endl;
    std::cout << "   • EMI minimization through geometric separation" << std::endl;
    std::cout << "   • Heat distribution modeling in 3D" << std::endl;
    std::cout << std::endl;
    
    // ===== PERFORMANCE CHARACTERISTICS =====
    std::cout << "8. ALGORITHM PERFORMANCE" << std::endl;
    
    std::cout << "Algorithm                 | Time Complexity | Space | Applications" << std::endl;
    std::cout << "========================= | =============== | ===== | ==================" << std::endl;
    std::cout << "Convex Hull (Andrew)      | O(n log n)      | O(n)  | Speaker placement" << std::endl;
    std::cout << "Closest Pair (D&C)       | O(n log n)      | O(n)  | Collision detection" << std::endl;
    std::cout << "Ray-Triangle Intersect    | O(1)            | O(1)  | 3D audio occlusion" << std::endl;
    std::cout << "Great Circle Distance     | O(1)            | O(1)  | Global positioning" << std::endl;
    std::cout << "Spherical Trigonometry    | O(1)            | O(1)  | Satellite audio" << std::endl;
    std::cout << "Wave Interference         | O(n²)           | O(n²) | Room acoustics" << std::endl;
    std::cout << "3D Audio Calculation      | O(1)            | O(1)  | Real-time spatial" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== ADVANCED TRIGONOMETRY AND GEOMETRY MODULE COMPLETE ===" << std::endl;
    std::cout << "🎵 Ready for advanced 3D audio processing and geometric synthesis! 🎵" << std::endl;
    
    return 0;
}