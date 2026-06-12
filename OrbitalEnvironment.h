#pragma once
#include <cmath>

namespace Orbital {
    struct Vec3 {
        double x, y, z;

        Vec3() : x(0), y(0), z(0) {
        }

        Vec3(double x, double y, double z) : x(x), y(y), z(z) {
        }

        Vec3 operator+(const Vec3 &o) const {
            return {x + o.x, y + o.y, z + o.z};
        }

        Vec3 operator-(const Vec3 &o) const {
            return {x - o.x, y - o.y, z - o.z};
        }

        Vec3 operator*(const double &o) const {
            return {x * o, y * o, z * o};
        }

        Vec3 operator/(const double &o) const {
            return {x / o, y / o, z / o};
        }

        [[nodiscard]] double dot(const Vec3 &o) const {
            return (x * o.x) + (y * o.y) + (z * o.z);
        }

        [[nodiscard]] Vec3 cross(const Vec3 &o) const {
            return {(y * o.z) - (z * o.y), (z * o.x) - (x * o.z), (x * o.y) - (y * o.x)};
        }

        [[nodiscard]] double norm() const {
            return sqrt(dot(*this));
        }

        [[nodiscard]] Vec3 normalize() const {
            return *this / norm();
        }
    };
    struct Matrix3 {
        double m[9]{};

        Matrix3() = default;

        Matrix3(double m00, double m01, double m02,
                double m10, double m11, double m12,
                double m20, double m21, double m22) {
            m[0] = m00;
            m[1] = m01;
            m[2] = m02;
            m[3] = m10;
            m[4] = m11;
            m[5] = m12;
            m[6] = m20;
            m[7] = m21;
            m[8] = m22;
        }

        Vec3 operator*(const Vec3 &o) const {
            return Vec3{
                (m[0] * o.x) + (m[1] * o.y) + (m[2] * o.z),
                (m[3] * o.x) + (m[4] * o.y) + (m[5] * o.z),
                (m[6] * o.x) + (m[7] * o.y) + (m[8] * o.z)
            };
        }

        Matrix3 operator*(const Matrix3 &o) const {
            return {
                (m[0] * o.m[0]) + (m[1] * o.m[3]) + (m[2] * o.m[6]),
                (m[0] * o.m[1]) + (m[1] * o.m[4]) + (m[2] * o.m[7]),
                (m[0] * o.m[2]) + (m[1] * o.m[5]) + (m[2] * o.m[8]),
                (m[3] * o.m[0]) + (m[4] * o.m[3]) + (m[5] * o.m[6]),
                (m[3] * o.m[1]) + (m[4] * o.m[4]) + (m[5] * o.m[7]),
                (m[3] * o.m[2]) + (m[4] * o.m[5]) + (m[5] * o.m[8]),
                (m[6] * o.m[0]) + (m[7] * o.m[3]) + (m[8] * o.m[6]),
                (m[6] * o.m[1]) + (m[7] * o.m[4]) + (m[8] * o.m[7]),
                (m[6] * o.m[2]) + (m[7] * o.m[5]) + (m[8] * o.m[8])
            };
        }

        [[nodiscard]] Matrix3 transpose() const {
            return {
                m[0], m[3], m[6],
                m[1], m[4], m[7],
                m[2], m[5], m[8]
            };
        }

        static Matrix3 skew(const Vec3 &o) {
            return {
                0, -o.z, o.y,
                o.z, 0, -o.x,
                -o.y, o.x, 0
            };
        }
    };
    struct Quat {
        double w, x, y, z;

        Quat() : w(0), x(0), y(0), z(0) {
        }

        Quat(double w, double x, double y, double z) : w(w), x(x), y(y), z(z) {
        }

        Quat operator*(const double s) const {
            return {w * s, x * s, y * s, z * s};
        }

        Quat operator*(const Quat &o) const {
            return {
                (w * o.w) - (x * o.x) - (y * o.y) - (z * o.z),
                (w * o.x) + (x * o.w) + (y * o.z) - (z * o.y),
                (w * o.y) - (x * o.z) + (y * o.w) + (z * o.x),
                (w * o.z) + (x * o.y) - (y * o.x) + (z * o.w)
            };
        }

        [[nodiscard]] Quat conjugate() const {
            return {w, -x, -y, -z};
        }

        [[nodiscard]] Vec3 rotate(const Vec3 &o) const {
            Quat p = {0, o.x, o.y, o.z};
            Quat q_i = (*this) * p;
            Quat q = q_i * conjugate();

            return {q.x, q.y, q.z};
        }

        [[nodiscard]] Quat derivative(const Vec3 &omega) const {
            Quat p = {0, omega.x, omega.y, omega.z};
            return ((*this) * p) * 0.5;
        }

        [[nodiscard]] Quat normalize() const {
            double norm = sqrt((w * w) + (x * x) + (y * y) + (z * z));
            return (*this) * (1 / norm);
        }

        [[nodiscard]] Matrix3 toMatrix() const {
            return {
                1 - 2 * ((y * y) + (z * z)), 2 * ((x * y) - (w * z)), 2 * ((x * z) + (w * y)),
                2 * ((x * y) + (w * z)), 1 - 2 * ((x * x) + (z * z)), 2 * ((y * z) - (w * x)),
                2 * ((x * z) - (w * y)), 2 * ((y * z) + (w * x)), 1 - 2 * ((x * x) + (y * y))
            };
        }
    };

    struct OrbitalElements {
        double a;       // semi-major axis [m]
        double e;       // eccentricity [dimensionless]
        double i;       // inclination [rad]
        double raan;    // right ascension of ascending node [rad]
        double argp;    // argument of perigee [rad]
        double nu;      // true anomaly [rad] -- changes with time
    };

    struct CartesianState {
        Vec3 position;  // [m] in ECI
        Vec3 velocity;  // [m/s] in ECI
    };

    enum class EclipseState { SUNLIT, PENUMBRA, UMBRA };
    struct EnvironmentState {
        CartesianState orbit;
        Vec3 B_field;
        Vec3 sun_vector;
        EclipseState eclipse{};
        double jd{};
    };

    double solveKepler(double M, double e, double tol = 1e-10);

    CartesianState keplerPropagate(const OrbitalElements& oe, double dt);

    Vec3 dipoleMagneticField(const Vec3& r_ECI, double t_seconds);

    Vec3 sunVectorEci(double jd);

    // shadow models
    bool isEclipse(const Vec3& rECI, const Vec3& sunECI);   // cylindrical
    EclipseState eclipseState (const Vec3& rECI, const Vec3& rSunECI);  // conical
};
