#include "OrbitalEnvironment.h"
#include <cmath>


namespace {
    using namespace Orbital;

    constexpr double MU = 3.986004418e14; // std grav param of Earth
    constexpr double R_EARTH = 6371000.0;
    constexpr double R_SUN = 696999999.0;
    constexpr double AU = 1.496e11;

    constexpr double toRad(double deg) {
        return (deg * M_PI) / 180;
    }

    Matrix3 R1(double angle) {
        return {
            1, 0, 0,
            0, cos(angle), -sin(angle),
            0, sin(angle), cos(angle)
        };
    }

    Matrix3 R3(double angle) {
        return {
            cos(angle), -sin(angle), 0,
            sin(angle), cos(angle), 0,
            0, 0, 1
        };
    }

    double computeGST(double t_seconds) {
        constexpr double GST_J2000 = toRad(280.46061837); // radians
        constexpr double OMEGA_EARTH = 7.2921150e-5;

        return GST_J2000 + OMEGA_EARTH * t_seconds;
    }

    Matrix3 eciToECEF(double gst) {
        return R3(gst);
    }

    Matrix3 ecefToDipole() {
        constexpr double angle = toRad(11.5);
        return R1(angle);
    }

    [[maybe_unused]] double julianDate(int year, int month, int day, double hour, double min, double sec) {
        return 367 * year
               - (int) (7 * (year + (int) ((month + 9) / 12)) / 4)
               + (int) (275 * month / 9)
               + day
               + 1721013.5
               + (hour + min / 60.0 + sec / 3600.0) / 24.0;
    }

}

namespace Orbital {
    double solveKepler(double M, double e, double tol) {
        double E = M;

        for (int i = 0; i < 100; i++) {
            double f = E - e * sin(E) - M;
            if (std::abs(f) < tol) break;
            double df = 1 - e * cos(E);
            E -= f / df;
        }

        return E;
    }

    CartesianState keplerPropagate(const OrbitalElements &oe, double dt) {
        // solve for M
        double E0 = 2 * atan2(sqrt(1 - oe.e) * sin(oe.nu / 2), sqrt(1 + oe.e) * cos(oe.nu / 2)); // eccentric anomaly

        // compute mean anomaly
        double M0 = E0 - oe.e * sin(E0); // mean anomaly

        double n = sqrt(MU / pow(oe.a, 3));
        double M = M0 + n * dt;

        double E = solveKepler(M, oe.e);

        double nu = 2 * atan2(sqrt(1 + oe.e) * sin(E / 2), sqrt(1 - oe.e) * cos(E / 2));

        double p = oe.a * (1 - pow(oe.e, 2));
        double r = p / (1 + oe.e * cos(nu));

        Vec3 rPQW = Vec3{cos(nu), sin(nu), 0} * r;
        Vec3 vPQW = Vec3{-sin(nu), oe.e + cos(nu), 0} * sqrt(MU / p);

        Matrix3 R = R3(oe.raan) * R1(oe.i) * R3(oe.argp);

        Vec3 rECI = R * rPQW;
        Vec3 vECI = R * vPQW;

        return CartesianState{rECI, vECI};
    }

    Vec3 dipoleMagneticField(const Vec3 &r_ECI, double t_seconds) {
        constexpr double MU_M = 8.0e22;

        double gst = computeGST(t_seconds);
        Vec3 r_ECEF = eciToECEF(gst) * r_ECI;
        Vec3 r_dipole = ecefToDipole() * r_ECEF;

        const double x = r_dipole.x;
        const double y = r_dipole.y;
        const double z = r_dipole.z;

        const double r = r_dipole.norm();
        double lambda_m = asin(z / r);
        double rho = sqrt(pow(x, 2) + pow(y, 2));

        double B_r = -2 * (MU_M / pow(r, 3)) * sin(lambda_m);
        double B_lambda = -(MU_M / pow(r, 3)) * cos(lambda_m);
        // double B_phi = 0;

        double B_x = (B_r * cos(lambda_m) - B_lambda * sin(lambda_m)) * (x / rho);
        double B_y = (B_r * cos(lambda_m) - B_lambda * sin(lambda_m)) * (y / rho);
        double B_z = B_r * sin(lambda_m) + B_lambda * cos(lambda_m);

        Vec3 B_dipole = {B_x, B_y, B_z};
        Vec3 B_ECEF = ecefToDipole().transpose() * B_dipole;
        Vec3 B_ECI = eciToECEF(gst).transpose() * B_ECEF;

        return B_ECI;
    }

    Vec3 sunVectorEci(double jd) {
        double T_uti = (jd - 2451545.0) / 36525;
        double lambda_sun = toRad(280.460 + 36000.771 * T_uti);
        double M_sun = toRad(357.529 + 35999.050 * T_uti);
        double lambda_ecl = lambda_sun + toRad(1.915 * sin(M_sun) + 0.02 * sin(2 * M_sun));
        double epsilon = toRad(23.439 - 0.012 * T_uti);

        Vec3 sunVec = { cos(lambda_ecl), sin(lambda_ecl) * cos(epsilon), sin(lambda_ecl) * sin(epsilon) };
        return sunVec * AU;
    }

    // shadow models
    bool isEclipse(const Vec3 &rECI, const Vec3 &sunECI) {
        Vec3 sunUnit = sunECI.normalize();
        double proj = rECI.dot(sunUnit);
        if (proj > 0) return false;
        Vec3 perpVec = rECI - sunUnit * proj;
        double perp = perpVec.norm();
        return perp < R_EARTH;
    }

    EclipseState eclipseState(const Vec3 &rECI, const Vec3 &rSunECI) {
        double rSunDist = rSunECI.norm();
        Vec3 sunUnit = (rSunECI * -1).normalize();

        double sat_proj = rECI.dot(sunUnit);

        if (sat_proj < 0) return EclipseState::SUNLIT;

        Vec3 perp_vec = rECI - sunUnit * sat_proj;
        double perp = perp_vec.norm();

        double f_umbra = asin((R_SUN - R_EARTH) / rSunDist);
        double f_penumbra = asin((R_SUN + R_EARTH) / rSunDist);

        double r_umbra = R_EARTH - sat_proj * tan(f_umbra);
        double r_penumbra = R_EARTH + sat_proj * tan(f_penumbra);

        if (perp <= r_umbra) return EclipseState::UMBRA;
        if (perp <= r_penumbra) return EclipseState::PENUMBRA;
        return EclipseState::SUNLIT;
    }
}
