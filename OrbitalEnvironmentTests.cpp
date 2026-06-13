#include "OrbitalEnvironment.h"
#include <iostream>
#include <cmath>

using namespace Orbital;

static int passed = 0;
static int failed = 0;

void check(bool condition, const char* name) {
    if (condition) {
        std::cout << "  PASS  " << name << "\n";
        ++passed;
    } else {
        std::cout << "  FAIL  " << name << "\n";
        ++failed;
    }
}

bool near(double a, double b, double tol = 1e-6) {
    return std::abs(a - b) < tol;
}

bool nearVec(const Vec3& a, const Vec3& b, double tol = 1e-3) {
    return (a - b).norm() < tol;
}

void testSolveKepler() {
    std::cout << "\n-- solveKepler --\n";

    check(near(solveKepler(1.0, 0.0), 1.0),
          "circular orbit E == M when e=0");
    check(near(solveKepler(0.0, 0.5), 0.0),
          "M=0 gives E=0");
    check(near(solveKepler(M_PI, 0.3), M_PI, 1e-9),
          "M=PI gives E=PI");

    double M = 1.2, e = 0.5;
    double E = solveKepler(M, e);
    check(std::abs(E - e * std::sin(E) - M) < 1e-10,
          "residual below tolerance for e=0.5");
}

void testKeplerPropagate() {
    std::cout << "\n-- keplerPropagate --\n";

    constexpr double MU = 3.986004418e14;

    OrbitalElements circ;
    circ.a = 6771000.0; circ.e = 0.0;
    circ.i = circ.raan = circ.argp = circ.nu = 0.0;
    double T = 2 * M_PI * std::sqrt(std::pow(circ.a, 3) / MU);

    check(near(keplerPropagate(circ, 0).position.norm(), circ.a, 1.0),
          "circular: radius at t=0 equals a");
    check(near(keplerPropagate(circ, T/4).position.norm(), circ.a, 1.0),
          "circular: radius at T/4 equals a");

    OrbitalElements oe;
    oe.a = 7000000.0; oe.e = 0.01;
    oe.i = 0.5; oe.raan = 0.3; oe.argp = 0.8; oe.nu = 1.0;
    double Tp = 2 * M_PI * std::sqrt(std::pow(oe.a, 3) / MU);
    auto s0 = keplerPropagate(oe, 0);
    auto sT = keplerPropagate(oe, Tp);

    check(nearVec(s0.position, sT.position, 10.0),
          "position returns to start after one period");
    check(nearVec(s0.velocity, sT.velocity, 0.01),
          "velocity returns to start after one period");

    double v_expected = std::sqrt(MU / circ.a);
    check(near(keplerPropagate(circ, 0).velocity.norm(), v_expected, 0.01),
          "circular speed matches sqrt(mu/a)");
}

void testSunVectorEci() {
    std::cout << "\n-- sunVectorEci --\n";

    constexpr double AU = 1.496e11;

    Vec3 sv = sunVectorEci(2451545.0);
    check(near(sv.norm(), AU, AU * 0.02),
          "sun vector magnitude within 2% of 1 AU");

    Vec3 sv2 = sunVectorEci(2451545.0 + 182.5);
    check(sv.normalize().dot(sv2.normalize()) < -0.8,
          "sun vectors 6 months apart are roughly opposite");

    check(std::abs(sv.normalize().z) < 0.5,
          "sun vector z component small near J2000");
}

void testDipoleMagneticField() {
    std::cout << "\n-- dipoleMagneticField --\n";

    Vec3 r = {6771000.0, 0.0, 0.0};
    double mag = dipoleMagneticField(r, 0.0).norm();
    check(mag > 1e-6 && mag < 1e-3,
          "field magnitude at LEO is plausible");

    double B_low  = dipoleMagneticField({6571000, 0, 0}, 0).norm();
    double B_high = dipoleMagneticField({7371000, 0, 0}, 0).norm();
    check(B_low > B_high, "field stronger at lower altitude");

    Vec3 r1 = {6771000, 0, 0};
    Vec3 r2 = r1 * 2.0;
    double ratio = dipoleMagneticField(r1, 0).norm() /
                   dipoleMagneticField(r2, 0).norm();
    check(near(ratio, 8.0, 1.5), "field scales as 1/r^3");

    double B_pole = dipoleMagneticField({0, 0, 6771000}, 0).norm();
    double B_eq   = dipoleMagneticField({6771000, 0, 0}, 0).norm();
    check(B_pole > B_eq, "field stronger at pole than equator");
}

void testIsEclipse() {
    std::cout << "\n-- isEclipse --\n";

    constexpr double AU      = 1.496e11;
    constexpr double R_EARTH = 6371000.0;
    Vec3 sun = {AU, 0, 0};

    check(isEclipse({-7000000, 0, 0}, sun),
          "satellite directly behind Earth is eclipsed");
    check(!isEclipse({7000000, 0, 0}, sun),
          "satellite on day side is not eclipsed");
    check(!isEclipse({-7000000, 10000000, 0}, sun),
          "satellite far from shadow axis is not eclipsed");
    check(isEclipse({-7000000, R_EARTH * 0.99, 0}, sun),
          "satellite just inside cylinder edge is eclipsed");
    check(!isEclipse({-7000000, R_EARTH * 1.01, 0}, sun),
          "satellite just outside cylinder edge is not eclipsed");
}

void testEclipseState() {
    std::cout << "\n-- eclipseState --\n";

    constexpr double AU = 1.496e11;
    Vec3 sun = {AU, 0, 0};

    check(eclipseState({7000000, 0, 0}, sun) == EclipseState::SUNLIT,
          "day side is SUNLIT");
    check(eclipseState({-7000000, 0, 0}, sun) == EclipseState::UMBRA,
          "directly behind Earth is UMBRA");
    check(eclipseState({-7000000, 20000000, 0}, sun) == EclipseState::SUNLIT,
          "far outside shadow is SUNLIT");

    Vec3 rDay   = {7000000, 0, 0};
    Vec3 rNight = {-7000000, 0, 0};
    check(!isEclipse(rDay, sun) &&
           eclipseState(rDay, sun) == EclipseState::SUNLIT,
          "both models agree: day side");
    check(isEclipse(rNight, sun) &&
           eclipseState(rNight, sun) != EclipseState::SUNLIT,
          "both models agree: night side in shadow");
}

int main() {
    std::cout << "=== OrbitalEnvironment Test Suite ===\n";

    testSolveKepler();
    testKeplerPropagate();
    testSunVectorEci();
    testDipoleMagneticField();
    testIsEclipse();
    testEclipseState();

    std::cout << "\n=====================================\n";
    std::cout << "Passed: " << passed << "\n";
    std::cout << "Failed: " << failed << "\n";
    std::cout << "=====================================\n";

    return failed > 0 ? 1 : 0;
}