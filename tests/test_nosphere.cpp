// Tests of the non-spherical Lagrangian particles (problemID 22, LagrangeNoSphereParticles): resistance of
// an ellipsoid, Jeffery orbits, the sphere limit (identical to problemID 2), settling of an inclined spheroid,
// rotation in a shear flow (uniform and refined mesh) and the contact of an ellipsoid with a wall.
// Lattice units equal physical units (latticeParams).

#include "TestFramework.h"
#include "TestCase.h"
#include <iostream>

using namespace std;

namespace {
	const double PI = 3.14159265358979323846;
	typedef LagrangeNoSphereParticles NSP;

	void report(const string& what, double value) {
		cout << "       " << what << " = " << value << endl;
	}

	// quaternion of the rotation by angle about the axis y
	array<double, 4> rotY(double angle) {
		return { cos(angle / 2), 0, sin(angle / 2), 0 };
	}

	vector<double> ellipsoid(double x, double y, double z, double a1, double a2, double a3, double rho, array<double, 4> q) {
		return { x, y, z, a1, a2, a3, rho, q[0], q[1], q[2], q[3] };
	}

	// elliptic integrals by quadrature (lambda = exp(s)), independent of the Carlson functions
	void integrals(const double a[3], double& chi0, double alpha[3]) {
		double abc = a[0] * a[1] * a[2];
		chi0 = 0;
		alpha[0] = alpha[1] = alpha[2] = 0;
		const int n = 40000;
		const double s0 = -40, s1 = 60, h = (s1 - s0) / n;
		for (int k = 0; k <= n; k++) {
			double l = exp(s0 + k * h);
			double D = sqrt((a[0] * a[0] + l) * (a[1] * a[1] + l) * (a[2] * a[2] + l));
			double wgt = (k == 0 || k == n ? 0.5 : 1.0) * h * l; // trapezoid, dlambda = lambda ds
			chi0 += wgt / D;
			for (int i = 0; i < 3; i++) {
				alpha[i] += wgt / ((a[i] * a[i] + l) * D);
			}
		}
		chi0 *= abc;
		for (int i = 0; i < 3; i++) {
			alpha[i] *= abc;
		}
	}

	// slip velocity of a settling ellipsoid in a quiescent fluid: K(|w|) w = (m - V) g, fixed point iteration
	array<double, 3> settlingSlip(const double a[3], const array<double, 4>& q, double rhoP, const double g[3], double nu) {
		double KT[3], KR[3], R[9];
		NSP::resistance(a, KT, KR);
		NSP::rotationMatrix(q.data(), R);
		double V = 4.0 / 3.0 * PI * a[0] * a[1] * a[2];
		double d = 2 * cbrt(a[0] * a[1] * a[2]);
		double F[3] = { (rhoP - 1) * V * g[0], (rhoP - 1) * V * g[1], (rhoP - 1) * V * g[2] };
		array<double, 3> w = { 0, 0, 0 };
		for (int it = 0; it < 200; it++) {
			double wm = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
			double Re = d * wm / nu;
			double corr = 1 + 0.15 * pow(Re, 0.687);
			// w = R diag(1 / k) R^T F
			double Fb[3];
			for (int i = 0; i < 3; i++) {
				Fb[i] = (R[i] * F[0] + R[3 + i] * F[1] + R[6 + i] * F[2]) / (nu * KT[i] * corr);
			}
			for (int r = 0; r < 3; r++) {
				w[r] = R[3 * r] * Fb[0] + R[3 * r + 1] * Fb[1] + R[3 * r + 2] * Fb[2];
			}
		}
		return w;
	}

	CaseSpec closedBox(const string& name, double nu, double g) {
		CaseSpec s;
		s.name = name;
		s.nBi = s.nBj = s.nBk = 2;
		s.ni = s.nj = s.nk = 8;
		s.params = latticeParams(22, 8, nu);
		s.flag = [](int, int, int) { return 1; };
		ostringstream os;
		os.precision(17);
		os << -g;
		s.params["gz"] = os.str();
		s.params["initializeWithHydroPress"] = "1";
		return s;
	}

	// Couette flow u = U z / 16 between a wall at z = 0 and the moving lid at z = 16, periodic in x and y
	CaseSpec couette(const string& name, double nu, double U, bool refined) {
		CaseSpec s;
		s.name = name;
		s.ni = s.nj = s.nk = 8;
		s.octree = octreeBlocks(2, 1, 2, refined ? 1 : 0, [](int l, int i, int, int) { return l == 0 && i == 0; });
		s.params = latticeParams(22, 8, nu);
		s.params["xPeriodicity"] = "1";
		s.params["yPeriodicity"] = "1";
		s.params["movingLidTop"] = "1";
		ostringstream os;
		os.precision(17);
		os << U;
		s.params["uwall"] = os.str();
		s.flagL = [](int, int, int, int) { return 1; };
		s.initL = [U](int l, int, int, int gk) {
			double z = (gk + 0.5) / (1 << l);
			return array<double, 4>{ 1, U * z / 16, 0, 0 };
		};
		return s;
	}
}

// Stokes resistance of ellipsoids: Carlson integrals against quadrature, closed forms for a sphere and for
// a prolate spheroid (Chwang & Wu 1975)
TEST(nosphere_resistance_of_ellipsoids) {
	double KT[3], KR[3];
	double sphere[3] = { 1.3, 1.3, 1.3 };
	NSP::resistance(sphere, KT, KR);
	for (int i = 0; i < 3; i++) {
		CHECK_NEAR(KT[i], 6 * PI * 1.3, 1e-12 * KT[i]);
		CHECK_NEAR(KR[i], 8 * PI * pow(1.3, 3), 1e-12 * KR[i]);
	}
	// prolate spheroid, the symmetry axis is the body axis 1
	double a = 3, b = 1;
	double pro[3] = { a, b, b };
	NSP::resistance(pro, KT, KR);
	double e = sqrt(1 - b * b / (a * a)), L = log((1 + e) / (1 - e));
	double axial = 16 * PI * a * pow(e, 3) / ((1 + e * e) * L - 2 * e);
	double broad = 32 * PI * a * pow(e, 3) / (2 * e + (3 * e * e - 1) * L);
	report("prolate 3:1, axial KT / exact", KT[0] / axial);
	report("prolate 3:1, broadside KT / exact", KT[1] / broad);
	CHECK_NEAR(KT[0], axial, 1e-10 * axial);
	CHECK_NEAR(KT[1], broad, 1e-10 * broad);
	CHECK_NEAR(KT[2], broad, 1e-10 * broad);
	CHECK(KT[0] < KT[1]); // moving along the axis is easier
	// general ellipsoid: the same formulas with the integrals evaluated by quadrature
	double ell[3] = { 1.0, 0.6, 0.35 };
	NSP::resistance(ell, KT, KR);
	double chi0, alpha[3];
	integrals(ell, chi0, alpha);
	double abc = ell[0] * ell[1] * ell[2];
	double worst = 0;
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3, k = (i + 2) % 3;
		double kt = 16 * PI * abc / (chi0 + ell[i] * ell[i] * alpha[i]);
		double kr = 16 * PI * abc * (ell[j] * ell[j] + ell[k] * ell[k]) / (3 * (ell[j] * ell[j] * alpha[j] + ell[k] * ell[k] * alpha[k]));
		worst = max(worst, max(fabs(KT[i] / kt - 1), fabs(KR[i] / kr - 1)));
	}
	report("triaxial ellipsoid, largest relative difference to quadrature", worst);
	CHECK(worst < 1e-8);
}

// Jeffery orbit of a prolate spheroid (aspect ratio 3) in the simple shear u = gamma z: the axis turns by
// 2 pi in T = 2 pi / gamma (r + 1/r); without inertia the integrator follows the analytic orbit
TEST(nosphere_jeffery_period) {
	const double r = 3, gamma = 1e-3, dt = 0.2;
	double a[3] = { r * 0.2, 0.2, 0.2 }, KT[3], KR[3];
	NSP::resistance(a, KT, KR);
	double I[3] = { 0, 0, 0 };
	double grad[9] = { 0, 0, gamma, 0, 0, 0, 0, 0, 0 }; // du/dz
	double q[4] = { 1, 0, 0, 0 }, w[3] = { 0, 0, 0 };
	double angle = 0, prev = 0, t = 0;
	double T = 2 * PI / gamma * (r + 1 / r);
	while (angle < 2 * PI && t < 2 * T) {
		NSP::rotationStep(q, w, a, KR, I, 0.1, grad, dt);
		t += dt;
		double R[9];
		NSP::rotationMatrix(q, R);
		double phi = atan2(-R[6], R[0]); // angle of the axis in the x-z plane, positive towards -z
		double d = phi - prev;
		if (d < -PI) d += 2 * PI;
		if (d > PI) d -= 2 * PI;
		angle += d;
		prev = phi;
	}
	report("period / Jeffery period", t / T);
	CHECK_NEAR(t, T, 0.002 * T);
	// the axis stays in the x-z plane (the orbit constant of an orbit in the shear plane)
	double R[9];
	NSP::rotationMatrix(q, R);
	CHECK(fabs(R[3]) < 1e-6);
}

// a sphere (a1 = a2 = a3) of problemID 22 behaves exactly as the particle of problemID 2
TEST(nosphere_sphere_limit_equals_problem2) {
	const double d = 0.8, rhoP = 20, g = 1e-4, nu = 0.02;
	CaseSpec s2 = closedBox("ns_sphere2", nu, g);
	s2.params["problemID"] = "2";
	s2.particles.push_back({ 8.3, 7.6, 12.4, d, rhoP });
	CaseSpec s22 = closedBox("ns_sphere22", nu, g);
	s22.particleRows.push_back(ellipsoid(8.3, 7.6, 12.4, d / 2, d / 2, d / 2, rhoP, { 1, 0, 0, 0 }));
	auto m2 = buildMesh(s2);
	auto m22 = buildMesh(s22);
	runSteps(*m2, 600);
	runSteps(*m22, 600);
	LagrangeParticles& P = *m2->particles;
	LagrangeParticles& Q = *m22->particles;
	double dz = fabs(P.partZ[0] - Q.partZ[0]), dv = fabs(P.veloZ[0] - Q.veloZ[0]);
	report("|z(22) - z(2)|", dz);
	report("|vz(22) - vz(2)| / |vz|", dv / fabs(P.veloZ[0]));
	CHECK(dz < 1e-9);
	CHECK(dv < 1e-8 * fabs(P.veloZ[0]));
	CHECK(fabs(P.partX[0] - Q.partX[0]) < 1e-9);
}

// a prolate spheroid inclined by 45 degrees glides sideways: the slip velocity is K^-1 (m - V) g in direction
// and magnitude, the orientation is kept (no torque on a settling ellipsoid in creeping flow)
TEST(nosphere_inclined_settling) {
	const double rhoP = 20, g = 1e-4, nu = 0.02;
	double a[3] = { 0.6, 0.3, 0.3 };
	array<double, 4> q = rotY(PI / 4);
	CaseSpec s = closedBox("ns_inclined", nu, g);
	s.particleRows.push_back(ellipsoid(8.3, 7.6, 12.4, a[0], a[1], a[2], rhoP, q));
	auto mesh = buildMesh(s);
	NSP& P = *(NSP*)mesh->particles;
	runSteps(*mesh, 1200);
	double gv[3] = { 0, 0, -g };
	array<double, 3> w = settlingSlip(a, q, rhoP, gv, nu);
	double slip[3] = { P.veloX[0] - P.veloFX[0], P.veloY[0] - P.veloFY[0], P.veloZ[0] - P.veloFZ[0] };
	double nw = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
	double ns = sqrt(slip[0] * slip[0] + slip[1] * slip[1] + slip[2] * slip[2]);
	double cosA = (slip[0] * w[0] + slip[1] * w[1] + slip[2] * w[2]) / (nw * ns);
	double angleDeg = acos(min(cosA, 1.0)) * 180 / PI;
	double R[9];
	NSP::rotationMatrix(&P.quat[0], R);
	double tilt = atan2(-R[6], R[0]) * 180 / PI;
	report("expected glide angle from the vertical [deg]", atan2(fabs(w[0]), fabs(w[2])) * 180 / PI);
	report("slip / expected (magnitude)", ns / nw);
	report("angle between slip and expected [deg]", angleDeg);
	report("tilt of the axis after the run [deg] (initially 45)", tilt);
	// the particle moves with slip + fluid velocity; the fluid at a point particle is dragged down by its own
	// (vertical) force, so the path is steeper than the slip direction (self-induced velocity, see CHANGES.md)
	double drift = (P.partX[0] - 8.3) / (12.4 - P.partZ[0]);
	double expected = (w[0] + P.veloFX[0]) / -(w[2] + P.veloFZ[0]);
	report("self-induced fluid velocity at the particle / vertical slip", P.veloFZ[0] / w[2]);
	report("lateral drift / fall", drift);
	report("lateral drift / fall, expected (slip + fluid velocity)", expected);
	CHECK_NEAR(ns, nw, 0.02 * nw);
	CHECK(angleDeg < 1.0);
	CHECK(fabs(tilt - 45) < 2.0);
	CHECK_NEAR(drift, expected, 0.05 * fabs(expected)); // it glided sideways as predicted
}

// a neutrally buoyant prolate spheroid (r = 3) aligned with the flow and a sphere in a Couette flow: the axis
// turns with gamma / (r^2 + 1), the sphere spins with gamma / 2; uniform mesh and the spheroid in the refined
// half (the sphere in the coarse half)
TEST(nosphere_rotation_in_shear_flow) {
	const double nu = 0.1, U = 0.02, gamma = U / 16, r = 3;
	for (bool refined : { false, true }) {
		CaseSpec s = couette(refined ? "ns_couette_ref" : "ns_couette", nu, U, refined);
		s.particleRows.push_back(ellipsoid(4.3, 4.2, 8.0, 0.2 * r, 0.2, 0.2, 1.0, { 1, 0, 0, 0 }));
		s.particleRows.push_back(ellipsoid(12.3, 4.2, 8.0, 0.3, 0.3, 0.3, 1.0, { 1, 0, 0, 0 }));
		auto mesh = buildMesh(s);
		NSP& P = *(NSP*)mesh->particles;
		runSteps(*mesh, 40);
		double w0[3], w1[3];
		P.omegaLab(0, w0);
		P.omegaLab(1, w1);
		string tag = refined ? "refined: " : "uniform: ";
		// Jeffery angular velocity with the interpolated gradient: checks the particle model independently
		// of the flow error; omega_y = 0.1 du/dz - 0.9 dw/dx amplifies the error of dw/dx for r = 3
		double wJ[3], wJlab[3], Rq[9];
		NSP::jefferyOmega(&P.quat[0], &P.semiAxes[0], &P.gradF[0], wJ);
		NSP::rotationMatrix(&P.quat[0], Rq);
		for (int i = 0; i < 3; i++) {
			wJlab[i] = Rq[3 * i] * wJ[0] + Rq[3 * i + 1] * wJ[1] + Rq[3 * i + 2] * wJ[2];
		}
		double others = 0;
		for (int c = 0; c < 9; c++) {
			others = max(others, c == 2 ? 0.0 : fabs(P.gradF[c]));
		}
		report(tag + "du/dz at the spheroid / gamma", P.gradF[2] / gamma);
		report(tag + "largest other gradient component / gamma", others / gamma);
		report(tag + "omega_y of the spheroid / (gamma / (r^2 + 1))", w0[1] / (gamma / (r * r + 1)));
		report(tag + "omega_y of the spheroid / Jeffery with the interpolated gradient", w0[1] / wJlab[1]);
		report(tag + "omega_y of the sphere / (gamma / 2)", w1[1] / (gamma / 2));
		CHECK_NEAR(P.gradF[2], gamma, 0.02 * gamma);
		CHECK(others < 0.02 * gamma);
		CHECK_NEAR(w0[1], wJlab[1], 1e-3 * fabs(wJlab[1]));
		CHECK_NEAR(w1[1], gamma / 2, 0.02 * gamma / 2);
		CHECK(fabs(w0[0]) + fabs(w0[2]) < 1e-3 * gamma);
		if (!refined) {
			CHECK_NEAR(w0[1], gamma / (r * r + 1), 0.01 * gamma / (r * r + 1));
		}
	}
}

// an ellipsoid lying on the floor rests at the height of its short semi-axis
TEST(nosphere_wall_contact) {
	const double rhoP = 20, g = 1e-4, nu = 0.02;
	CaseSpec s = closedBox("ns_wall", nu, g);
	s.particleRows.push_back(ellipsoid(8.3, 7.6, 2.0, 0.6, 0.35, 0.2, rhoP, { 1, 0, 0, 0 }));
	auto mesh = buildMesh(s);
	NSP& P = *(NSP*)mesh->particles;
	runSteps(*mesh, 2500);
	double R[9];
	NSP::rotationMatrix(&P.quat[0], R);
	report("z - a3 (short semi-axis)", P.partZ[0] - 0.2);
	report("z component of the long axis", R[6]);
	CHECK(P.partZ[0] > 0.2 - 0.02 && P.partZ[0] < 0.2 + 0.05);
	CHECK(fabs(R[6]) < 0.05);
}
