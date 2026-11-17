// Tests of the Eulerian non-spherical particles (problemID 33, BlockFreeSurfaceEulerNoSphereParticlesFVM):
// orientation equations and closures, the sphere limit (problemID 3), settling with the averaged mobility,
// Jeffery orbits and rotary diffusion of the orientation tensor in a Couette flow, transport of the orientation
// with the concentration (uniform flow, free surface, refined mesh).
// Lattice units equal physical units (latticeParams).

#include "TestFramework.h"
#include "TestCase.h"
#include <iostream>
#include <random>

using namespace std;

namespace {
	const double PI = 3.14159265358979323846;
	typedef BlockFreeSurfaceEulerNoSphereParticlesFVM ENS;
	typedef BlockFreeSurfaceEulerParticlesFVM FVM;
	typedef array<double, 6> Sym; // xx, yy, zz, xy, xz, yz

	void report(const string& what, double value) {
		cout << "       " << what << " = " << value << endl;
	}

	string num(double v) {
		ostringstream os;
		os.precision(17);
		os << v;
		return os.str();
	}

	Sym aligned(double px, double py, double pz) {
		double n = sqrt(px * px + py * py + pz * pz);
		px /= n;
		py /= n;
		pz /= n;
		return { px * px, py * py, pz * pz, px * py, px * pz, py * pz };
	}

	const Sym ISO = { 1.0 / 3, 1.0 / 3, 1.0 / 3, 0, 0, 0 };

	void setParticles(CaseSpec& s, double d, double rhoP, double r) {
		s.params["particleDiameter"] = num(d);
		s.params["particleDensity"] = num(rhoP);
		s.params["particleAspectRatio"] = num(r);
	}

	// closed box of nBi x nBj x nBk blocks of 8^3 cells full of liquid
	CaseSpec box(const string& name, int problemID, int nBi, int nBj, int nBk, double nu) {
		CaseSpec s;
		s.name = name;
		s.nBi = nBi;
		s.nBj = nBj;
		s.nBk = nBk;
		s.params = latticeParams(problemID, 8, nu);
		s.flag = [](int, int, int) { return 1; };
		return s;
	}

	// Couette flow u = U z / 16 between a wall at z = 0 and the moving lid at z = 16, periodic in x and y
	CaseSpec couette(const string& name, double nu, double U) {
		CaseSpec s;
		s.name = name;
		s.octree = octreeBlocks(2, 1, 2, 0, [](int, int, int, int) { return false; });
		s.params = latticeParams(33, 8, nu);
		s.params["xPeriodicity"] = "1";
		s.params["yPeriodicity"] = "1";
		s.params["movingLidTop"] = "1";
		s.params["uwall"] = num(U);
		s.flagL = [](int, int, int, int) { return 1; };
		s.initL = [U](int l, int, int, int gk) {
			double z = (gk + 0.5) / (1 << l);
			return array<double, 4>{ 1, U * z / 16, 0, 0 };
		};
		return s;
	}

	double maxDiff(const double* a, const Sym& b) {
		double m = 0;
		for (int k = 0; k < 6; k++) {
			m = max(m, fabs(a[k] - b[k]));
		}
		return m;
	}

	// smallest eigenvalue and |trace - 1| of the orientation tensors of the liquid cells
	array<double, 2> orientationValidity(Mesh& mesh) {
		array<double, 2> r = { 1, 0 };
		for (int b = 0; b < mesh.nBlock; b++) {
			ENS* bl = (ENS*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				if (bl->boundary[ind] || (bl->flagMap[ind] != 1 && bl->flagMap[ind] != 3)) {
					continue;
				}
				const double* A = &bl->orient[6 * ind];
				double val[3], vec[9];
				ENS::eigenSym(A, val, vec);
				r[0] = min(r[0], val[2]);
				r[1] = max(r[1], fabs(A[0] + A[1] + A[2] - 1));
			}
		}
		return r;
	}

	// sum of the concentration (all = false: liquid cells only) and its centroid
	array<double, 4> particleMoments(Mesh& mesh, bool all = true) {
		array<double, 4> m = { 0, 0, 0, 0 };
		for (int b = 0; b < mesh.nBlock; b++) {
			FVM* bl = (FVM*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				if (bl->boundary[ind] || (!all && bl->flagMap[ind] != 1 && bl->flagMap[ind] != 3)) {
					continue;
				}
				double c = bl->concentration[ind];
				int gi, gj, gk; // cells of the block's level
				bl->globalCoords(ind, gi, gj, gk);
				m[0] += c;
				m[1] += c * (gi + 0.5);
				m[2] += c * (gj + 0.5);
				m[3] += c * (gk + 0.5);
			}
		}
		for (int d = 1; d < 4; d++) {
			m[d] /= m[0];
		}
		return m;
	}

	// semi-axes of the spheroid of aspect ratio r with the volume of the sphere of diameter d
	void spheroid(double d, double r, double a[3]) {
		double b = d / 2 / cbrt(r);
		a[0] = r * b;
		a[1] = a[2] = b;
	}

	// settling slip velocity of a spheroid with the axis p: w = M F / g(Re), M = I / KTperp + (1 / KTpar - 1 / KTperp) p p
	array<double, 3> slipVelocity(double d, double r, double rhoP, const Sym& A, const double g[3], double nu) {
		double a[3], KT[3], KR[3];
		spheroid(d, r, a);
		LagrangeNoSphereParticles::resistance(a, KT, KR);
		double V = PI * d * d * d / 6;
		double F[3] = { (rhoP - 1) * V * g[0], (rhoP - 1) * V * g[1], (rhoP - 1) * V * g[2] };
		double M[9] = { A[0], A[3], A[4], A[3], A[1], A[5], A[4], A[5], A[2] };
		for (int k = 0; k < 9; k++) {
			M[k] *= 1 / KT[0] - 1 / KT[1];
		}
		M[0] += 1 / KT[1];
		M[4] += 1 / KT[1];
		M[8] += 1 / KT[1];
		array<double, 3> w = { 0, 0, 0 };
		for (int it = 0; it < 200; it++) {
			double Re = d * sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]) / nu;
			double corr = 1 + 0.15 * pow(Re, 0.687);
			for (int i = 0; i < 3; i++) {
				w[i] = (M[3 * i] * F[0] + M[3 * i + 1] * F[1] + M[3 * i + 2] * F[2]) / (nu * corr);
			}
		}
		return w;
	}
}

// the equations of one cell: drag tensor, closures, Jeffery rate, eigen decomposition, projection
TEST(euler_nosphere_orientation_equations) {
	const double d = 1.0, r = 3;
	CaseSpec s = box("ens_unit", 33, 1, 1, 1, 0.02);
	setParticles(s, d, 2.0, r);
	auto mesh = buildMesh(s);
	ENS& B = *(ENS*)mesh->blocks[mesh->blockId[0]];

	// resistances of the prolate spheroid against Chwang & Wu (1975)
	double a[3];
	spheroid(d, r, a);
	double e = sqrt(1 - a[1] * a[1] / (a[0] * a[0])), L = log((1 + e) / (1 - e));
	double axial = 16 * PI * a[0] * pow(e, 3) / ((1 + e * e) * L - 2 * e);
	double broad = 32 * PI * a[0] * pow(e, 3) / (2 * e + (3 * e * e - 1) * L);
	CHECK_NEAR(B.KTpar, axial, 1e-12 * axial);
	CHECK_NEAR(B.KTperp, broad, 1e-12 * broad);
	// drag tensor: aligned particles have the resistances of one particle, random ones the harmonic mean
	double K[9];
	B.dragTensor(aligned(1, 0, 0).data(), K);
	CHECK_NEAR(K[0], axial, 1e-12 * axial);
	CHECK_NEAR(K[4], broad, 1e-12 * broad);
	CHECK(fabs(K[1]) + fabs(K[2]) + fabs(K[5]) < 1e-12 * axial);
	B.dragTensor(ISO.data(), K);
	double kIso = 3 / (1 / axial + 2 / broad);
	CHECK_NEAR(K[0], kIso, 1e-12 * kIso);
	CHECK_NEAR(K[8], kIso, 1e-12 * kIso);
	report("K random / K across (r = 3)", kIso / broad);

	mt19937 gen(7);
	uniform_real_distribution<double> U(-1, 1);
	double worstJeffery = 0, worstIso = 0, worstTrace = 0, worstEigen = 0;
	for (int t = 0; t < 50; t++) {
		double G[9];
		for (int k = 0; k < 9; k++) {
			G[k] = 1e-3 * U(gen);
		}
		double tr = (G[0] + G[4] + G[8]) / 3; // incompressible flow
		G[0] -= tr;
		G[4] -= tr;
		G[8] -= tr;
		double D[9], W[9];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				D[3 * i + j] = (G[3 * i + j] + G[3 * j + i]) / 2;
				W[3 * i + j] = (G[3 * i + j] - G[3 * j + i]) / 2;
			}
		}
		// aligned particles: the closure is exact, dA/dt = p' p + p p' with the Jeffery rate of p
		double p[3] = { U(gen), U(gen), U(gen) };
		double n = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
		for (int i = 0; i < 3; i++) {
			p[i] /= n;
		}
		double Dp[3], Wp[3], pDp = 0;
		for (int i = 0; i < 3; i++) {
			Dp[i] = D[3 * i] * p[0] + D[3 * i + 1] * p[1] + D[3 * i + 2] * p[2];
			Wp[i] = W[3 * i] * p[0] + W[3 * i + 1] * p[1] + W[3 * i + 2] * p[2];
			pDp += p[i] * Dp[i];
		}
		double dp[3];
		for (int i = 0; i < 3; i++) {
			dp[i] = Wp[i] + B.lambdaJ * (Dp[i] - pDp * p[i]);
		}
		Sym A = aligned(p[0], p[1], p[2]);
		double dA[6];
		B.orientationRate(A.data(), G, dA);
		const int I[6] = { 0, 1, 2, 0, 0, 1 }, J[6] = { 0, 1, 2, 1, 2, 2 };
		for (int k = 0; k < 6; k++) {
			double ex = dp[I[k]] * p[J[k]] + p[I[k]] * dp[J[k]];
			worstJeffery = max(worstJeffery, fabs(dA[k] - ex));
		}
		// random orientation: the hybrid closure is the linear one, exact here, dA/dt = 2 lambda D / 5
		B.orientationRate(ISO.data(), G, dA);
		for (int k = 0; k < 6; k++) {
			worstIso = max(worstIso, fabs(dA[k] - 2 * B.lambdaJ * D[3 * I[k] + J[k]] / 5));
		}
		// any orientation: the trace is kept
		Sym R = { fabs(U(gen)), fabs(U(gen)), fabs(U(gen)), 0.1 * U(gen), 0.1 * U(gen), 0.1 * U(gen) };
		ENS::projectOrientation(R.data());
		B.orientationRate(R.data(), G, dA);
		worstTrace = max(worstTrace, fabs(dA[0] + dA[1] + dA[2]));
		// eigen decomposition reproduces the tensor
		double val[3], vec[9];
		ENS::eigenSym(R.data(), val, vec);
		CHECK(val[0] >= val[1] && val[1] >= val[2]);
		for (int k = 0; k < 6; k++) {
			double m = 0;
			for (int c = 0; c < 3; c++) {
				m += vec[3 * I[k] + c] * val[c] * vec[3 * J[k] + c];
			}
			worstEigen = max(worstEigen, fabs(m - R[k]));
		}
	}
	report("aligned particles: |dA/dt - Jeffery| / |grad u|", worstJeffery / 1e-3);
	report("random orientation: |dA/dt - 2 lambda D / 5| / |grad u|", worstIso / 1e-3);
	report("|trace dA/dt| / |grad u|", worstTrace / 1e-3);
	report("|V diag(val) V^T - A|", worstEigen);
	CHECK(worstJeffery < 1e-15);
	CHECK(worstIso < 1e-15);
	CHECK(worstTrace < 1e-15);
	CHECK(worstEigen < 1e-14);

	// projection: a valid tensor is kept, an indefinite one gets the nearest positive semidefinite tensor of trace 1
	Sym ok = aligned(1, 2, 3);
	Sym kept = ok;
	ENS::projectOrientation(kept.data());
	CHECK(maxDiff(kept.data(), ok) < 1e-15);
	Sym bad = { 1.2, -0.1, -0.1, 0, 0, 0 };
	ENS::projectOrientation(bad.data());
	CHECK(maxDiff(bad.data(), aligned(1, 0, 0)) < 1e-15);
}

// r = 1: the same result as problemID 3 (to round-off, the drag is computed as a tensor)
TEST(euler_nosphere_sphere_limit) {
	const double d = 1.0, rhoP = 20, g = 1e-4, nu = 0.02;
	unique_ptr<Mesh> m[2];
	for (int c = 0; c < 2; c++) {
		CaseSpec s = box(c ? "ens_sphere33" : "ens_sphere3", c ? 33 : 3, 1, 1, 4, nu);
		s.params["gz"] = num(-g);
		s.params["gx"] = num(0.3 * g);
		s.params["initializeWithHydroPress"] = "1";
		setParticles(s, d, rhoP, 1.0);
		s.concentration = [](int i, int, int k) { return (k >= 18 && k < 26 && i > 2) ? 0.05 : 0.0; };
		m[c] = buildMesh(s);
		runSteps(*m[c], 400);
	}
	double dc = 0, cMax = 0, dv = 0, vMax = 0, du = 0;
	for (int b = 0; b < m[0]->nBlock; b++) {
		FVM* p = (FVM*)m[0]->blocks[m[0]->blockId[b]];
		FVM* q = (FVM*)m[1]->blocks[m[1]->blockId[b]];
		for (int ind = 0; ind < p->nAll; ind++) {
			dc = max(dc, fabs(p->concentration[ind] - q->concentration[ind]));
			cMax = max(cMax, fabs(p->concentration[ind]));
			dv = max(dv, fabs(p->WP[ind] - q->WP[ind]) + fabs(p->UP[ind] - q->UP[ind]));
			vMax = max(vMax, fabs(p->WP[ind]));
			du = max(du, fabs(p->W[ind] - q->W[ind]));
		}
	}
	report("max |c(33) - c(3)| / c max", dc / cMax);
	report("max |v(33) - v(3)| / v max", dv / vMax);
	report("max |u(33) - u(3)| / v max", du / vMax);
	CHECK(dc < 1e-10 * cMax);
	CHECK(dv < 1e-10 * vMax);
	CHECK(du < 1e-10 * vMax);
}

// randomly oriented spheroids (r = 3) settle with the averaged mobility, 1/3 (1/KT_par + 2/KT_perp)
TEST(euler_nosphere_isotropic_settling) {
	const double d = 1.0, rhoP = 20, g = 1e-4, nu = 0.02, r = 3;
	const int n = 1500;
	CaseSpec s = box("ens_iso", 33, 1, 1, 6, nu);
	s.params["gz"] = num(-g);
	s.params["initializeWithHydroPress"] = "1";
	setParticles(s, d, rhoP, r);
	s.concentration = [](int, int, int k) { return (k >= 28 && k < 40) ? 0.02 : 0.0; };
	auto mesh = buildMesh(s);
	array<double, 4> m0 = particleMoments(*mesh);
	runSteps(*mesh, n);
	array<double, 4> m1 = particleMoments(*mesh);
	double gv[3] = { 0, 0, -g };
	array<double, 3> w = slipVelocity(d, r, rhoP, ISO, gv, nu);
	double vt = -w[2];
	double mass = rhoP * PI * d * d * d / 6;
	double tau = mass * vt / ((rhoP - 1) * PI * d * d * d / 6 * g); // m / k at the terminal velocity
	double shift = vt * (n - tau);
	array<double, 3> wSphere = slipVelocity(d, 1.0, rhoP, ISO, gv, nu);
	report("terminal velocity / sphere", vt / -wSphere[2]);
	report("settling distance (analytic)", shift);
	report("settling distance", m0[3] - m1[3]);
	report("relative change of the number of particles", (m1[0] - m0[0]) / m0[0]);
	array<double, 2> v = orientationValidity(*mesh);
	report("smallest eigenvalue of A", v[0]);
	CHECK_NEAR(m1[0], m0[0], 1e-12 * m0[0]);
	CHECK_NEAR(m0[3] - m1[3], shift, 0.03 * shift);
	CHECK(v[0] > -1e-12 && v[1] < 1e-12);
}

// aligned spheroids (r = 2) inclined by 45 degrees glide sideways with the slip of one particle
// (the same as the single Lagrangian ellipsoid of problemID 22)
TEST(euler_nosphere_aligned_inclined_settling) {
	const double d = 0.75, rhoP = 20, g = 1e-4, nu = 0.02, r = 2;
	CaseSpec s = box("ens_inclined", 33, 2, 2, 2, nu);
	s.params["gz"] = num(-g);
	s.params["initializeWithHydroPress"] = "1";
	setParticles(s, d, rhoP, r);
	Sym A0 = aligned(1, 0, 1);
	s.orientation = [A0](int, int, int) { return A0; };
	s.concentration = [](int i, int j, int k) { return (i >= 6 && i < 10 && j >= 6 && j < 10 && k >= 9 && k < 13) ? 0.02 : 0.0; };
	auto mesh = buildMesh(s);
	array<double, 4> m0 = particleMoments(*mesh);
	runSteps(*mesh, 400);
	array<double, 4> m1 = particleMoments(*mesh);
	double gv[3] = { 0, 0, -g };
	array<double, 3> w = slipVelocity(d, r, rhoP, A0, gv, nu);
	int ind;
	ENS* b = (ENS*)blockAtPoint(*mesh, 8.5, 8.5, 9.5, ind);
	double slip[3] = { b->UP[ind] - b->UF[ind], b->VP[ind] - b->VF[ind], b->WP[ind] - b->WF[ind] };
	double nw = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
	double ns = sqrt(slip[0] * slip[0] + slip[1] * slip[1] + slip[2] * slip[2]);
	double angle = acos(min((slip[0] * w[0] + slip[1] * w[1] + slip[2] * w[2]) / (nw * ns), 1.0)) * 180 / PI;
	double val[3], vec[9];
	ENS::eigenSym(&b->orient[6 * ind], val, vec);
	report("expected glide angle from the vertical [deg]", atan2(fabs(w[0]), fabs(w[2])) * 180 / PI);
	report("slip / expected (magnitude)", ns / nw);
	report("angle between slip and expected [deg]", angle);
	report("alignment (largest eigenvalue of A) in the cloud", val[0]);
	// the cloud drags the liquid down, so its path is somewhat steeper than the slip
	double drift = (m1[1] - m0[1]) / (m0[3] - m1[3]);
	report("lateral drift / fall of the cloud", drift);
	report("lateral drift / fall of the cloud / (w_x / -w_z)", drift / (w[0] / -w[2]));
	CHECK_NEAR(ns, nw, 0.02 * nw);
	CHECK(angle < 1.0);
	CHECK(val[0] > 0.99);
	CHECK(drift / (w[0] / -w[2]) > 0.7 && drift / (w[0] / -w[2]) < 1.05);
}

// aligned spheroids (r = 3, C_I = 0) in a Couette flow: A = p p with p on the Jeffery orbit,
// tan(phi) = -tan(gamma t / (r + 1/r)) / r
TEST(euler_nosphere_jeffery_orbit) {
	const double nu = 0.1, U = 0.05, r = 3;
	const int n = 1000;
	CaseSpec s = couette("ens_jeffery", nu, U);
	setParticles(s, 0.5, 1.0, r);
	s.orientationL = [](int, int, int, int) { return aligned(1, 0, 0); };
	auto mesh = buildMesh(s);
	runSteps(*mesh, n);
	int ind;
	ENS* b = (ENS*)blockAtPoint(*mesh, 4.5, 4.5, 8.5, ind);
	double gamma = b->gradU[9 * ind + 2];
	double phi = atan(-tan(gamma * n / (r + 1 / r)) / r);
	Sym ex = aligned(cos(phi), 0, sin(phi));
	double err = maxDiff(&b->orient[6 * ind], ex);
	double val[3], vec[9];
	ENS::eigenSym(&b->orient[6 * ind], val, vec);
	report("du/dz / (U / 16)", gamma / (U / 16));
	report("angle of the axis (analytic) [deg]", phi * 180 / PI);
	report("angle of the axis [deg]", atan2(vec[6], vec[0]) * 180 / PI);
	report("max |A - p p|", err);
	report("1 - largest eigenvalue", 1 - val[0]);
	CHECK_NEAR(gamma, U / 16, 1e-3 * U / 16);
	CHECK(err < 1e-3);
	CHECK(1 - val[0] < 1e-6);
}

// rotary diffusion (r = 1, C_I > 0): A = I/3 + exp(-6 C_I gamma t) Q (A0 - I/3) Q^T, Q the rotation by the vorticity
TEST(euler_nosphere_rotary_diffusion) {
	const double nu = 0.1, U = 0.05, CI = 0.01;
	const int n = 1000;
	CaseSpec s = couette("ens_diffusion", nu, U);
	setParticles(s, 0.5, 1.0, 1.0);
	s.params["orientationDiffusion"] = num(CI);
	s.orientationL = [](int, int, int, int) { return aligned(1, 0, 0); };
	auto mesh = buildMesh(s);
	runSteps(*mesh, n);
	int ind;
	ENS* b = (ENS*)blockAtPoint(*mesh, 4.5, 4.5, 8.5, ind);
	double gamma = b->gradU[9 * ind + 2];
	double th = gamma * n / 2, decay = exp(-6 * CI * gamma * n);
	// p = (cos th, 0, -sin th): (A0 - I/3) rotated = p p - (I - p p) / 2 ... written out for A0 = e_x e_x
	Sym P = aligned(cos(th), 0, -sin(th));
	Sym ex;
	for (int k = 0; k < 6; k++) {
		double dev = P[k] - ISO[k]; // Q (A0 - I/3) Q^T = p p - I/3
		ex[k] = ISO[k] + decay * dev;
	}
	double err = maxDiff(&b->orient[6 * ind], ex);
	report("decay exp(-6 C_I gamma t)", decay);
	report("A_xx (analytic)", ex[0]);
	report("A_xx", b->orient[6 * ind]);
	report("max |A - analytic|", err);
	CHECK(err < 1e-3);
}

// in a uniform flow the orientation of the particles is carried with them unchanged, also where the
// concentration varies (the flux of c A uses A of the upwind cell)
TEST(euler_nosphere_uniform_orientation_is_kept) {
	CaseSpec s;
	s.name = "ens_uniform";
	s.octree = octreeBlocks(2, 2, 2, 0, [](int, int, int, int) { return false; });
	s.params = latticeParams(33, 8, 0.02);
	s.params["xPeriodicity"] = "1";
	s.params["yPeriodicity"] = "1";
	s.params["zPeriodicity"] = "1";
	setParticles(s, 0.5, 1.0, 4.0);
	s.flagL = [](int, int, int, int) { return 1; };
	const double ux = 0.03;
	s.initL = [ux](int, int, int, int) { return array<double, 4>{ 1, ux, -0.01, 0.02 }; };
	Sym A0 = aligned(0.3, -0.5, 0.8);
	for (int k = 0; k < 6; k++) {
		A0[k] = 0.7 * A0[k] + 0.3 * ISO[k];
	}
	s.orientationL = [A0](int, int, int, int) { return A0; };
	s.concentrationL = [](int, int i, int j, int k) {
		double r2 = (i - 6.5) * (i - 6.5) + (j - 8.5) * (j - 8.5) + (k - 6.5) * (k - 6.5);
		return r2 < 16 ? 0.1 * (1 - r2 / 16) : 0.0;
	};
	auto mesh = buildMesh(s);
	array<double, 4> m0 = particleMoments(*mesh);
	const int n = 100; // the cloud stays inside the periodic box (the centroid would jump)
	runSteps(*mesh, n);
	array<double, 4> m1 = particleMoments(*mesh);
	double worst = 0;
	for (int b = 0; b < mesh->nBlock; b++) {
		ENS* bl = (ENS*)mesh->blocks[mesh->blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			worst = max(worst, maxDiff(&bl->orient[6 * ind], A0));
		}
	}
	report("shift of the centroid in x / (u t)", (m1[1] - m0[1]) / (ux * n));
	report("relative change of the number of particles", (m1[0] - m0[0]) / m0[0]);
	report("max |A - A0|", worst);
	CHECK_NEAR(m1[1] - m0[1], ux * n, 0.02 * ux * n);
	CHECK_NEAR(m1[0], m0[0], 1e-12 * m0[0]);
	CHECK(worst < 1e-12);
}

// dam break: the particles of the cells that become gas go to the liquid together with their orientation
TEST(euler_nosphere_free_surface) {
	CaseSpec s = box("ens_dam", 33, 2, 1, 2, 0.05);
	s.params["gz"] = "-2e-4";
	setParticles(s, 0.2, 2.0, 2.0);
	s.flag = [](int i, int, int k) { return (i < 7 && k < 11) ? 1 : 2; };
	s.concentration = [](int i, int, int k) { return (i < 7 && k < 11) ? 5.0 : 0.0; };
	s.orientation = [](int, int, int) { return aligned(0, 0, 1); };
	auto mesh = buildMesh(s);
	double n0 = particleMoments(*mesh, false)[0];
	runSteps(*mesh, 800);
	double n1 = particleMoments(*mesh, false)[0];
	array<double, 2> v = orientationValidity(*mesh);
	report("relative change of the particles in the liquid", (n1 - n0) / n0);
	report("smallest eigenvalue of A", v[0]);
	report("max |trace A - 1|", v[1]);
	CHECK(mesh->getUmax() > 1e-3);
	CHECK_LT(fabs(n1 - n0) / n0, 1e-10);
	CHECK(v[0] > -1e-12 && v[1] < 1e-12);
}

// refined box 16^3, the half x < 8 refined: the particles (r = 1, aligned along y) settle across the
// level interface. Their number is kept and A stays aligned (spheres only turn with the liquid).
TEST(refinement_euler_nosphere_conservation) {
	CaseSpec s;
	s.name = "ref_ens";
	s.octree = octreeBlocks(2, 2, 2, 1, [](int l, int i, int, int) { return l == 0 && i == 0; });
	s.params = latticeParams(33, 8, 0.02);
	s.params["gx"] = "-1e-4";
	setParticles(s, 1.0, 20, 1.0);
	s.flagL = [](int, int, int, int) { return 1; };
	s.concentrationL = [](int l, int i, int j, int k) {
		return (l == 0 && i >= 9 && i < 15 && j >= 4 && j < 12 && k >= 4 && k < 12) ? 0.05 : 0.0;
	};
	s.orientationL = [](int, int, int, int) { return aligned(0, 1, 0); };
	auto mesh = buildMesh(s);
	double n0 = particleMoments(*mesh)[0];
	double worst = 0;
	for (int k = 0; k < 5; k++) {
		runSteps(*mesh, 100, 100 * k);
		worst = max(worst, fabs(particleMoments(*mesh)[0] - n0) / n0);
	}
	double fine = 0, minAlign = 1;
	for (int b = 0; b < mesh->nBlock; b++) {
		ENS* bl = (ENS*)mesh->blocks[mesh->blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			if (bl->boundary[ind] || bl->concentration[ind] < 1e-4) {
				continue;
			}
			if (bl->level > 0) {
				fine += bl->concentration[ind];
			}
			double val[3], vec[9];
			ENS::eigenSym(&bl->orient[6 * ind], val, vec);
			minAlign = min(minAlign, val[0]);
		}
	}
	array<double, 2> v = orientationValidity(*mesh);
	report("relative change of the number of particles", worst);
	report("part of the particles in the finer half", fine / n0);
	report("smallest alignment of the cells with particles", minAlign);
	report("smallest eigenvalue of A", v[0]);
	CHECK(fine / n0 > 0.2);
	CHECK_LT(worst, 1e-11);
	CHECK(minAlign > 0.95);
	CHECK(v[0] > -1e-12 && v[1] < 1e-12);
}
