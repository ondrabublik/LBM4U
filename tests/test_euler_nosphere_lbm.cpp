// Tests of the Eulerian non-spherical particles transported by LBM (problemID 42,
// BlockFreeSurfaceEulerNoSphereParticlesLBM): the sphere limit (problemID 4), settling with the averaged mobility,
// Jeffery orbits and rotary diffusion in a Couette flow, transport of the orientation with the populations
// (uniform flow, conservation, free surface, refined mesh).
// Lattice units equal physical units (latticeParams).

#include "TestFramework.h"
#include "TestCase.h"
#include <iostream>

using namespace std;

namespace {
	const double PI = 3.14159265358979323846;
	typedef BlockFreeSurfaceEulerNoSphereParticlesLBM NSL;
	typedef BlockFreeSurfaceEulerParticlesLBM LBM;
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

	void setParticles(CaseSpec& s, double d, double rhoP, double r, double diffusivity) {
		s.params["particleDiameter"] = num(d);
		s.params["particleDensity"] = num(rhoP);
		s.params["particleAspectRatio"] = num(r);
		s.params["particleDiffusivity"] = num(diffusivity);
		s.params["particleInterForces"] = "0";
		s.params["overpackDiffusivity"] = "0";
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
		s.params = latticeParams(42, 8, nu);
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

	bool liquid(BlockFreeSurface* b, int ind) {
		return !b->boundary[ind] && (b->flagMap[ind] == 1 || b->flagMap[ind] == 3);
	}

	// smallest eigenvalue and |trace - 1| of the orientation tensors of the liquid cells
	array<double, 2> orientationValidity(Mesh& mesh) {
		array<double, 2> r = { 1, 0 };
		for (int b = 0; b < mesh.nBlock; b++) {
			NSL* bl = (NSL*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				if (!liquid(bl, ind)) {
					continue;
				}
				const double* A = &bl->orient[6 * ind];
				double val[3], vec[9];
				SpheroidOrientation::eigenSym(A, val, vec);
				r[0] = min(r[0], val[2]);
				r[1] = max(r[1], fabs(A[0] + A[1] + A[2] - 1));
			}
		}
		return r;
	}

	// sum of the concentration and its centroid (cells of the block's level)
	array<double, 4> particleMoments(Mesh& mesh) {
		array<double, 4> m = { 0, 0, 0, 0 };
		for (int b = 0; b < mesh.nBlock; b++) {
			LBM* bl = (LBM*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				if (bl->boundary[ind]) {
					continue;
				}
				double c = bl->concentration[ind];
				int gi, gj, gk;
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

	// conserved quantities of the transport: sum of the populations S (theta times the number of particles) and
	// the measure sum(S A) of all computed cells (the populations hold numbers of particles on every level)
	array<double, 7> populationMeasure(Mesh& mesh) {
		array<double, 7> m = { 0, 0, 0, 0, 0, 0, 0 };
		for (int b = 0; b < mesh.nBlock; b++) {
			NSL* bl = (NSL*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				if (bl->boundary[ind]) {
					continue;
				}
				double S = 0;
				for (int p = 0; p < bl->nf; p++) {
					S += bl->gc[bl->nf * ind + p];
				}
				m[0] += S;
				for (int k = 0; k < 6; k++) {
					m[1 + k] += S * bl->orient[6 * ind + k];
				}
			}
		}
		return m;
	}

	double measureChange(const array<double, 7>& a, const array<double, 7>& b) {
		double d = 0;
		for (int k = 0; k < 7; k++) {
			d = max(d, fabs(a[k] - b[k]));
		}
		return d / a[0];
	}

	// semi-axes of the spheroid of aspect ratio r with the volume of the sphere of diameter d
	void spheroid(double d, double r, double a[3]) {
		double b = d / 2 / cbrt(r);
		a[0] = r * b;
		a[1] = a[2] = b;
	}

	// settling slip velocity of spheroids with the orientation tensor A: w = <M> F / g(Re)
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

// r = 1: the same result as problemID 4 (to round-off, the drag is computed as a tensor)
TEST(euler_nosphere_lbm_sphere_limit) {
	const double d = 1.0, rhoP = 20, g = 1e-4, nu = 0.02;
	unique_ptr<Mesh> m[2];
	for (int c = 0; c < 2; c++) {
		CaseSpec s = box(c ? "nsl_sphere42" : "nsl_sphere4", c ? 42 : 4, 1, 1, 4, nu);
		s.params["gz"] = num(-g);
		s.params["gx"] = num(0.3 * g);
		s.params["initializeWithHydroPress"] = "1";
		setParticles(s, d, rhoP, 1.0, 0.005);
		s.concentration = [](int i, int, int k) { return (k >= 18 && k < 26 && i > 2) ? 0.05 : 0.0; };
		m[c] = buildMesh(s);
		runSteps(*m[c], 400);
	}
	double dc = 0, cMax = 0, dv = 0, vMax = 0, du = 0;
	for (int b = 0; b < m[0]->nBlock; b++) {
		LBM* p = (LBM*)m[0]->blocks[m[0]->blockId[b]];
		LBM* q = (LBM*)m[1]->blocks[m[1]->blockId[b]];
		for (int ind = 0; ind < p->nAll; ind++) {
			dc = max(dc, fabs(p->concentration[ind] - q->concentration[ind]));
			cMax = max(cMax, fabs(p->concentration[ind]));
			dv = max(dv, fabs(p->WP[ind] - q->WP[ind]) + fabs(p->UP[ind] - q->UP[ind]));
			vMax = max(vMax, fabs(p->WP[ind]));
			du = max(du, fabs(p->W[ind] - q->W[ind]));
		}
	}
	report("max |c(42) - c(4)| / c max", dc / cMax);
	report("max |v(42) - v(4)| / v max", dv / vMax);
	report("max |u(42) - u(4)| / v max", du / vMax);
	CHECK(dc < 1e-10 * cMax);
	CHECK(dv < 1e-10 * vMax);
	CHECK(du < 1e-10 * vMax);
}

// randomly oriented spheroids (r = 3) settle with the averaged mobility, 1/3 (1/KT_par + 2/KT_perp)
TEST(euler_nosphere_lbm_isotropic_settling) {
	const double d = 1.0, rhoP = 20, g = 1e-4, nu = 0.02, r = 3;
	const int n = 1500;
	CaseSpec s = box("nsl_iso", 42, 1, 1, 6, nu);
	s.params["gz"] = num(-g);
	s.params["initializeWithHydroPress"] = "1";
	setParticles(s, d, rhoP, r, 0.005);
	s.concentration = [](int, int, int k) { return (k >= 28 && k < 40) ? 0.02 : 0.0; };
	auto mesh = buildMesh(s);
	array<double, 4> m0 = particleMoments(*mesh);
	array<double, 7> q0 = populationMeasure(*mesh);
	runSteps(*mesh, n);
	array<double, 4> m1 = particleMoments(*mesh);
	array<double, 7> q1 = populationMeasure(*mesh);
	double gv[3] = { 0, 0, -g };
	array<double, 3> w = slipVelocity(d, r, rhoP, ISO, gv, nu);
	double vt = -w[2];
	double mass = rhoP * PI * d * d * d / 6;
	double tau = mass * vt / ((rhoP - 1) * PI * d * d * d / 6 * g); // m / k at the terminal velocity
	double shift = vt * (n - tau);
	array<double, 3> wSphere = slipVelocity(d, 1.0, rhoP, ISO, gv, nu);
	array<double, 2> v = orientationValidity(*mesh);
	report("terminal velocity / sphere", vt / -wSphere[2]);
	report("settling distance (analytic)", shift);
	report("settling distance", m0[3] - m1[3]);
	report("relative change of the sum of the populations", (q1[0] - q0[0]) / q0[0]);
	report("relative change of the number of particles", (m1[0] - m0[0]) / m0[0]);
	report("smallest eigenvalue of A", v[0]);
	CHECK_NEAR(q1[0], q0[0], 1e-12 * q0[0]);
	CHECK_NEAR(m1[0], m0[0], 1e-6 * m0[0]);
	CHECK_NEAR(m0[3] - m1[3], shift, 0.03 * shift);
	CHECK(v[0] > -1e-12 && v[1] < 1e-12);
}

// aligned spheroids (r = 2) inclined by 45 degrees glide sideways with the slip of one particle
TEST(euler_nosphere_lbm_aligned_inclined_settling) {
	const double d = 0.75, rhoP = 20, g = 1e-4, nu = 0.02, r = 2;
	CaseSpec s = box("nsl_inclined", 42, 2, 2, 2, nu);
	s.params["gz"] = num(-g);
	s.params["initializeWithHydroPress"] = "1";
	setParticles(s, d, rhoP, r, 0.005);
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
	NSL* b = (NSL*)blockAtPoint(*mesh, 8.5, 8.5, 9.5, ind);
	double slip[3] = { b->UP[ind] - b->UF[ind], b->VP[ind] - b->VF[ind], b->WP[ind] - b->WF[ind] };
	double nw = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
	double ns = sqrt(slip[0] * slip[0] + slip[1] * slip[1] + slip[2] * slip[2]);
	double angle = acos(min((slip[0] * w[0] + slip[1] * w[1] + slip[2] * w[2]) / (nw * ns), 1.0)) * 180 / PI;
	double val[3], vec[9];
	SpheroidOrientation::eigenSym(&b->orient[6 * ind], val, vec);
	double drift = (m1[1] - m0[1]) / (m0[3] - m1[3]);
	report("expected glide angle from the vertical [deg]", atan2(fabs(w[0]), fabs(w[2])) * 180 / PI);
	report("slip / expected (magnitude)", ns / nw);
	report("angle between slip and expected [deg]", angle);
	report("alignment (largest eigenvalue of A) in the cloud", val[0]);
	report("lateral drift / fall of the cloud / (w_x / -w_z)", drift / (w[0] / -w[2]));
	CHECK_NEAR(ns, nw, 0.02 * nw);
	CHECK(angle < 1.0);
	CHECK(val[0] > 0.99);
	CHECK(drift / (w[0] / -w[2]) > 0.7 && drift / (w[0] / -w[2]) < 1.05);
}

// aligned spheroids (r = 3, C_I = 0) in a Couette flow: A = p p with p on the Jeffery orbit,
// tan(phi) = -tan(gamma t / (r + 1/r)) / r; the particles fill the channel (transport of a uniform A)
TEST(euler_nosphere_lbm_jeffery_orbit) {
	const double nu = 0.1, U = 0.05, r = 3;
	const int n = 1000;
	CaseSpec s = couette("nsl_jeffery", nu, U);
	setParticles(s, 0.5, 1.0, r, 0.01);
	s.concentrationL = [](int, int, int, int) { return 0.01; };
	s.orientationL = [](int, int, int, int) { return aligned(1, 0, 0); };
	auto mesh = buildMesh(s);
	runSteps(*mesh, n);
	int ind;
	NSL* b = (NSL*)blockAtPoint(*mesh, 4.5, 4.5, 8.5, ind);
	double gamma = b->gradU[9 * ind + 2];
	double phi = atan(-tan(gamma * n / (r + 1 / r)) / r);
	Sym ex = aligned(cos(phi), 0, sin(phi));
	double err = maxDiff(&b->orient[6 * ind], ex);
	double val[3], vec[9];
	SpheroidOrientation::eigenSym(&b->orient[6 * ind], val, vec);
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
TEST(euler_nosphere_lbm_rotary_diffusion) {
	const double nu = 0.1, U = 0.05, CI = 0.01;
	const int n = 1000;
	CaseSpec s = couette("nsl_diffusion", nu, U);
	setParticles(s, 0.5, 1.0, 1.0, 0.01);
	s.params["orientationDiffusion"] = num(CI);
	s.concentrationL = [](int, int, int, int) { return 0.01; };
	s.orientationL = [](int, int, int, int) { return aligned(1, 0, 0); };
	auto mesh = buildMesh(s);
	runSteps(*mesh, n);
	int ind;
	NSL* b = (NSL*)blockAtPoint(*mesh, 4.5, 4.5, 8.5, ind);
	double gamma = b->gradU[9 * ind + 2];
	double th = gamma * n / 2, decay = exp(-6 * CI * gamma * n);
	Sym P = aligned(cos(th), 0, -sin(th));
	Sym ex;
	for (int k = 0; k < 6; k++) {
		ex[k] = ISO[k] + decay * (P[k] - ISO[k]); // Q (A0 - I/3) Q^T = p p - I/3
	}
	double err = maxDiff(&b->orient[6 * ind], ex);
	report("decay exp(-6 C_I gamma t)", decay);
	report("A_xx (analytic)", ex[0]);
	report("A_xx", b->orient[6 * ind]);
	report("max |A - analytic|", err);
	CHECK(err < 1e-3);
}

namespace {
	// periodic box 16^3 with a uniform flow, a cloud of neutrally buoyant particles
	CaseSpec uniformFlow(const string& name, double ux) {
		CaseSpec s;
		s.name = name;
		s.octree = octreeBlocks(2, 2, 2, 0, [](int, int, int, int) { return false; });
		s.params = latticeParams(42, 8, 0.02);
		s.params["xPeriodicity"] = "1";
		s.params["yPeriodicity"] = "1";
		s.params["zPeriodicity"] = "1";
		setParticles(s, 0.5, 1.0, 4.0, 0.01);
		s.flagL = [](int, int, int, int) { return 1; };
		s.initL = [ux](int, int, int, int) { return array<double, 4>{ 1, ux, -0.01, 0.02 }; };
		s.concentrationL = [](int, int i, int j, int k) {
			double r2 = (i - 6.5) * (i - 6.5) + (j - 8.5) * (j - 8.5) + (k - 6.5) * (k - 6.5);
			return r2 < 16 ? 0.1 * (1 - r2 / 16) : 0.0;
		};
		return s;
	}
}

// in a uniform flow the orientation is carried with the particles unchanged, also where the concentration varies.
// Small particles: the fluid fraction 1 - c V of the cloud makes the flow slightly non-uniform (|u - u0| ~ 2e-9 here)
// and the velocity gradient turns the spheroids (Folgar-Tucker), with d = 0.5 |A - A0| would be 2e-4.
TEST(euler_nosphere_lbm_uniform_orientation_is_kept) {
	const double ux = 0.03;
	CaseSpec s = uniformFlow("nsl_uniform", ux);
	s.params["particleDiameter"] = num(0.02);
	Sym A0 = aligned(0.3, -0.5, 0.8);
	for (int k = 0; k < 6; k++) {
		A0[k] = 0.7 * A0[k] + 0.3 * ISO[k];
	}
	s.orientationL = [A0](int, int, int, int) { return A0; };
	auto mesh = buildMesh(s);
	array<double, 4> m0 = particleMoments(*mesh);
	const int n = 100; // the cloud stays inside the periodic box (the centroid would jump)
	runSteps(*mesh, n);
	array<double, 4> m1 = particleMoments(*mesh);
	double worst = 0;
	for (int b = 0; b < mesh->nBlock; b++) {
		NSL* bl = (NSL*)mesh->blocks[mesh->blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			worst = max(worst, maxDiff(&bl->orient[6 * ind], A0));
		}
	}
	report("shift of the centroid in x / (u t)", (m1[1] - m0[1]) / (ux * n));
	report("relative change of the number of particles", (m1[0] - m0[0]) / m0[0]);
	report("max |A - A0|", worst);
	CHECK_NEAR(m1[1] - m0[1], ux * n, 0.02 * ux * n);
	CHECK_NEAR(m1[0], m0[0], 1e-9 * m0[0]);
	CHECK(worst < 1e-7);
}

namespace {
	// liquid at rest, Gaussian cloud of neutrally buoyant particles diffusing (D = 0.1), aligned along x for
	// x < 8 and along z for x > 8; the velocity gradient is zero, so A changes by the transport only
	void restCloud(CaseSpec& s) {
		setParticles(s, 0.5, 1.0, 3.0, 0.1);
		s.concentrationL = [](int l, int i, int j, int k) {
			double x = (i + 0.5) / (1 << l), y = (j + 0.5) / (1 << l), z = (k + 0.5) / (1 << l);
			double r2 = (x - 7.3) * (x - 7.3) + (y - 8.4) * (y - 8.4) + (z - 7.8) * (z - 7.8);
			return 0.05 * exp(-r2 / (2 * 2.0 * 2.0)) / (1 << (3 * l)); // number of particles of the cell
		};
		s.orientationL = [](int l, int i, int, int) { return (i >> l) < 8 ? aligned(1, 0, 0) : aligned(0, 0, 1); };
	}

	// A_xx part sum(S A_xx) / sum(S) of the particles in the cells with x > 8
	double rightAxx(Mesh& mesh) {
		double S = 0, Sxx = 0;
		for (int b = 0; b < mesh.nBlock; b++) {
			NSL* bl = (NSL*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				int gi, gj, gk;
				bl->globalCoords(ind, gi, gj, gk);
				if (bl->boundary[ind] || (gi >> bl->level) < 8) {
					continue;
				}
				double sum = 0;
				for (int p = 0; p < bl->nf; p++) {
					sum += bl->gc[bl->nf * ind + p];
				}
				S += sum;
				Sxx += sum * bl->orient[6 * ind];
			}
		}
		return Sxx / S;
	}

	// The measure sum(S A) is conserved exactly where the populations are non-negative; the over-relaxed first
	// collisions of the initial Gaussian give small negative populations in its tails (there A is mixed from the
	// positive ones only), afterwards the measure has to be kept to round-off.
	void measureConservation(Mesh& mesh, const string& what) {
		array<double, 7> q0 = populationMeasure(mesh);
		double right0 = rightAxx(mesh);
		runSteps(mesh, 20);
		array<double, 7> q1 = populationMeasure(mesh);
		runSteps(mesh, 100, 20);
		array<double, 7> q2 = populationMeasure(mesh);
		array<double, 2> v = orientationValidity(mesh);
		report(what + ": fluid velocity (rest)", mesh.getUmax());
		report(what + ": relative change of the sum of the populations", (q2[0] - q0[0]) / q0[0]);
		report(what + ": relative change of sum(S A), steps 0 - 20", measureChange(q0, q1));
		report(what + ": relative change of sum(S A), steps 20 - 120", measureChange(q1, q2));
		report(what + ": A_xx part of the particles with x > 8 (initial)", right0);
		report(what + ": A_xx part of the particles with x > 8", rightAxx(mesh));
		report(what + ": smallest eigenvalue of A", v[0]);
		CHECK_NEAR(q2[0], q0[0], 1e-12 * q0[0]);
		CHECK_LT(measureChange(q0, q1), 1e-4);
		CHECK_LT(measureChange(q1, q2), 1e-12);
		CHECK(rightAxx(mesh) > 0.1); // particles aligned along x crossed x = 8
		CHECK(v[0] > -1e-12 && v[1] < 1e-12);
	}
}

TEST(euler_nosphere_lbm_orientation_measure_is_conserved) {
	CaseSpec s = uniformFlow("nsl_measure", 0.0);
	s.initL = [](int, int, int, int) { return array<double, 4>{ 1, 0, 0, 0 }; };
	restCloud(s);
	auto mesh = buildMesh(s);
	measureConservation(*mesh, "uniform mesh");
}

// The orientation moves only with a net flux of particles: at rest with a uniform concentration a step of A stays
// sharp (an exchange of the gross populations would smear it over ~sqrt(2 t / (6 theta)) cells, 4 cells here),
// with a uniform flow the step moves with the flow and is smeared only by the upwind error.
TEST(euler_nosphere_lbm_orientation_needs_particle_flux) {
	for (double ux : { 0.0, 0.02 }) {
		CaseSpec s = uniformFlow(ux == 0 ? "nsl_step_rest" : "nsl_step_flow", ux);
		s.initL = [ux](int, int, int, int) { return array<double, 4>{ 1, ux, 0, 0 }; };
		s.params["particleDiameter"] = num(0.02);
		s.concentrationL = [](int, int, int, int) { return 0.01; };
		s.orientationL = [](int, int i, int, int) { return (i >= 4 && i < 8) ? aligned(0, 0, 1) : aligned(1, 0, 0); };
		auto mesh = buildMesh(s);
		const int n = 100;
		runSteps(*mesh, n);
		// A_zz along x (y = z = 8), exact: the step 4 <= x < 8 shifted by ux t
		double err = 0, width = 0;
		for (int i = 0; i < 16; i++) {
			int ind;
			NSL* b = (NSL*)mesh->cellAt(i, 8, 8, ind);
			double x = i + 0.5 - ux * n;
			double ex = (x >= 4 && x < 8) ? 1.0 : 0.0;
			double azz = b->orient[6 * ind + 2];
			err = max(err, fabs(azz - ex));
			width += azz * (1 - azz); // cells with mixed A
		}
		report(string(ux == 0 ? "rest" : "flow 0.02") + ": sum A_zz (1 - A_zz) along x (mixed cells)", width);
		if (ux == 0) {
			report("rest: max |A_zz - step|", err);
			CHECK(err < 1e-12);
		}
		else {
			// first order upwind: two fronts with sigma = sqrt(u (1 - u) t) = 1.4 cells give about 1.6, the exchange of the
			// gross populations (sigma ~ 6.5 cells) about 7
			CHECK(width < 2.5);
		}
	}
}

// dam break: the particles of the cells that become gas go to the liquid together with their orientation
TEST(euler_nosphere_lbm_free_surface) {
	CaseSpec s = box("nsl_dam", 42, 2, 1, 2, 0.05);
	s.params["gz"] = "-2e-4";
	setParticles(s, 0.2, 2.0, 2.0, 0.01);
	s.flag = [](int i, int, int k) { return (i < 7 && k < 11) ? 1 : 2; };
	s.concentration = [](int i, int, int k) { return (i < 7 && k < 11) ? 5.0 : 0.0; };
	s.orientation = [](int, int, int) { return aligned(0, 0, 1); };
	auto mesh = buildMesh(s);
	array<double, 7> q0 = populationMeasure(*mesh);
	runSteps(*mesh, 800);
	array<double, 7> q1 = populationMeasure(*mesh);
	array<double, 2> v = orientationValidity(*mesh);
	report("relative change of the sum of the populations", (q1[0] - q0[0]) / q0[0]);
	report("smallest eigenvalue of A", v[0]);
	report("max |trace A - 1|", v[1]);
	CHECK(mesh->getUmax() > 1e-3);
	CHECK_LT(fabs(q1[0] - q0[0]) / q0[0], 1e-12);
	CHECK(v[0] > -1e-12 && v[1] < 1e-12);
}

namespace {
	// refined box 16^3, the half x < 8 refined
	CaseSpec refinedBox(const string& name) {
		CaseSpec s;
		s.name = name;
		s.octree = octreeBlocks(2, 2, 2, 1, [](int l, int i, int, int) { return l == 0 && i == 0; });
		s.params = latticeParams(42, 8, 0.02);
		s.flagL = [](int, int, int, int) { return 1; };
		return s;
	}
}

// the same on the refined mesh (x < 8 refined): the populations exchanged with the other level carry the
// orientation of the cell they leave
TEST(refinement_euler_nosphere_lbm_measure) {
	CaseSpec s = refinedBox("ref_nsl_measure");
	restCloud(s);
	auto mesh = buildMesh(s);
	measureConservation(*mesh, "refined mesh");
}

// the particles (r = 1, aligned along y) settle across the level interface. Their number is kept and A stays
// aligned (spheres only turn with the liquid).
TEST(refinement_euler_nosphere_lbm_settling) {
	CaseSpec s = refinedBox("ref_nsl");
	s.params["gx"] = "-1e-4";
	setParticles(s, 1.0, 20, 1.0, 0.005);
	s.concentrationL = [](int l, int i, int j, int k) {
		return (l == 0 && i >= 9 && i < 15 && j >= 4 && j < 12 && k >= 4 && k < 12) ? 0.05 : 0.0;
	};
	s.orientationL = [](int, int, int, int) { return aligned(0, 1, 0); };
	auto mesh = buildMesh(s);
	double n0 = populationMeasure(*mesh)[0];
	double worst = 0;
	for (int k = 0; k < 5; k++) {
		runSteps(*mesh, 100, 100 * k);
		worst = max(worst, fabs(populationMeasure(*mesh)[0] - n0) / n0);
	}
	double fine = 0, total = 0, minAlign = 1;
	for (int b = 0; b < mesh->nBlock; b++) {
		NSL* bl = (NSL*)mesh->blocks[mesh->blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			if (bl->boundary[ind]) {
				continue;
			}
			total += bl->concentration[ind];
			if (bl->concentration[ind] < 1e-4) {
				continue;
			}
			if (bl->level > 0) {
				fine += bl->concentration[ind];
			}
			double val[3], vec[9];
			SpheroidOrientation::eigenSym(&bl->orient[6 * ind], val, vec);
			minAlign = min(minAlign, val[0]);
		}
	}
	array<double, 2> v = orientationValidity(*mesh);
	report("relative change of the sum of the populations", worst);
	report("part of the particles in the finer half", fine / total);
	report("smallest alignment of the cells with particles", minAlign);
	report("smallest eigenvalue of A", v[0]);
	CHECK(fine / total > 0.2);
	CHECK_LT(worst, 1e-11);
	CHECK(minAlign > 0.95);
	CHECK(v[0] > -1e-12 && v[1] < 1e-12);
}
