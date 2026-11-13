// Tests of the particle models: Lagrangian particles (problemID 2) and the Eulerian concentration
// solved with FVM (problemID 3) and with LBM (problemID 4).
// Geometry: closed or periodic box of 8^3 blocks, lattice units equal physical units unless stated.

#include "TestFramework.h"
#include "TestCase.h"
#include <omp.h>
#include <iostream>

using namespace std;

namespace {
	typedef BlockFreeSurfaceLagrangeParticles LagBlock;
	typedef BlockFreeSurfaceEulerParticlesFVM FvmBlock;
	typedef BlockFreeSurfaceEulerParticlesLBM LbmBlock;

	void report(const string& what, double value) {
		cout << "       " << what << " = " << value << endl;
	}

	// terminal slip velocity of a sphere with the Schiller-Naumann drag (lattice units, fluid density 1)
	double terminalVelocity(double d, double rhoP, double g, double nu) {
		double v = (rhoP - 1) * g * d * d / (18 * nu);
		for (int it = 0; it < 200; it++) {
			v = (rhoP - 1) * g * d * d / (18 * nu * (1 + 0.15 * pow(d * v / nu, 0.687)));
		}
		return v;
	}

	// relaxation time of the particle velocity at the terminal velocity
	double relaxationTime(double d, double rhoP, double g, double nu) {
		double v = terminalVelocity(d, rhoP, g, nu);
		return rhoP * d * d / (18 * nu * (1 + 0.15 * pow(d * v / nu, 0.687)));
	}

	// box of nBi x nBj x nBk blocks of 8^3 cells full of liquid, the exterior of the domain is a wall
	CaseSpec box(const string& name, int problemID, int nBi, int nBj, int nBk, double nu) {
		CaseSpec s;
		s.name = name;
		s.nBi = nBi;
		s.nBj = nBj;
		s.nBk = nBk;
		s.ni = s.nj = s.nk = 8;
		s.params = latticeParams(problemID, 8, nu);
		s.flag = [](int, int, int) { return 1; };
		return s;
	}

	void setGravity(CaseSpec& s, double gz) {
		ostringstream os;
		os.precision(17);
		os << gz;
		s.params["gz"] = os.str();
		s.params["initializeWithHydroPress"] = "1";
	}

	void setEulerParticles(CaseSpec& s, double d, double rhoP, double diffusivity) {
		ostringstream os;
		os.precision(17);
		os << d;
		s.params["particleDiameter"] = os.str();
		os.str("");
		os << rhoP;
		s.params["particleDensity"] = os.str();
		os.str("");
		os << diffusivity;
		s.params["particleDiffusivity"] = os.str();
		s.params["particleInterForces"] = "0";
		s.params["overpackDiffusivity"] = "0";
	}

	// momentum held by the post-collision populations
	array<double, 3> populationMomentum(Mesh& mesh) {
		array<double, 3> m = { 0, 0, 0 };
		for (int b = 0; b < mesh.nBlock; b++) {
			Block* bl = mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				for (int p = 0; p < bl->nf; p++) {
					double f = bl->fc[bl->nf * ind + p];
					m[0] += bl->ex[p] * f;
					m[1] += bl->ey[p] * f;
					m[2] += bl->ez[p] * f;
				}
			}
		}
		return m;
	}

	array<double, 3> spreadForce(Mesh& mesh) {
		array<double, 3> f = { 0, 0, 0 };
		for (int b = 0; b < mesh.nBlock; b++) {
			LagBlock* bl = (LagBlock*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				f[0] += bl->particleForceX[ind];
				f[1] += bl->particleForceY[ind];
				f[2] += bl->particleForceZ[ind];
			}
		}
		return f;
	}

	// sum, centre and variance (in z) of the concentration and its extremes
	struct Moments {
		double sum = 0, zc = 0, var = 0, cMin = 1e300, cMax = -1e300;
	};

	template <class B> Moments moments(Mesh& mesh) {
		Moments m;
		for (int b = 0; b < mesh.nBlock; b++) {
			B* bl = (B*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				if (bl->boundary[ind]) {
					continue;
				}
				int gi, gj, gk;
				bl->globalCoords(ind, gi, gj, gk);
				double c = bl->concentration[ind];
				double z = gk + 0.5;
				m.sum += c;
				m.zc += c * z;
				m.var += c * z * z;
				m.cMin = min(m.cMin, c);
				m.cMax = max(m.cMax, c);
			}
		}
		m.zc /= m.sum;
		m.var = m.var / m.sum - m.zc * m.zc;
		return m;
	}

	// sum of the scalar populations of the liquid cells (the quantity conserved by the scalar LBM)
	double scalarSum(Mesh& mesh) {
		double s = 0;
		for (int b = 0; b < mesh.nBlock; b++) {
			LbmBlock* bl = (LbmBlock*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				if (bl->flagMap[ind] == 1 || bl->flagMap[ind] == 3) {
					for (int p = 0; p < bl->nf; p++) {
						s += bl->gc[bl->nf * ind + p];
					}
				}
			}
		}
		return s;
	}

	double totalMass(Mesh& m, int iter) {
		m.evalBlocks(iter, 105);
		return m.massTot;
	}
}

// ---------------------------------------------------------------------------------------------
// Lagrangian particles
// ---------------------------------------------------------------------------------------------

// a single sphere settles in a closed box: the slip velocity has to be the terminal velocity
TEST(lagrange_terminal_velocity) {
	const double d = 0.8, rhoP = 20, g = 1e-4, nu = 0.02;
	CaseSpec s = box("lag_terminal", 2, 2, 2, 2, nu);
	setGravity(s, -g);
	s.particles.push_back({ 8.3, 7.6, 12.4, d, rhoP });
	auto mesh = buildMesh(s);
	LagrangeParticles& P = *mesh->particles;
	runSteps(*mesh, 1200);
	double vt = terminalVelocity(d, rhoP, g, nu);
	double slip = P.veloFZ[0] - P.veloZ[0];
	report("terminal velocity (analytic)", vt);
	report("slip velocity", slip);
	report("fluid velocity at the particle / vt", P.veloFZ[0] / vt);
	report("z", P.partZ[0]);
	CHECK_NEAR(slip, vt, 0.01 * vt);
	// the particle crossed the block boundary at z = 8 and stays on its vertical line
	CHECK(P.partZ[0] < 8.0 && P.partZ[0] > 4.0);
	CHECK_NEAR(P.partX[0], 8.3, 0.02);
	CHECK_NEAR(P.partY[0], 7.6, 0.02);
	CHECK(mesh->getUmax() < 2 * vt);
}

// periodic box without gravity: the momentum lost by the particle is gained by the fluid
TEST(lagrange_momentum_conservation) {
	for (int nIntIter : { 1, 4 }) { // sub-steps of the particles in a LBM step
		CaseSpec s = box("lag_momentum", 2, 2, 2, 2, 0.05);
		s.params["xPeriodicity"] = "1";
		s.params["yPeriodicity"] = "1";
		s.params["zPeriodicity"] = "1";
		s.params["nIntIter"] = to_string(nIntIter);
		s.particles.push_back({ 8.3, 7.6, 8.1, 0.8, 5.0 });
		auto mesh = buildMesh(s);
		LagrangeParticles& P = *mesh->particles;
		P.veloX[0] = 0.03;
		P.veloY[0] = -0.02;
		P.veloZ[0] = 0.01;
		double m = P.partMass[0];
		double p0[3] = { m * 0.03, -m * 0.02, m * 0.01 };
		double worst = 0;
		for (int it = 0; it < 100; it++) {
			mesh->timeStep(it);
			array<double, 3> pf = populationMomentum(*mesh);
			double pp[3] = { m * P.veloX[0], m * P.veloY[0], m * P.veloZ[0] };
			for (int dd = 0; dd < 3; dd++) {
				worst = max(worst, fabs(pf[dd] + pp[dd] - p0[dd]) / fabs(p0[0]));
			}
		}
		report("relative momentum error with " + to_string(nIntIter) + " sub-step(s)", worst);
		CHECK(fabs(P.veloX[0]) < 0.003); // the particle was stopped by the drag
		CHECK_LT(worst, 1e-9);
	}
}

// a particle leaving a periodic box enters it at the opposite side
TEST(lagrange_periodic_box) {
	CaseSpec s = box("lag_periodic", 2, 2, 2, 2, 0.05);
	s.params["xPeriodicity"] = "1";
	s.params["yPeriodicity"] = "1";
	s.params["zPeriodicity"] = "1";
	s.params["gx"] = "1e-4";
	s.particles.push_back({ 15.9, 7.6, 8.1, 0.5, 1.0 }); // neutrally buoyant: it follows the accelerating fluid
	auto mesh = buildMesh(s);
	LagrangeParticles& P = *mesh->particles;
	runSteps(*mesh, 200);
	report("x after crossing the periodic boundary", P.partX[0]);
	CHECK(!P.inGas[0]);
	CHECK(P.partX[0] > 0.5 && P.partX[0] < 8.0);
	CHECK(P.veloX[0] > 0.015 && P.veloX[0] < 0.021); // g t = 0.02
}

// the force of a particle next to a wall is spread to the fluid cells only, nothing may be lost
TEST(lagrange_force_spreading_near_wall) {
	CaseSpec s = box("lag_spreading", 2, 2, 2, 2, 0.02);
	setGravity(s, -1e-4);
	s.particles.push_back({ 8.3, 7.6, 0.6, 0.4, 20 });  // above the floor
	s.particles.push_back({ 0.7, 0.6, 8.2, 0.4, 20 });  // in a vertical edge of the box
	auto mesh = buildMesh(s);
	LagrangeParticles& P = *mesh->particles;
	runSteps(*mesh, 1);
	P.applyForce();
	array<double, 3> f = spreadForce(*mesh);
	double fz = P.Gz[0] + P.Gz[1];
	report("spread force / particle force", f[2] / fz);
	CHECK(fabs(fz) > 1e-8);
	CHECK_NEAR(f[2], fz, 1e-10 * fabs(fz));
}

// writing the results must not influence the computation
TEST(lagrange_output_does_not_change_result) {
	vector<double> ref;
	for (int withOutput = 0; withOutput < 2; withOutput++) {
		CaseSpec s = box("lag_output", 2, 2, 2, 2, 0.02);
		// physical units: cell 1 mm, time step 0.1 ms
		s.params["H"] = "0.008";
		s.params["uMax"] = "1";
		s.params["vKin"] = "2e-4";
		s.params["density"] = "1000";
		setGravity(s, -9.81);
		// two particles closer than the range of the repulsive force
		s.particles.push_back({ 8.2e-3, 8.2e-3, 8.3e-3, 5e-4, 2500 });
		s.particles.push_back({ 8.8e-3, 8.2e-3, 8.3e-3, 5e-4, 2500 });
		auto mesh = buildMesh(s);
		LagrangeParticles& P = *mesh->particles;
		for (int it = 0; it < 60; it++) {
			mesh->timeStep(it);
			if (withOutput && it % 10 == 0) {
				mesh->writeResultsToVTK(it);
			}
		}
		vector<double> v = { P.partX[0], P.partX[1], P.partZ[0], P.partZ[1] };
		CHECK(P.partX[1] - P.partX[0] > 0.6001); // they repel each other
		if (ref.empty()) {
			ref = v;
		}
		else {
			report("distance without output", ref[1] - ref[0]);
			report("distance with output", v[1] - v[0]);
			for (size_t i = 0; i < v.size(); i++) {
				CHECK_NEAR(v[i], ref[i], 1e-12);
			}
		}
	}
}

// a particle rests on the floor, a particle touching a vertical wall keeps falling along it
TEST(lagrange_wall_contact) {
	const double d = 0.8, rhoP = 20, g = 1e-4, nu = 0.02;
	const int n = 2500;
	CaseSpec s = box("lag_wall", 2, 2, 2, 2, nu);
	setGravity(s, -g);
	s.particles.push_back({ 8.3, 7.6, 1.5, d, rhoP });
	s.particles.push_back({ 0.38, 4.4, 13.0, d, rhoP });
	auto mesh = buildMesh(s);
	LagrangeParticles& P = *mesh->particles;
	runSteps(*mesh, n);
	double vt = terminalVelocity(d, rhoP, g, nu);
	report("particle on the floor: z - d/2", P.partZ[0] - d / 2);
	report("particle on the floor: vz / vt", P.veloZ[0] / vt);
	report("particle at the side wall: x - d/2", P.partX[1] - d / 2);
	report("particle at the side wall: fall / (vt n)", (13.0 - P.partZ[1]) / (vt * n));
	CHECK(P.partZ[0] > d / 2 - 0.02 && P.partZ[0] < d / 2 + 0.05);
	CHECK(P.partX[1] > d / 2 - 0.02 && P.partX[1] < d / 2 + 0.1);
	CHECK(P.partZ[1] > d / 2 - 0.02);
	CHECK(13.0 - P.partZ[1] > 0.5 * vt * n);
}

// free fall in the gas above the liquid, then the particle dives and settles
TEST(lagrange_free_fall_in_gas) {
	const double g = 1e-3;
	CaseSpec s = box("lag_gas", 2, 2, 2, 2, 0.02);
	setGravity(s, -g);
	s.flag = [](int, int, int k) { return k < 8 ? 1 : 2; };
	s.particles.push_back({ 8.3, 7.6, 13.5, 0.5, 5.0 });
	auto mesh = buildMesh(s);
	LagrangeParticles& P = *mesh->particles;
	double m0 = totalMass(*mesh, 0);
	runSteps(*mesh, 40);
	CHECK(P.inGas[0]);
	CHECK_NEAR(P.veloZ[0], -40 * g, 1e-12);
	CHECK_NEAR(P.partZ[0], 13.5 - g * 40 * 41 / 2.0, 1e-10);
	// nothing acts on the liquid yet
	P.applyForce();
	array<double, 3> f = spreadForce(*mesh);
	CHECK(P.Gz[0] == 0 && f[0] == 0 && f[1] == 0 && f[2] == 0);
	runSteps(*mesh, 3000, 40);
	double uMax = mesh->getUmax();
	report("z at the end", P.partZ[0]);
	report("uMax at the end", uMax);
	CHECK(!P.inGas[0]);
	CHECK(uMax == uMax && uMax < 0.05);
	CHECK(P.partZ[0] > 0.2 && P.partZ[0] < 7.5);
	CHECK_LT(fabs(totalMass(*mesh, 3040) - m0) / m0, 1e-9);
}

// gravity rotates around the z axis only on request (rotating mould), its z component is kept
TEST(particles_gravity_rotation) {
	for (int problemID : { 2, 4 }) {
		for (double omega : { 0.0, 0.01 }) {
			CaseSpec s = box("gravity_rotation", problemID, 2, 2, 1, 0.02);
			s.params["gx"] = "1e-4";
			s.params["gz"] = "-2e-4";
			if (omega != 0) {
				s.params["gravityRotationOmega"] = "0.01";
			}
			if (problemID == 2) {
				s.particles.push_back({ 8.3, 7.6, 4.1, 0.5, 2.0 });
			}
			else {
				setEulerParticles(s, 0.5, 2.0, 0.01);
				s.concentration = [](int, int, int) { return 1.0; };
			}
			auto mesh = buildMesh(s);
			const int n = 50;
			runSteps(*mesh, n);
			double a = omega * (n - 1); // Ct = 1, the last step has iter = n - 1
			Block* b = mesh->blocks[mesh->blockId[0]];
			CHECK_NEAR(b->G[0], 1e-4 * cos(a), 1e-15);
			CHECK_NEAR(b->G[1], 1e-4 * sin(a), 1e-15);
			CHECK_NEAR(b->G[2], -2e-4, 1e-15);
			if (problemID == 2) {
				CHECK_NEAR(mesh->particles->G[0], 1e-4 * cos(a), 1e-15);
				CHECK_NEAR(mesh->particles->G[1], 1e-4 * sin(a), 1e-15);
				CHECK_NEAR(mesh->particles->G[2], -2e-4, 1e-15);
			}
		}
	}
}

// a particle outside of the domain must not be mapped to another block
TEST(lagrange_particle_outside_of_domain) {
	CaseSpec s = box("lag_outside", 2, 2, 2, 1, 0.02);
	setGravity(s, -1e-4);
	s.particles.push_back({ -0.5, 12.0, 4.0, 0.5, 5.0 });
	auto mesh = buildMesh(s);
	LagrangeParticles& P = *mesh->particles;
	runSteps(*mesh, 3);
	P.applyForce();
	array<double, 3> f = spreadForce(*mesh);
	report("force applied to the fluid", f[2]);
	CHECK(f[0] == 0 && f[1] == 0 && f[2] == 0);
}

// ---------------------------------------------------------------------------------------------
// Eulerian particles, FVM
// ---------------------------------------------------------------------------------------------

// a layer of particles settles in a column: conservation, settling velocity and boundedness
TEST(euler_fvm_settling) {
	const double d = 1.0, rhoP = 20, g = 1e-4, nu = 0.02, c0 = 0.02;
	const int n = 1500;
	CaseSpec s = box("fvm_settling", 3, 1, 1, 6, nu);
	setGravity(s, -g);
	setEulerParticles(s, d, rhoP, 0);
	s.concentration = [](int, int, int k) { return (k >= 28 && k < 40) ? 0.02 : 0.0; };
	auto mesh = buildMesh(s);
	Moments m0 = moments<FvmBlock>(*mesh);
	runSteps(*mesh, n);
	Moments m1 = moments<FvmBlock>(*mesh);
	double vt = terminalVelocity(d, rhoP, g, nu);
	double shift = vt * (n - relaxationTime(d, rhoP, g, nu));
	report("settling distance (analytic)", shift);
	report("settling distance", m0.zc - m1.zc);
	report("relative change of the number of particles", (m1.sum - m0.sum) / m0.sum);
	report("c min / c0", m1.cMin / c0);
	report("c max / c0", m1.cMax / c0);
	CHECK_NEAR(m1.sum, m0.sum, 1e-12 * m0.sum);
	CHECK_NEAR(m0.zc - m1.zc, shift, 0.03 * shift);
	CHECK(m1.cMin > -1e-12 && m1.cMax < 1.02 * c0);
}

namespace {
	// a cloud of particles in a corner of a half filled box drives a flow; the result must not depend on the threads
	template <class B> double threadDifference(int problemID) {
		vector<double> ref;
		double worst = 0;
		for (int threads : { 1, 8 }) {
			CaseSpec s = box(problemID == 3 ? "fvm_threads" : "lbm_threads", problemID, 2, 2, 2, 0.02);
			setGravity(s, -1e-4);
			setEulerParticles(s, 1.0, 20, 0.01);
			s.flag = [](int, int, int k) { return k < 13 ? 1 : 2; };
			s.concentration = [](int i, int j, int k) { return (i < 9 && j > 5 && k >= 6 && k < 13) ? 0.05 : 0.0; };
			auto mesh = buildMesh(s);
			omp_set_num_threads(threads);
			runSteps(*mesh, 300);
			vector<double> v;
			for (int b = 0; b < mesh->nBlock; b++) {
				B* bl = (B*)mesh->blocks[mesh->blockId[b]];
				for (int ind = 0; ind < bl->nAll; ind += 7) {
					v.push_back(bl->concentration[ind]);
					v.push_back(bl->W[ind]);
				}
			}
			if (ref.empty()) {
				ref = v;
				double wMax = 0;
				for (size_t i = 1; i < v.size(); i += 2) {
					wMax = max(wMax, fabs(v[i]));
				}
				CHECK(wMax > 1e-5); // the particles really move the fluid
			}
			else {
				for (size_t i = 0; i < v.size(); i++) {
					worst = max(worst, fabs(v[i] - ref[i]));
				}
			}
		}
		omp_set_num_threads(4);
		return worst;
	}

	// sum of the force of the particles before and after it is spread to the liquid cells
	template <class B> array<double, 2> forceBalance(Mesh& mesh) {
		array<double, 2> f = { 0, 0 };
		for (int b = 0; b < mesh.nBlock; b++) {
			B* bl = (B*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				f[0] += bl->Gz[ind];
			}
		}
		for (int b = 0; b < mesh.nBlock; b++) {
			((B*)mesh.blocks[mesh.blockId[b]])->normalizeForce();
		}
		for (int b = 0; b < mesh.nBlock; b++) {
			B* bl = (B*)mesh.blocks[mesh.blockId[b]];
			bl->smoothForce();
			for (int ind = 0; ind < bl->nAll; ind++) {
				f[1] += bl->GzSmooth[ind] * bl->vof[ind]; // cellForce applies GzSmooth * vof / rho
			}
		}
		return f;
	}

	CaseSpec halfFilledBox(const string& name, int problemID) {
		CaseSpec s = box(name, problemID, 2, 2, 2, 0.02);
		setGravity(s, -1e-4);
		setEulerParticles(s, 0.5, 5, 0.01);
		s.flag = [](int, int, int k) { return k < 9 ? 1 : 2; };
		s.concentration = [](int i, int, int k) { return k < 9 ? 1.0 + 0.1 * i : 0.0; };
		return s;
	}
}

TEST(euler_thread_determinism) {
	double fvm = threadDifference<FvmBlock>(3);
	double lbm = threadDifference<LbmBlock>(4);
	report("FVM: largest difference 1 vs 8 threads", fvm);
	report("LBM: largest difference 1 vs 8 threads", lbm);
	CHECK_LT(fvm, 1e-12);
	CHECK_LT(lbm, 1e-12);
}

// the force of the particles reaches the fluid completely, also next to the walls and to the free surface
TEST(euler_force_spreading_is_conservative) {
	{
		auto mesh = buildMesh(halfFilledBox("fvm_force", 3));
		runSteps(*mesh, 20);
		mesh->evalBlocks(20, 30);
		for (int b = 0; b < mesh->nBlock; b++) {
			((FvmBlock*)mesh->blocks[mesh->blockId[b]])->computeNewParticleVelocityField();
		}
		array<double, 2> f = forceBalance<FvmBlock>(*mesh);
		report("FVM: spread force / force of the particles", f[1] / f[0]);
		CHECK(fabs(f[0]) > 1e-6);
		CHECK_NEAR(f[1], f[0], 1e-11 * fabs(f[0]));
	}
	{
		auto mesh = buildMesh(halfFilledBox("lbm_force", 4));
		runSteps(*mesh, 20);
		mesh->evalBlocks(20, 40);
		for (int b = 0; b < mesh->nBlock; b++) {
			Block* bl = mesh->blocks[mesh->blockId[b]];
			((LbmBlock*)bl)->computeNewParticleVelocityField(bl->wallDistance, bl->wallDistanceGradient);
		}
		array<double, 2> f = forceBalance<LbmBlock>(*mesh);
		report("LBM: spread force / force of the particles", f[1] / f[0]);
		CHECK(fabs(f[0]) > 1e-6);
		CHECK_NEAR(f[1], f[0], 1e-11 * fabs(f[0]));
	}
}

// ---------------------------------------------------------------------------------------------
// Eulerian particles, LBM
// ---------------------------------------------------------------------------------------------

TEST(euler_lbm_settling) {
	const double d = 1.0, rhoP = 20, g = 1e-4, nu = 0.02;
	const int n = 1500;
	CaseSpec s = box("lbm_settling", 4, 1, 1, 8, nu);
	setGravity(s, -g);
	setEulerParticles(s, d, rhoP, 0.005);
	s.concentration = [](int, int, int k) { return 0.02 * exp(-(k + 0.5 - 40) * (k + 0.5 - 40) / (2 * 3.0 * 3.0)); };
	auto mesh = buildMesh(s);
	Moments m0 = moments<LbmBlock>(*mesh);
	double g0 = scalarSum(*mesh);
	runSteps(*mesh, n);
	Moments m1 = moments<LbmBlock>(*mesh);
	double vt = terminalVelocity(d, rhoP, g, nu);
	double shift = vt * (n - relaxationTime(d, rhoP, g, nu));
	report("settling distance (analytic)", shift);
	report("settling distance", m0.zc - m1.zc);
	report("relative change of the number of particles", (m1.sum - m0.sum) / m0.sum);
	report("relative change of the sum of populations", (scalarSum(*mesh) - g0) / g0);
	report("c min / c0", m1.cMin / 0.02);
	CHECK_NEAR(scalarSum(*mesh), g0, 1e-12 * g0);
	CHECK_NEAR(m1.sum, m0.sum, 1e-6 * m0.sum);
	CHECK_NEAR(m0.zc - m1.zc, shift, 0.03 * shift);
	CHECK(m1.cMin > -1e-3 * 0.02); // the scalar LBM is not positivity preserving
}

// neutrally buoyant particles without gravity: the variance grows as 2 D t
TEST(euler_lbm_diffusion) {
	const double D = 0.02;
	const int n = 1000;
	CaseSpec s = box("lbm_diffusion", 4, 1, 1, 8, 0.02);
	setEulerParticles(s, 0.1, 1.0, D);
	s.concentration = [](int, int, int k) { return exp(-(k + 0.5 - 32) * (k + 0.5 - 32) / (2 * 3.0 * 3.0)); };
	auto mesh = buildMesh(s);
	Moments m0 = moments<LbmBlock>(*mesh);
	runSteps(*mesh, n);
	Moments m1 = moments<LbmBlock>(*mesh);
	double Deff = (m1.var - m0.var) / (2 * n);
	report("diffusivity (prescribed)", D);
	report("diffusivity (measured)", Deff);
	CHECK_NEAR(m1.sum, m0.sum, 1e-9 * m0.sum);
	CHECK_NEAR(m1.zc, m0.zc, 1e-6);
	CHECK_NEAR(Deff, D, 0.02 * D);
}

// the number of particles is conserved on every lattice
TEST(euler_lbm_lattices) {
	for (int model : { 319, 315, 29 }) {
		CaseSpec s = box("lbm_lattice", 4, 2, 2, model == 29 ? 1 : 2, 0.02);
		if (model == 29) {
			s.nk = 1;
		}
		s.params["model"] = to_string(model);
		setEulerParticles(s, 0.1, 1.0, 0.02);
		s.concentration = [](int i, int, int) { return 1.0 + 0.5 * sin(0.4 * i); };
		auto mesh = buildMesh(s);
		Moments m0 = moments<LbmBlock>(*mesh);
		runSteps(*mesh, 50);
		Moments m1 = moments<LbmBlock>(*mesh);
		report("model " + to_string(model) + ": relative change of the number of particles", (m1.sum - m0.sum) / m0.sum);
		CHECK_NEAR(m1.sum, m0.sum, 1e-9 * m0.sum);
	}
}

// gradient of the particle velocity used by the convective term
TEST(euler_lbm_velocity_gradient) {
	CaseSpec s = box("lbm_gradient", 4, 2, 2, 2, 0.02);
	setEulerParticles(s, 0.1, 1.0, 0.02);
	s.concentration = [](int, int, int) { return 1.0; };
	auto mesh = buildMesh(s);
	for (int b = 0; b < mesh->nBlock; b++) {
		LbmBlock* bl = (LbmBlock*)mesh->blocks[mesh->blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			int gi, gj, gk;
			bl->globalCoords(ind, gi, gj, gk);
			bl->UP[ind] = 0.01 * gi + 0.02 * gj;
			bl->VP[ind] = -0.03 * gk;
			bl->WP[ind] = 0.005 * gi;
		}
	}
	for (int b = 0; b < mesh->nBlock; b++) {
		((LbmBlock*)mesh->blocks[mesh->blockId[b]])->computeGradU();
	}
	// cell at a corner of four blocks, far from the walls
	int ind;
	LbmBlock* bl = (LbmBlock*)mesh->cellAt(8, 7, 8, ind);
	report("dUp/dx (exact 0.01)", bl->gUx[ind]);
	CHECK_NEAR(bl->gUx[ind], 0.01, 1e-12);
	CHECK_NEAR(bl->gUy[ind], 0.02, 1e-12);
	CHECK_NEAR(bl->gUz[ind], 0.0, 1e-12);
	CHECK_NEAR(bl->gVz[ind], -0.03, 1e-12);
	CHECK_NEAR(bl->gWx[ind], 0.005, 1e-12);
}

// the fluid velocity seen by the particles is defined in the interface cells too
TEST(euler_lbm_fluid_velocity_in_interface_cells) {
	CaseSpec s = box("lbm_interface", 4, 2, 2, 2, 0.02);
	setEulerParticles(s, 0.1, 2.0, 0.02);
	s.flag = [](int, int, int k) { return k < 8 ? 1 : 2; };
	s.concentration = [](int, int, int k) { return k < 8 ? 1.0 : 0.0; };
	auto mesh = buildMesh(s);
	for (int b = 0; b < mesh->nBlock; b++) {
		LbmBlock* bl = (LbmBlock*)mesh->blocks[mesh->blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			if (bl->flagMap[ind] == 1 || bl->flagMap[ind] == 3) {
				bl->U[ind] = 0.01;
			}
		}
	}
	for (int b = 0; b < mesh->nBlock; b++) {
		((LbmBlock*)mesh->blocks[mesh->blockId[b]])->interpolateFluidVelocity();
	}
	int ind;
	LbmBlock* bl = (LbmBlock*)mesh->cellAt(8, 7, 7, ind);
	CHECK(bl->flagMap[ind] == 3);
	report("fluid velocity in an interface cell (exact 0.01)", bl->UF[ind]);
	CHECK_NEAR(bl->UF[ind], 0.01, 1e-12);
	bl = (LbmBlock*)mesh->cellAt(8, 7, 4, ind);
	CHECK_NEAR(bl->UF[ind], 0.01, 1e-12);
}

// ---------------------------------------------------------------------------------------------
// both Eulerian models with a moving free surface (dam break in a closed box)
// ---------------------------------------------------------------------------------------------

namespace {
	template <class B> double liquidParticles(Mesh& mesh) {
		double s = 0;
		for (int b = 0; b < mesh.nBlock; b++) {
			B* bl = (B*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				if (bl->flagMap[ind] == 1 || bl->flagMap[ind] == 3) {
					s += bl->concentration[ind];
				}
			}
		}
		return s;
	}

	template <class B> double damBreakParticleLoss(int problemID) {
		CaseSpec s = box(problemID == 3 ? "fvm_dam" : "lbm_dam", problemID, 2, 1, 2, 0.05);
		s.params["gz"] = "-2e-4";
		setEulerParticles(s, 0.2, 2.0, 0.01);
		s.flag = [](int i, int, int k) { return (i < 7 && k < 11) ? 1 : 2; };
		s.concentration = [](int i, int, int k) { return (i < 7 && k < 11) ? 5.0 : 0.0; };
		auto mesh = buildMesh(s);
		double n0 = liquidParticles<B>(*mesh);
		runSteps(*mesh, 800);
		CHECK(mesh->getUmax() > 1e-3);
		return (liquidParticles<B>(*mesh) - n0) / n0;
	}
}

TEST(euler_particles_follow_free_surface) {
	double lossFvm = damBreakParticleLoss<FvmBlock>(3);
	double lossLbm = damBreakParticleLoss<LbmBlock>(4);
	report("FVM: relative change of the particles in the liquid", lossFvm);
	report("LBM: relative change of the particles in the liquid", lossLbm);
	CHECK_LT(fabs(lossFvm), 1e-10);
	CHECK_LT(fabs(lossLbm), 1e-4); // the scalar LBM conserves the populations, the concentration follows them
}

// ---------------------------------------------------------------------------------------------
// particles on a refined mesh: box 16^3 (2 x 2 x 2 blocks of the level 0), the half x < 8 is refined
// ---------------------------------------------------------------------------------------------

namespace {
	CaseSpec refinedBox(const string& name, int problemID, double nu) {
		CaseSpec s;
		s.name = name;
		s.ni = s.nj = s.nk = 8;
		s.octree = octreeBlocks(2, 2, 2, 1, [](int l, int i, int, int) { return l == 0 && i == 0; });
		s.params = latticeParams(problemID, 8, nu);
		s.params["maxParticlesInCell"] = "20";
		s.flagL = [](int, int, int, int) { return 1; };
		return s;
	}

	// number of the Eulerian particles in the computed cells (FVM: concentration, LBM: sum of the populations / theta)
	template <class B> array<double, 2> particleCount(Mesh& mesh, bool populations) {
		array<double, 2> n = { 0, 0 }; // total, finer level
		for (int b = 0; b < mesh.nBlock; b++) {
			B* bl = (B*)mesh.blocks[mesh.blockId[b]];
			for (int ind = 0; ind < bl->nAll; ind++) {
				if (bl->boundary[ind]) {
					continue;
				}
				double c = 0;
				if (populations) {
					for (int p = 0; p < bl->nf; p++) {
						c += ((BlockFreeSurfaceEulerParticlesLBM*)bl)->gc[bl->nf * ind + p];
					}
				}
				else {
					c = ((BlockFreeSurfaceEulerParticlesFVM*)bl)->concentration[ind];
				}
				n[0] += c;
				if (bl->level > 0) {
					n[1] += c;
				}
			}
		}
		return n;
	}

	template <class B> void refinedConservation(int problemID) {
		CaseSpec s = refinedBox(problemID == 3 ? "ref_fvm" : "ref_lbm", problemID, 0.02);
		s.params["gx"] = "-1e-4"; // the particles settle from the coarser half into the finer one
		setEulerParticles(s, 1.0, 20, 0.005);
		s.concentrationL = [](int l, int i, int j, int k) {
			return (l == 0 && i >= 9 && i < 15 && j >= 4 && j < 12 && k >= 4 && k < 12) ? 0.05 : 0.0;
		};
		auto mesh = buildMesh(s);
		bool pop = problemID == 4;
		array<double, 2> n0 = particleCount<B>(*mesh, pop);
		double worst = 0;
		for (int n = 0; n < 5; n++) {
			runSteps(*mesh, 100, 100 * n);
			array<double, 2> c = particleCount<B>(*mesh, pop);
			worst = max(worst, fabs(c[0] - n0[0]) / n0[0]);
		}
		array<double, 2> c = particleCount<B>(*mesh, pop);
		report(string(problemID == 3 ? "FVM" : "LBM") + ": relative change of the number of particles", worst);
		report(string(problemID == 3 ? "FVM" : "LBM") + ": part of the particles in the finer half", c[1] / c[0]);
		CHECK(mesh->getUmax() == mesh->getUmax());
		CHECK(c[1] / c[0] > 0.2); // the particles crossed the level interface
		CHECK_LT(worst, 1e-11);
	}
}

TEST(refinement_euler_particles_conservation) {
	refinedConservation<FvmBlock>(3);
	refinedConservation<LbmBlock>(4);
}

// the impulse of the particles reaches the fluid in both levels in the right units
TEST(refinement_lagrange_force_and_settling) {
	const double d = 0.5, rhoP = 20, g = 1e-4, nu = 0.02;
	CaseSpec s = refinedBox("ref_lag", 2, nu);
	s.params["gz"] = to_string(-g);
	s.params["initializeWithHydroPress"] = "1";
	s.particles.push_back({ 4.3, 7.6, 12.4, d, rhoP });   // finer half
	s.particles.push_back({ 12.3, 7.6, 12.4, d, rhoP });  // coarser half
	auto mesh = buildMesh(s);
	LagrangeParticles& P = *mesh->particles;
	runSteps(*mesh, 600);
	double vt = terminalVelocity(d, rhoP, g, nu);
	for (int p = 0; p < 2; p++) {
		double slip = P.veloFZ[p] - P.veloZ[p];
		report(string(p == 0 ? "finer" : "coarser") + " half: slip / vt", slip / vt);
		report(string(p == 0 ? "finer" : "coarser") + " half: fluid velocity at the particle / vt", P.veloFZ[p] / vt);
		CHECK_NEAR(slip, vt, 0.02 * vt);
		CHECK(P.veloFZ[p] < 0 && P.veloFZ[p] > -0.5 * vt); // the fluid is dragged down, but not more than in a uniform mesh
	}
	// momentum given to the fluid during one step of the level 0 (sum over the levels in the units of the level 0)
	P.applyForce();
	double fz = 0;
	for (int b = 0; b < mesh->nBlock; b++) {
		LagBlock* bl = (LagBlock*)mesh->blocks[mesh->blockId[b]];
		double scale = (double)(1 << bl->level) / (1 << (3 * bl->level)); // 2^L steps, momentum unit 8^-L
		for (int ind = 0; ind < bl->nAll; ind++) {
			fz += bl->particleForceZ[ind] * bl->vof[ind] * scale;
		}
	}
	report("momentum given to the fluid / impulse of the particles", fz / (P.Gz[0] + P.Gz[1]));
	CHECK_NEAR(fz, P.Gz[0] + P.Gz[1], 1e-10 * fabs(P.Gz[0] + P.Gz[1]));
	// wall distance across the level interface: the exact distance is 7.7, the interpolated cell centred distance is
	// about half a cell shorter (7.23 on a uniform mesh of the level 0), more accurate on the finer side
	double La = P.wallDistance(0, 7.95, 8.3, 8.3), Lb = P.wallDistance(0, 8.05, 8.3, 8.3);
	report("wall distance on the finer side of the interface (exact 7.7)", La);
	report("wall distance on the coarser side of the interface (exact 7.7)", Lb);
	CHECK_NEAR(La, 7.7, 0.6);
	CHECK_NEAR(Lb, 7.7, 0.6);
	CHECK_NEAR(La, Lb, 0.5);
}
