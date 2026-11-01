// Tests of the single phase solver (problemID 0)

#include "TestFramework.h"
#include "TestCase.h"
#include <omp.h>

using namespace std;

namespace {
	// periodic channel in x and y with walls at the domain boundary in z, driven by gravity in x
	CaseSpec poiseuilleCase(const string& name, int collisionModel, int nBk, double g) {
		CaseSpec s;
		s.name = name;
		s.ni = 4;
		s.nj = 4;
		s.nk = 20 / nBk;
		s.nBk = nBk;
		s.params = latticeParams(0, 4, 0.1);
		s.params["collisionModel"] = to_string(collisionModel);
		s.params["xPeriodicity"] = "1";
		s.params["yPeriodicity"] = "1";
		s.params["gx"] = to_string(g);
		return s;
	}

	double poiseuilleError(Mesh& mesh, double g, double nu, int H) {
		double maxErr = 0, uMaxA = 0;
		for (int k = 0; k < H; k++) {
			double z = k + 0.5;
			double ua = g / (2 * nu) * z * (H - z);
			uMaxA = max(uMaxA, ua);
			maxErr = max(maxErr, fabs(velocityAt(mesh, 1, 2, k)[0] - ua));
		}
		return maxErr / uMaxA;
	}
}

TEST(guo_forcing_momentum) {
	for (int cm : { COLLISION_BGK, COLLISION_TRT, COLLISION_MRT }) {
		CaseSpec s;
		s.name = "forcing";
		s.ni = s.nj = s.nk = 4;
		s.params = latticeParams(0, 4, 0.1);
		s.params["collisionModel"] = to_string(cm);
		s.params["xPeriodicity"] = "1";
		s.params["yPeriodicity"] = "1";
		s.params["zPeriodicity"] = "1";
		s.params["gx"] = "1e-5";
		auto mesh = buildMesh(s);
		int n = 100;
		runSteps(*mesh, n);
		// velocity of the n-th collision includes half of the force: (n - 1/2) g
		for (int i = 0; i < 4; i++) {
			CHECK_NEAR(velocityAt(*mesh, i, i, 3 - i)[0], (n - 0.5) * 1e-5, (n - 0.5) * 1e-5 * POP_TOL);
		}
	}
}

TEST(poiseuille_profile) {
	double g = 1e-5, nu = 0.1;
	for (int cm : { COLLISION_TRT, COLLISION_MRT, COLLISION_BGK }) {
		auto mesh = buildMesh(poiseuilleCase("poiseuille", cm, 1, g));
		runSteps(*mesh, 20000);
		double err = poiseuilleError(*mesh, g, nu, 20);
		// TRT with magic 3/16 and MRT with the related energy flux rates place the wall exactly halfway
		CHECK_LT(err, (cm == COLLISION_BGK ? 0.02 : 1e-4));
	}
}

TEST(multiblock_equals_single_block) {
	auto one = buildMesh(poiseuilleCase("pois_one", COLLISION_TRT, 1, 1e-5));
	auto four = buildMesh(poiseuilleCase("pois_four", COLLISION_TRT, 4, 1e-5));
	runSteps(*one, 2000);
	runSteps(*four, 2000);
	for (int k = 0; k < 20; k++) {
		for (int j = 0; j < 4; j++) {
			auto a = velocityAt(*one, 2, j, k);
			auto b = velocityAt(*four, 2, j, k);
			CHECK_NEAR(a[0], b[0], 1e-16);
			CHECK_NEAR(a[2], b[2], 1e-16);
		}
	}
}

TEST(wall_distance) {
	CaseSpec s = poiseuilleCase("distance", COLLISION_TRT, 2, 0);
	s.nk = 5; // 10 cells in z
	auto mesh = buildMesh(s);
	for (int k = 0; k < 10; k++) {
		int ind;
		Block* b = mesh->cellAt(1, 1, k, ind);
		CHECK_NEAR(b->wallDistance[ind], min(k + 0.5, 10 - k - 0.5), 1e-12);
	}
}

TEST(inlet_outlet_mass_flux) {
	CaseSpec s;
	s.name = "inlet_outlet";
	s.nBi = 4;
	s.ni = s.nj = s.nk = 8;
	s.params = latticeParams(0, 8, 0.05);
	s.params["uInlet"] = "0.02";
	s.cellType = [](int i, int j, int k) { return i == 0 ? CT_INLET : (i == 31 ? CT_OUTLET : CT_FLUID); };
	auto mesh = buildMesh(s);
	runSteps(*mesh, 6000);
	auto flux = [&](int i) {
		double q = 0;
		for (int j = 0; j < 8; j++)
			for (int k = 0; k < 8; k++) {
				int ind;
				Block* b = mesh->cellAt(i, j, k, ind);
				q += b->Rho[ind] * velocityAt(*mesh, i, j, k)[0];
			}
		return q;
	};
	double q1 = flux(8), q2 = flux(24);
	CHECK_LT(fabs(q1 - q2) / q1, 0.01);
	CHECK_LT(fabs(q1 - 0.02 * 64) / (0.02 * 64), 0.1);
}

TEST(moving_lid_drives_flow) {
	CaseSpec s;
	s.name = "cavity";
	s.ni = s.nj = s.nk = 16;
	s.params = latticeParams(0, 16, 0.05);
	s.params["movingLidTop"] = "1";
	s.params["uwall"] = "0.05";
	auto mesh = buildMesh(s);
	runSteps(*mesh, 2000);
	double uTop = velocityAt(*mesh, 8, 8, 15)[0];
	double uBottom = velocityAt(*mesh, 8, 8, 2)[0];
	CHECK(uTop > 0.01);
	CHECK(uBottom < 0);
}

TEST(thread_count_does_not_change_result) {
	vector<double> ref;
	for (int threads : { 1, 6 }) {
		omp_set_num_threads(threads);
		CaseSpec s = poiseuilleCase("threads_sp", COLLISION_MRT, 4, 1e-5);
		s.nBi = 2;
		s.nBj = 2;
		auto mesh = buildMesh(s);
		omp_set_num_threads(threads);
		runSteps(*mesh, 500);
		vector<double> u;
		for (int k = 0; k < 20; k++) {
			u.push_back(velocityAt(*mesh, 3, 5, k)[0]);
		}
		if (ref.empty()) {
			ref = u;
		}
		else {
			for (size_t i = 0; i < u.size(); i++) {
				CHECK(u[i] == ref[i]);
			}
		}
	}
	omp_set_num_threads(4);
}

TEST(smagorinsky_increases_viscosity) {
	double uLam, uLes;
	for (int les = 0; les < 2; les++) {
		CaseSpec s = poiseuilleCase("les", COLLISION_TRT, 1, 1e-4);
		s.params["vKin"] = "0.01";
		s.params["smagorinsky"] = les ? "0.5" : "0";
		auto mesh = buildMesh(s);
		runSteps(*mesh, 4000);
		(les ? uLes : uLam) = velocityAt(*mesh, 1, 1, 10)[0];
	}
	CHECK_LT(uLes, 0.95 * uLam);
	CHECK_LT(0.0, uLes);
}
