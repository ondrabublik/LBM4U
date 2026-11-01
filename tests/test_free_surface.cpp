// Tests of the free surface solver (problemID 1)

#include "TestFramework.h"
#include "TestCase.h"
#include <omp.h>

using namespace std;

namespace {
	// dam break in a closed box 16 x 8 x 16 split into 2 x 1 x 2 blocks
	CaseSpec damBreak(const string& name, int collisionModel) {
		CaseSpec s;
		s.name = name;
		s.nBi = 2;
		s.nBk = 2;
		s.ni = s.nj = s.nk = 8;
		s.params = latticeParams(1, 8, 0.05);
		s.params["collisionModel"] = to_string(collisionModel);
		s.params["gz"] = "-2e-4";
		s.flag = [](int i, int j, int k) { return (i < 7 && k < 11) ? 1 : 2; };
		return s;
	}

	double totalMass(Mesh& m, int iter) {
		m.evalBlocks(iter, 105);
		return m.massTot;
	}

	double maxVelocity(Mesh& m) {
		return m.getUmax();
	}
}

TEST(free_surface_mass_conservation) {
	for (int cm : { COLLISION_TRT, COLLISION_MRT, COLLISION_BGK }) {
		auto mesh = buildMesh(damBreak("fs_mass", cm));
		double m0 = totalMass(*mesh, 0);
		for (int n = 0; n < 6; n++) {
			runSteps(*mesh, 100, n * 100);
			double m = totalMass(*mesh, (n + 1) * 100);
			CHECK_LT(fabs(m - m0) / m0, POP_TOL);
		}
		CHECK(maxVelocity(*mesh) > 1e-3); // the column really collapses
	}
}

TEST(free_surface_thread_determinism) {
	vector<double> ref;
	for (int threads : { 1, 8 }) {
		auto mesh = buildMesh(damBreak("fs_threads", COLLISION_TRT));
		omp_set_num_threads(threads);
		runSteps(*mesh, 300);
		vector<double> v;
		for (int i = 0; i < 16; i++) {
			for (int k = 0; k < 16; k++) {
				int ind;
				Block* b = mesh->cellAt(i, 4, k, ind);
				v.push_back(((BlockFreeSurface*)b)->mass[ind]);
				v.push_back(velocityAt(*mesh, i, 4, k)[0]);
			}
		}
		if (ref.empty()) {
			ref = v;
		}
		else {
			for (size_t i = 0; i < v.size(); i++) {
				CHECK_NEAR(v[i], ref[i], 1e-12);
			}
		}
	}
	omp_set_num_threads(4);
}

TEST(hydrostatic_initialisation_is_at_rest) {
	double u[2];
	for (int hydro = 0; hydro < 2; hydro++) {
		CaseSpec s;
		s.name = "hydrostatic";
		s.ni = s.nj = 8;
		s.nk = 16;
		s.params = latticeParams(1, 8, 0.05);
		s.params["gz"] = "-1e-4";
		s.params["initializeWithHydroPress"] = to_string(hydro);
		s.flag = [](int i, int j, int k) { return k < 10 ? 1 : 2; };
		auto mesh = buildMesh(s);
		u[hydro] = 0;
		for (int n = 0; n < 300; n += 10) {
			runSteps(*mesh, 10, n);
			u[hydro] = max(u[hydro], maxVelocity(*mesh));
		}
	}
	// without the hydrostatic initialisation the column falls with g until the pressure builds up
	CHECK_LT(u[1], 2e-4);
	CHECK_LT(u[1], 0.2 * u[0]);
}

namespace {
	// sphere of radius R: volume fractions sampled by 24^3 sub points
	unique_ptr<Mesh> sphereMesh(const string& name, int n, double R, double sigma, int curvatureModel) {
		CaseSpec s;
		s.name = name;
		s.ni = s.nj = s.nk = n;
		s.params = latticeParams(1, n, 1.0 / 6);
		s.params["surfaceTension"] = to_string(sigma);
		s.params["curvatureModel"] = to_string(curvatureModel);
		double c = n / 2.0 + 0.3;
		auto frac = [=](int i, int j, int k) {
			double dx = i + 0.5 - c, dy = j + 0.5 - c, dz = k + 0.5 - c;
			double d = sqrt(dx * dx + dy * dy + dz * dz) - R;
			if (d > 0.9) {
				return 0.0;
			}
			if (d < -0.9) {
				return 1.0;
			}
			int in = 0, m = 24;
			for (int a = 0; a < m; a++)
				for (int b = 0; b < m; b++)
					for (int d = 0; d < m; d++) {
						double x = i + (a + 0.5) / m - c, y = j + (b + 0.5) / m - c, z = k + (d + 0.5) / m - c;
						in += (x * x + y * y + z * z < R * R);
					}
			return in / double(m * m * m);
		};
		s.flag = [=](int i, int j, int k) { return frac(i, j, k) > 0 ? 1 : 2; };
		auto mesh = buildMesh(s);
		// every partially filled cell is an interface cell with the exact volume fraction
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < n; k++) {
					int ind;
					BlockFreeSurface* b = (BlockFreeSurface*)mesh->cellAt(i, j, k, ind);
					double v = frac(i, j, k);
					if (b->flagMap[ind] == 3 || (v > 0 && v < 1)) {
						b->flagMap[ind] = 3;
						b->vof[ind] = v;
						b->mass[ind] = v * b->Rho[ind];
					}
				}
		return mesh;
	}

	// mean and rms relative error of the stored curvature (-1/R expected)
	void curvatureError(Mesh& m, double R, double& meanErr, double& rmsErr) {
		for (int level = 1; level <= 4; level++) {
			m.forBlocks([level](Block* b) { ((BlockFreeSurface*)b)->findSurfaceLevels(level); });
		}
		m.evalBlocks(0, 5);
		m.evalBlocks(0, 6);
		m.evalBlocks(0, 7);
		double sum = 0, sum2 = 0;
		int cnt = 0;
		BlockFreeSurface* b = (BlockFreeSurface*)m.blocks[m.blockId[0]];
		for (int ind = 0; ind < b->nAll; ind++) {
			if (b->flagMap[ind] == 3 && b->vof[ind] > 0.05 && b->vof[ind] < 0.95) {
				double e = (b->curvature[ind] - (-1.0 / R)) * R;
				sum += e;
				sum2 += e * e;
				cnt++;
			}
		}
		meanErr = sum / cnt;
		rmsErr = sqrt(sum2 / cnt);
	}
}

TEST(curvature_height_function_sphere) {
	double R = 8;
	double meanHF, rmsHF, meanDiv, rmsDiv;
	auto hf = sphereMesh("curv_hf", 32, R, 1e-3, 1);
	curvatureError(*hf, R, meanHF, rmsHF);
	auto dv = sphereMesh("curv_div", 32, R, 1e-3, 0);
	curvatureError(*dv, R, meanDiv, rmsDiv);
	CHECK_LT(fabs(meanHF), 0.05);
	CHECK_LT(rmsHF, 0.15);
	CHECK_LT(rmsHF, rmsDiv);
}

TEST(laplace_pressure_of_drop) {
	double R = 7, sigma = 2e-3;
	auto mesh = sphereMesh("laplace", 28, R, sigma, 0);
	runSteps(*mesh, 1500);
	CHECK_LT(mesh->getUmax(), 1e-3); // spurious currents
	// density in the centre of the drop
	int ind;
	Block* b = mesh->cellAt(14, 14, 14, ind);
	double dRho = b->Rho[ind] - 1.0;
	CHECK_LT(fabs(dRho - 6 * sigma / R) / (6 * sigma / R), 0.05);
}

TEST(non_newtonian_local_shear_rate) {
	// Newtonian power law (n = 1) in the gravity driven channel, local shear rate vs analytic |du/dz|
	CaseSpec s;
	s.name = "shear_rate";
	s.ni = s.nj = 4;
	s.nk = 20;
	double nu = 0.1, g = 1e-5;
	s.params = latticeParams(1, 4, nu);
	s.params["xPeriodicity"] = "1";
	s.params["yPeriodicity"] = "1";
	s.params["gx"] = to_string(g);
	s.params["viscosityModel"] = "1";
	s.params["fluidModelParameters"] = "0.1, 1";
	s.flag = [](int, int, int) { return 1; };
	auto mesh = buildMesh(s);
	runSteps(*mesh, 15000);
	double maxRate = g / nu * 10;
	for (int k = 2; k < 18; k++) {
		int ind;
		BlockFreeSurface* b = (BlockFreeSurface*)mesh->cellAt(1, 1, k, ind);
		double expected = g / nu * fabs(10 - (k + 0.5));
		CHECK_NEAR(b->shearRate[ind], expected, 0.03 * maxRate);
		CHECK_NEAR(b->viscosity[ind], nu, 1e-12);
	}
}

TEST(free_surface_inlet_filling) {
	CaseSpec s;
	s.name = "fs_inlet";
	s.ni = s.nj = 8;
	s.nk = 16;
	s.params = latticeParams(1, 8, 0.05);
	s.params["gz"] = "-1e-5";
	s.params["wInlet"] = "0.01";
	s.cellType = [](int i, int j, int k) { return k == 0 ? CT_INLET : CT_FLUID; };
	s.flag = [](int i, int j, int k) { return k == 1 ? 1 : 2; };
	auto mesh = buildMesh(s);
	double m0 = totalMass(*mesh, 0);
	runSteps(*mesh, 800);
	double dm = totalMass(*mesh, 800) - m0;
	// each fluid cell above the inlet receives 6 w rho u per inlet neighbour, edge cells miss the diagonals
	// pointing outside of the domain: 64 cells - 32 missing diagonals with 6 w = 1/6
	double expected = 0.01 * (64 - 32.0 / 6) * 800;
	CHECK_NEAR(dm, expected, 0.05 * expected);
}

TEST(free_surface_d2q9) {
	CaseSpec s;
	s.name = "fs_2d";
	s.nBi = 2;
	s.ni = 16;
	s.nj = 16;
	s.nk = 1;
	s.params = latticeParams(1, 16, 0.05);
	s.params["model"] = "29";
	s.params["gy"] = "-2e-4";
	s.flag = [](int i, int j, int k) { return (i < 10 && j < 10) ? 1 : 2; };
	auto mesh = buildMesh(s);
	double m0 = totalMass(*mesh, 0);
	runSteps(*mesh, 400);
	double m = totalMass(*mesh, 400);
	CHECK_LT(fabs(m - m0) / m0, POP_TOL);
	CHECK(maxVelocity(*mesh) > 1e-3);
}
