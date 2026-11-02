// Tests of the static octree grid refinement

#include "TestFramework.h"
#include "TestCase.h"
#include <omp.h>

using namespace std;

namespace {
	array<double, 3> velocityAtPoint(Mesh& mesh, double x, double y, double z) {
		int ind;
		Block* b = blockAtPoint(mesh, x, y, z, ind);
		array<double, 3> u = { 0, 0, 0 };
		if (b && !b->boundary[ind]) {
			b->outputVelocity(ind, u.data());
			for (double& c : u) {
				c *= b->Cu;
			}
		}
		return u;
	}

	// channel periodic in x and y with walls at the domain boundary in z (height 32 cells of the level 0),
	// driven by gravity in x; refine(level, bi, bj, bk) selects the refined blocks
	CaseSpec channel(const string& name, int nBi, function<bool(int, int, int, int)> refine, double g, double nu) {
		CaseSpec s;
		s.name = name;
		s.ni = s.nj = s.nk = 8;
		s.octree = octreeBlocks(nBi, 1, 4, 1, refine);
		s.params = latticeParams(0, 8, nu);
		s.params["xPeriodicity"] = "1";
		s.params["yPeriodicity"] = "1";
		s.params["gx"] = to_string(g);
		return s;
	}

	// maximum error of the velocity profile relative to the maximum velocity, sampled at the finest cells
	double poiseuilleError(Mesh& mesh, double x, double g, double nu) {
		double H = 32, maxErr = 0, uMaxA = g / (8 * nu) * H * H;
		for (int k = 0; k < 64; k++) {
			double z = (k + 0.5) / 2;
			int ind;
			Block* b = blockAtPoint(mesh, x, 1.0, z, ind);
			int gi, gj, gk;
			b->globalCoords(ind, gi, gj, gk);
			double zc = (gk + 0.5) / (1 << b->level);
			double ua = g / (2 * nu) * zc * (H - zc);
			maxErr = max(maxErr, fabs(velocityAtPoint(mesh, x, 1.0, z)[0] - ua));
		}
		return maxErr / uMaxA;
	}
}

TEST(refinement_rest_state_is_kept) {
	CaseSpec s;
	s.name = "ref_rest";
	s.ni = s.nj = s.nk = 8;
	s.octree = octreeBlocks(2, 2, 2, 2, [](int l, int i, int j, int k) { return (l == 0 && i == 0 && j == 0 && k == 0) || (l == 1 && i == 0 && j == 0 && k == 0); });
	s.params = latticeParams(0, 8, 0.1);
	auto mesh = buildMesh(s);
	CHECK(mesh->maxLevel == 2);
	double m0 = latticeMass(*mesh);
	runSteps(*mesh, 50);
	CHECK(mesh->getUmax() < 1e-12);
	CHECK_NEAR(latticeMass(*mesh), m0, 1e-12 * m0);
}

TEST(refinement_unbalanced_mesh_is_rejected) {
	CaseSpec s;
	s.name = "ref_unbalanced";
	s.ni = s.nj = s.nk = 8;
	// level 2 blocks next to a level 0 block
	s.octree = octreeBlocks(2, 1, 1, 2, [](int l, int i, int j, int k) { return (l == 0 && i == 1) || (l == 1 && i == 2); });
	s.params = latticeParams(0, 8, 0.1);
	CHECK_THROWS(buildMesh(s));
}

TEST(refinement_mass_conservation) {
	// flow through a refined segment of a channel, initial shear wave in z
	CaseSpec s = channel("ref_mass", 4, [](int l, int i, int j, int k) { return i == 1 || i == 2; }, 1e-5, 0.05);
	s.initL = [](int l, int gi, int gj, int gk) {
		double z = (gk + 0.5) / (1 << l);
		return array<double, 4>{ 1.0 + 0.01 * sin(z / 5), 0.02 * sin(6.283185307179586 * z / 32), 0.0, 0.01 * cos(gi * 0.3 / (1 << l)) };
	};
	auto mesh = buildMesh(s);
	double m0 = latticeMass(*mesh);
	runSteps(*mesh, 200);
	CHECK(mesh->getUmax() > 1e-3);
	CHECK_NEAR(latticeMass(*mesh), m0, POP_TOL * m0);
}

TEST(refinement_poiseuille_refined_walls) {
	// interfaces parallel to the flow: the blocks next to the walls are refined
	double g = 2e-6, nu = 0.2;
	auto mesh = buildMesh(channel("ref_pois_walls", 1, [](int l, int i, int j, int k) { return k == 0 || k == 3; }, g, nu));
	runSteps(*mesh, 6000);
	CHECK_LT(poiseuilleError(*mesh, 1.0, g, nu), 2e-3);
}

TEST(refinement_poiseuille_refined_segment) {
	// interfaces normal to the flow: a segment of the channel is refined
	double g = 2e-6, nu = 0.2;
	auto mesh = buildMesh(channel("ref_pois_segment", 4, [](int l, int i, int j, int k) { return i == 1 || i == 2; }, g, nu));
	runSteps(*mesh, 6000);
	for (double x : { 4.0, 12.0, 16.0, 20.0, 28.0 }) {
		CHECK_LT(poiseuilleError(*mesh, x, g, nu), 2.5e-2); // walls crossing the interfaces
	}
}

namespace {
	// dam break in a closed box of 2 x 1 x 2 blocks (16 x 8 x 16 cells of the level 0), refine selects the refined blocks
	CaseSpec damBreakRefined(const string& name, function<bool(int, int, int, int)> refine) {
		CaseSpec s;
		s.name = name;
		s.ni = s.nj = s.nk = 8;
		s.octree = octreeBlocks(2, 1, 2, 1, refine);
		s.params = latticeParams(1, 8, 0.05);
		s.params["gz"] = "-2e-4";
		// water column x < 7, z < 11 (cells of the level 0)
		s.flagL = [](int l, int gi, int gj, int gk) {
			double x = (gi + 0.5) / (1 << l), z = (gk + 0.5) / (1 << l);
			return (x < 7 && z < 11) ? 1 : 2;
		};
		return s;
	}

	double fsMass(Mesh& m, int iter) {
		m.evalBlocks(iter, 105);
		return m.massTot;
	}
}

TEST(refinement_free_surface_mass_conservation) {
	// the column and its surface cross the interface of the levels
	for (int variant = 0; variant < 2; variant++) {
		auto refine = variant == 0 ? function<bool(int, int, int, int)>([](int l, int i, int j, int k) { return i == 0 && k == 0; })
			: function<bool(int, int, int, int)>([](int l, int i, int j, int k) { return i == 1 || k == 1; });
		auto mesh = buildMesh(damBreakRefined("ref_fs_mass" + to_string(variant), refine));
		double m0 = fsMass(*mesh, 0);
		for (int n = 0; n < 6; n++) {
			runSteps(*mesh, 100, n * 100);
			double m = fsMass(*mesh, (n + 1) * 100);
			CHECK_LT(fabs(m - m0) / m0, 1e-9);
		}
		double u = mesh->getUmax();
		CHECK(u > 1e-3 && u < 0.3);
	}
}

TEST(refinement_free_surface_at_rest) {
	// water layer with the surface in the refined blocks stays at rest
	CaseSpec s = damBreakRefined("ref_fs_rest", [](int l, int i, int j, int k) { return k == 1; });
	s.params["initializeWithHydroPress"] = "1";
	s.flagL = [](int l, int gi, int gj, int gk) { return (gk + 0.5) / (1 << l) < 11 ? 1 : 2; };
	auto mesh = buildMesh(s);
	double m0 = fsMass(*mesh, 0);
	runSteps(*mesh, 400);
	CHECK_LT(mesh->getUmax(), 2e-3);
	CHECK_LT(fabs(fsMass(*mesh, 400) - m0) / m0, 1e-9);
}
