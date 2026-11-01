// Unit tests of the lattice, collision operators, parameters and neighbour tables

#include "TestFramework.h"
#include "TestCase.h"
#include "../Collision.h"
#include "../SimulationIO.h"
#include <random>
#include <fstream>
#include <filesystem>

using namespace std;

namespace {
	struct Lattice {
		int nf;
		int* ex, * ey, * ez, * opp;
		double* w;
	};

	Lattice lattice(int model) {
		// Mesh::model fills the arrays, a tiny mesh is the simplest way to get them
		CaseSpec s;
		s.name = "lattice" + to_string(model);
		s.ni = s.nj = 4;
		s.nk = (model == 29) ? 1 : 4;
		s.params = latticeParams(0, 4, 0.1);
		s.params["model"] = to_string(model);
		static vector<unique_ptr<Mesh> > keep;
		keep.push_back(buildMesh(s));
		Mesh& m = *keep.back();
		return Lattice{ m.nf, m.ex, m.ey, m.ez, m.opp, m.w };
	}
}

TEST(lattice_isotropy) {
	for (int model : { 29, 315, 319 }) {
		Lattice L = lattice(model);
		double s = 0, sx = 0, sxx = 0, syy = 0, sxy = 0;
		for (int p = 0; p < L.nf; p++) {
			s += L.w[p];
			sx += L.w[p] * L.ex[p];
			sxx += L.w[p] * L.ex[p] * L.ex[p];
			syy += L.w[p] * L.ey[p] * L.ey[p];
			sxy += L.w[p] * L.ex[p] * L.ey[p];
			CHECK(L.ex[L.opp[p]] == -L.ex[p]);
			CHECK(L.ey[L.opp[p]] == -L.ey[p]);
			CHECK(L.ez[L.opp[p]] == -L.ez[p]);
		}
		CHECK_NEAR(s, 1.0, 1e-14);
		CHECK_NEAR(sx, 0.0, 1e-14);
		CHECK_NEAR(sxx, 1.0 / 3, 1e-14);
		CHECK_NEAR(syy, 1.0 / 3, 1e-14);
		CHECK_NEAR(sxy, 0.0, 1e-14);
	}
}

TEST(mrt_matrix_inverse) {
	for (int model : { 29, 315, 319 }) {
		Lattice L = lattice(model);
		CollisionOperator c;
		c.init(COLLISION_MRT, L.nf, L.ex, L.ey, L.ez, L.opp, 3.0 / 16, true, 0.8);
		int n = L.nf;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double sum = 0;
				for (int k = 0; k < n; k++) {
					sum += c.M[n * i + k] * c.iM[n * k + j];
				}
				CHECK_NEAR(sum, i == j ? 1.0 : 0.0, 1e-12);
			}
		}
		// rows of M are orthogonal
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				double sum = 0;
				for (int k = 0; k < n; k++) {
					sum += c.M[n * i + k] * c.M[n * j + k];
				}
				CHECK_NEAR(sum, 0.0, 1e-12);
			}
		}
	}
}

TEST(collision_conserves_mass_and_momentum) {
	mt19937 gen(1);
	uniform_real_distribution<double> dist(0.0, 0.1);
	for (int model : { 315, 319 }) {
		Lattice L = lattice(model);
		for (int cm : { COLLISION_BGK, COLLISION_TRT, COLLISION_MRT }) {
			CollisionOperator c;
			c.init(cm, L.nf, L.ex, L.ey, L.ez, L.opp, 3.0 / 16, true, 0.7);
			for (double tau : { 0.7, 0.9, 1.7 }) {
				double f[27], feq[27], fneq[27], delta[27];
				for (int p = 0; p < L.nf; p++) {
					f[p] = L.w[p] + dist(gen) * L.w[p];
				}
				double rho = 0, j[3] = { 0, 0, 0 };
				for (int p = 0; p < L.nf; p++) {
					rho += f[p];
					j[0] += L.ex[p] * f[p];
					j[1] += L.ey[p] * f[p];
					j[2] += L.ez[p] * f[p];
				}
				double u[3] = { j[0] / rho, j[1] / rho, j[2] / rho };
				double u2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
				for (int p = 0; p < L.nf; p++) {
					double eu = L.ex[p] * u[0] + L.ey[p] * u[1] + L.ez[p] * u[2];
					feq[p] = L.w[p] * rho * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * u2);
					fneq[p] = f[p] - feq[p];
				}
				c.relax(fneq, tau, delta);
				double dm = 0, dj[3] = { 0, 0, 0 };
				for (int p = 0; p < L.nf; p++) {
					dm += delta[p];
					dj[0] += L.ex[p] * delta[p];
					dj[1] += L.ey[p] * delta[p];
					dj[2] += L.ez[p] * delta[p];
				}
				CHECK_NEAR(dm, 0.0, 1e-14);
				CHECK_NEAR(dj[0], 0.0, 1e-14);
				CHECK_NEAR(dj[1], 0.0, 1e-14);
				CHECK_NEAR(dj[2], 0.0, 1e-14);
			}
		}
	}
}

// all operators relax the deviatoric stress with 1/tau
TEST(collision_stress_relaxation_rate) {
	Lattice L = lattice(319);
	for (int cm : { COLLISION_BGK, COLLISION_TRT, COLLISION_MRT }) {
		CollisionOperator c;
		c.init(cm, L.nf, L.ex, L.ey, L.ez, L.opp, 3.0 / 16, true, 0.6);
		for (double tau : { 0.6, 1.3 }) {
			double fneq[27], delta[27];
			for (int p = 0; p < L.nf; p++) {
				fneq[p] = L.w[p] * L.ex[p] * L.ey[p]; // pure xy shear moment
			}
			c.relax(fneq, tau, delta);
			double pxy = 0, dpxy = 0;
			for (int p = 0; p < L.nf; p++) {
				pxy += L.ex[p] * L.ey[p] * fneq[p];
				dpxy += L.ex[p] * L.ey[p] * delta[p];
			}
			CHECK_NEAR(dpxy, pxy / tau, 1e-14);
		}
	}
}

TEST(parameter_parsing) {
	ParamMap p = parseParams("a = 1\n  b=2.5e-3 # comment\n% c = 7\nmaxIter = 1e5\nlist = 1.5, 2,3e-1\n");
	CHECK(p.size() == 4);
	CHECK(Block::getParInt("a", p, true, 0) == 1);
	CHECK_NEAR(Block::getParDouble("b", p, true, 0), 2.5e-3, 1e-18);
	CHECK(Block::getParInt("maxIter", p, true, 0) == 100000);
	CHECK(Block::getParInt("missing", p, false, 42) == 42);
	CHECK(!Block::hasPar("c", p));
	double* arr = Block::getParDoubleArray("list", p);
	CHECK_NEAR(arr[0], 1.5, 1e-15);
	CHECK_NEAR(arr[1], 2.0, 1e-15);
	CHECK_NEAR(arr[2], 0.3, 1e-15);
	delete[] arr;
	ParamMap big = parseParams("saveRate = 1e10");
	CHECK(Block::getParInt("saveRate", big, true, 0) == 2147483647);
}

TEST(boundary_file_cell_types) {
	CaseSpec s;
	s.name = "boundary_types";
	s.ni = s.nj = s.nk = 4;
	s.params = latticeParams(0, 4, 0.1);
	s.cellType = [](int i, int j, int k) { return (i + 2 * j + 3 * k) % 5; };
	auto mesh = buildMesh(s);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 4; k++) {
				int ind;
				Block* b = mesh->cellAt(i, j, k, ind);
				int t = (i + 2 * j + 3 * k) % 5;
				CHECK(b->cellType[ind] == t);
				CHECK(b->boundary[ind] == (t != 0));
			}

	// invalid values are reported
	s.name = "boundary_invalid";
	s.cellType = [](int i, int j, int k) { return 9; };
	CHECK_THROWS(buildMesh(s));
}

TEST(neighbour_tables_match_local2global) {
	CaseSpec s;
	s.name = "neighbours";
	s.nBi = 3;
	s.nBj = 2;
	s.nBk = 2;
	s.ni = 5;
	s.nj = 4;
	s.nk = 3;
	s.params = latticeParams(0, 5, 0.1);
	s.params["xPeriodicity"] = "1";
	s.blockExists = [](int bi, int bj, int bk) { return !(bi == 1 && bj == 1 && bk == 1); };
	auto mesh = buildMesh(s);
	for (int b = 0; b < mesh->nBlock; b++) {
		Block* bl = mesh->blocks[mesh->blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			for (int dir = 0; dir < bl->nf; dir++) {
				int eb, el;
				bl->local2globalShift(bl->indBlock, ind, eb, el, bl->ex[dir], bl->ey[dir], bl->ez[dir]);
				if (eb >= 0 && !mesh->blocks[eb]) {
					eb = -1;
				}
				CHECK(bl->nbBlock(ind, dir) == eb);
				if (eb >= 0) {
					CHECK(bl->nbLocal(ind, dir) == el);
				}
			}
		}
	}
}

TEST(missing_block_file_is_reported) {
	filesystem::create_directories("tests_tmp/nomesh/sim");
	ParamMap p = parseParams("ni=4\nnj=4\nnk=4\nproblemID=0\nmodel=319\nuMax=0.1\nH=4\ndensity=1\nvKin=0.1");
	CHECK_THROWS(Mesh("tests_tmp/nomesh/", "tests_tmp/nomesh/sim/", p));
}
