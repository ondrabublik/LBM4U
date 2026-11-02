#include "TestCase.h"
#include "../SimulationIO.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;

namespace {
	// silences std::cout while alive
	struct Silence {
		streambuf* old;
		ostringstream sink;
		Silence() { old = cout.rdbuf(sink.rdbuf()); }
		~Silence() { cout.rdbuf(old); }
	};
}

map<string, string> latticeParams(int problemID, int ni, double vKin) {
	map<string, string> p;
	p["problemID"] = to_string(problemID);
	p["model"] = "319";
	p["collisionModel"] = "2";
	p["threads"] = "4";
	p["uMax"] = "0.1";
	p["H"] = to_string(ni);
	p["density"] = "1";
	ostringstream v;
	v.precision(17);
	v << vKin;
	p["vKin"] = v.str();
	p["movingLidTop"] = "0";
	p["writeFields"] = "0";
	return p;
}

string caseRoot(const CaseSpec& spec) {
	return "tests_tmp/" + spec.name + "/";
}

unique_ptr<Mesh> buildMesh(const CaseSpec& spec) {
	string root = caseRoot(spec);
	filesystem::remove_all(root);
	string geom = root + "geom/";
	string sim = geom + "sim/";
	filesystem::create_directories(geom + "mesh/boundary");
	filesystem::create_directories(sim + "init");
	filesystem::create_directories(sim + "result");
	filesystem::create_directories(sim + "out");

	map<string, string> params = spec.params;
	params["ni"] = to_string(spec.ni);
	params["nj"] = to_string(spec.nj);
	params["nk"] = to_string(spec.nk);
	{
		ofstream pf(sim + "parameters.txt");
		for (auto& kv : params) {
			pf << kv.first << " = " << kv.second << "\n";
		}
	}
	ofstream bf(geom + "mesh/blocks.txt");
	if (!spec.octree.empty()) {
		int nBi = 0, nBj = 0, nBk = 0;
		for (auto& b : spec.octree) {
			nBi = max(nBi, (b[0] >> b[3]) + 1);
			nBj = max(nBj, (b[1] >> b[3]) + 1);
			nBk = max(nBk, (b[2] >> b[3]) + 1);
		}
		for (auto& b : spec.octree) {
			int l = b[3];
			bf << b[0] << " " << b[1] << " " << b[2] << " " << l << "\n";
			string tag = l == 0 ? to_string(nBi * nBj * b[2] + nBi * b[1] + b[0]) : "L" + to_string(l) + "_" + to_string(b[0]) + "_" + to_string(b[1]) + "_" + to_string(b[2]);
			int g0 = b[0] * spec.ni, g1 = b[1] * spec.nj, g2 = b[2] * spec.nk;
			if (spec.cellTypeL) {
				ofstream f(geom + "mesh/boundary/blockBoundary" + tag + ".txt");
				for (int i = 0; i < spec.ni; i++)
					for (int j = 0; j < spec.nj; j++)
						for (int k = 0; k < spec.nk; k++)
							f << spec.cellTypeL(l, g0 + i, g1 + j, g2 + k) << "\n";
			}
			if (spec.flagL || spec.initL) {
				ofstream f(sim + "init/blockInit" + tag + ".txt");
				f.precision(17);
				for (int i = 0; i < spec.ni; i++)
					for (int j = 0; j < spec.nj; j++)
						for (int k = 0; k < spec.nk; k++) {
							array<double, 4> st = { 1, 0, 0, 0 };
							if (spec.initL) {
								st = spec.initL(l, g0 + i, g1 + j, g2 + k);
							}
							f << st[0] << " " << st[1] << " " << st[2] << " " << st[3];
							if (spec.flagL) {
								f << " " << spec.flagL(l, g0 + i, g1 + j, g2 + k);
							}
							f << "\n";
						}
			}
		}
	}
	for (int bk = 0; bk < spec.nBk && spec.octree.empty(); bk++) {
		for (int bj = 0; bj < spec.nBj; bj++) {
			for (int bi = 0; bi < spec.nBi; bi++) {
				if (spec.blockExists && !spec.blockExists(bi, bj, bk)) {
					continue;
				}
				bf << bi << " " << bj << " " << bk << "\n";
				int id = spec.nBi * spec.nBj * bk + spec.nBi * bj + bi;
				// files are written in the order of Mesh::indAll (i, j, k with k fastest)
				if (spec.cellType) {
					ofstream f(geom + "mesh/boundary/blockBoundary" + to_string(id) + ".txt");
					for (int i = 0; i < spec.ni; i++)
						for (int j = 0; j < spec.nj; j++)
							for (int k = 0; k < spec.nk; k++)
								f << spec.cellType(bi * spec.ni + i, bj * spec.nj + j, bk * spec.nk + k) << "\n";
				}
				if (spec.flag) {
					ofstream f(sim + "init/blockInit" + to_string(id) + ".txt");
					for (int i = 0; i < spec.ni; i++)
						for (int j = 0; j < spec.nj; j++)
							for (int k = 0; k < spec.nk; k++)
								f << "1 0 0 0 " << spec.flag(bi * spec.ni + i, bj * spec.nj + j, bk * spec.nk + k) << "\n";
				}
			}
		}
	}
	bf.close();

	Silence s;
	ParamMap par = loadParamFile(sim + "parameters.txt");
	unique_ptr<Mesh> mesh(new Mesh(geom, sim, par));
	mesh->initBlocks();
	return mesh;
}

void runSteps(Mesh& mesh, int n, int startIter) {
	for (int it = startIter; it < startIter + n; it++) {
		mesh.timeStep(it);
	}
}

array<double, 3> velocityAt(Mesh& mesh, int gi, int gj, int gk) {
	int ind;
	Block* b = mesh.cellAt(gi, gj, gk, ind);
	array<double, 3> u = { 0, 0, 0 };
	if (b && !b->boundary[ind]) {
		b->outputVelocity(ind, u.data());
		for (double& x : u) {
			x *= b->Cu;
		}
	}
	return u;
}

vector<array<int, 4> > octreeBlocks(int nBi, int nBj, int nBk, int maxLevel, function<bool(int, int, int, int)> refine) {
	vector<array<int, 4> > out;
	function<void(int, int, int, int)> add = [&](int l, int i, int j, int k) {
		if (l < maxLevel && refine(l, i, j, k)) {
			for (int o = 0; o < 8; o++) {
				add(l + 1, 2 * i + (o & 1), 2 * j + ((o >> 1) & 1), 2 * k + (o >> 2));
			}
		}
		else {
			out.push_back({ i, j, k, l });
		}
	};
	for (int k = 0; k < nBk; k++)
		for (int j = 0; j < nBj; j++)
			for (int i = 0; i < nBi; i++)
				add(0, i, j, k);
	return out;
}

Block* blockAtPoint(Mesh& mesh, double x, double y, double z, int& ind) {
	for (int l = mesh.maxLevel; l >= 0; l--) {
		double s = (double)(1 << l);
		Block* b = mesh.cellAtLevel(l, (int)floor(x * s), (int)floor(y * s), (int)floor(z * s), ind);
		if (b && b->ghost == GHOST_NONE) {
			return b;
		}
	}
	return nullptr;
}

double latticeMass(Mesh& mesh) {
	double m = 0;
	for (int b = 0; b < mesh.nBlock; b++) {
		Block* bl = mesh.blocks[mesh.blockId[b]];
		double v = 1.0 / (1 << (3 * bl->level));
		for (int ind = 0; ind < bl->nAll; ind++) {
			if (!bl->boundary[ind]) {
				m += bl->Rho[ind] * v;
			}
		}
	}
	return m;
}
