#include "TestCase.h"
#include "../SimulationIO.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

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
	for (int bk = 0; bk < spec.nBk; bk++) {
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
