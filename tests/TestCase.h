/*
 * Helpers that write a complete simulation case (mesh, boundary, parameters, initial condition)
 * into a temporary directory and build the Mesh from it.
 */

#ifndef TESTCASE_H
#define TESTCASE_H

#include <string>
#include <map>
#include <functional>
#include <memory>
#include <array>
#include "../Mesh.h"

// relative tolerance of quantities affected by the precision of the populations
const double POP_TOL = sizeof(pop_t) == 4 ? 1e-6 : 1e-10;

struct CaseSpec {
	std::string name;
	int nBi = 1, nBj = 1, nBk = 1;
	int ni = 8, nj = 8, nk = 8;
	std::map<std::string, std::string> params;
	// cell type of a global cell (CT_FLUID if not set)
	std::function<int(int, int, int)> cellType;
	// free surface flag of a global cell (1 fluid, 2 gas), no init files if not set
	std::function<int(int, int, int)> flag;
	// optional block existence filter
	std::function<bool(int, int, int)> blockExists;
};

// base parameter set: lattice units equal physical units (Ch = Ct = Cu = 1)
std::map<std::string, std::string> latticeParams(int problemID, int ni, double vKin);

// writes the case and returns the Mesh (initialised); output of the solver is suppressed
std::unique_ptr<Mesh> buildMesh(const CaseSpec& spec);

// runs n time steps without console output
void runSteps(Mesh& mesh, int n, int startIter = 0);

// physical velocity of a global cell
std::array<double, 3> velocityAt(Mesh& mesh, int gi, int gj, int gk);

std::string caseRoot(const CaseSpec& spec);

#endif
