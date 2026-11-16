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
#include <vector>
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
	// Lagrangian particles (problemID 2): x, y, z, diameter, density in physical units -> init/particles.txt
	std::vector<std::array<double, 5> > particles;
	// any particle lines (e.g. ellipsoids of problemID 22: x y z a1 a2 a3 density qw qx qy qz [f1 f2 f3]),
	// used instead of particles when not empty
	std::vector<std::vector<double> > particleRows;
	// Eulerian particles (problemID 3, 4): number of particles in a global cell -> init/concentration<id>.txt
	std::function<double(int, int, int)> concentration;

	// refined mesh: blocks (bi, bj, bk, level) replace the nBi x nBj x nBk grid, the functions below
	// take the level and the global cell indexes of the level and replace cellType and flag
	std::vector<std::array<int, 4> > octree;
	std::function<int(int, int, int, int)> cellTypeL;
	std::function<int(int, int, int, int)> flagL;
	// initial density and velocity (lattice units of the level 0) of the single phase solver
	std::function<std::array<double, 4>(int, int, int, int)> initL;
	// number of particles in a cell of the level (Eulerian particles on a refined mesh)
	std::function<double(int, int, int, int)> concentrationL;
};

// leaves of an octree over nBi x nBj x nBk blocks, refine(level, bi, bj, bk) splits a block
std::vector<std::array<int, 4> > octreeBlocks(int nBi, int nBj, int nBk, int maxLevel, std::function<bool(int, int, int, int)> refine);

// computed block holding the point given in the cells of the level 0 (the finest one), nullptr outside
Block* blockAtPoint(Mesh& mesh, double x, double y, double z, int& ind);

// sum of rho * cell volume over the fluid cells of the computed blocks (volume of a level 0 cell is 1)
double latticeMass(Mesh& mesh);

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
