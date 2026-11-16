/*
 * File:   LagrangeParticles.h
 * Author: obublik
 *
 * Created on 9. ledna 2019, 14:27
 *
 * Lagrangian point particles coupled with the free surface LBM (problemID 2): drag (Schiller-Naumann),
 * gravity and buoyancy, repulsive forces between the particles, reflection from the walls, the reaction
 * of the drag is spread to the liquid cells around the particle.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "math.h"
#include <string>
#include <unordered_map>
#include <functional>

#include "Block.h"

#ifndef LAGRANGEPARTICLES_H
#define LAGRANGEPARTICLES_H

using namespace std;

class BlockFreeSurfaceLagrangeParticles;

class LagrangeParticles {
public:

	double Ch, Cu, Ct, Ca, Crho, Cvis, Cf;
	double vLB;
	double G[3], G0[3];
	double restitution; // of the normal velocity at the walls (1 - elastic reflection)
	Block** blocks;
	const BlockGrid* grid = nullptr; // block grid of the mesh (periodicity), set by the Mesh
	int nBlock, nBlockMax;
	int nBi, nBj, nBk;
	int ni, nj, nk;
	int nParticles = 0;
	double dMax = 0; // largest diameter (range of the repulsive force)
	// all quantities are in lattice units of the level 0, the fluid density is one
	vector<double> partD, partRho, partMass, partVolume, partCrossection;
	vector<double> partX, partY, partZ;
	vector<double> veloX, veloY, veloZ;
	vector<double> veloFX, veloFY, veloFZ;   // fluid velocity at the particle
	vector<double> pX, pY, pZ;               // pressure gradient at the particle (output only)
	vector<double> Gx, Gy, Gz;               // impulse given to the fluid during the LBM step
	vector<double> FinterX, FinterY, FinterZ;
	vector<double> volumeSolid, areaSolid;
	vector<int> nInter;
	vector<double> partWallNormalX, partWallNormalY, partWallNormalZ;
	vector<char> inGas;                      // the particle is not in the liquid (free flight, no coupling)
	int interpLevel;

	LagrangeParticles(Block** blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, const ParamMap& par);

	virtual ~LagrangeParticles() = default;
	virtual void readInitParticlesDistribution(string simulPath, double Ch, double Cu);
	void saveParticlesPositionsToVTK(string simulPath, int iter);
	virtual void interpolateFluidVelocity();
	void applyForce();
	void boxParticles();
	virtual void newParticlePosition(double dt);
	void particleParticleInteraction();
	double gammaPP(double L, double avgD, double uRel);
	double gammaWall(double wallDistance, double D);
	double delta(double x);
	double cd(double d, double uSlip);
	void getParticlesInfo();
	void particleBarrierForce(int iter);
	// distance of the point from the wall (-1 inside of a wall or outside of the domain), sets the wall normal of the particle
	double wallDistance(int p, double x, double y, double z);

	// cell holding a point: block slot, cell indexes and the coordinates of the point in the block
	struct Cell {
		int bl, i, j, k;
		int level = 0;   // refinement level of the block, the coordinates below are in its cells
		double x, y, z;
	};
	// false if no computed block holds the point (the finest one is taken on a refined mesh)
	bool locate(double x, double y, double z, Cell& c) const;

protected:
	void resize(int n);
	// runs fn(block, cell, weight) for the liquid cells in the support of the delta function of the point
	template <class F> void liquidCells(const Cell& c, F fn);
	// the same for the derived models (the template is defined in LagrangeParticles.cpp)
	void forLiquidCells(const Cell& c, const function<void(BlockFreeSurfaceLagrangeParticles*, int, double)>& fn);
	void wallContact(int p, double x0, double y0, double z0, double L0, const double* n0);
	// distance of the particle surface from its centre along the current wall normal (partWallNormal*)
	virtual double contactRadius(int p) const { return partD[p] / 2; }
	// additional particle data of the derived models in output_particles*.vtk (CELL_DATA section)
	virtual void writeExtraVTK(ofstream& f) {}
};

#endif /* LAGRANGEPARTICLES_H */
