/*
 * File:   BlockFreeSurfaceEulerParticlesFVM.h
 * Author: obublik
 *
 * Created on 10. ledna 2019, 12:31
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <unordered_map>
#include "math.h"
#include "BlockFreeSurface.h"

#ifndef BLOCKFREESURFACEEULERPARTICLES_H
#define BLOCKFREESURFACEEULERPARTICLES_H

class BlockFreeSurfaceEulerParticlesFVM : public BlockFreeSurface {
public:

	double* concentration, * fluxConcentration; // number of particles in the cell
	double* released;                           // particles of a cell without liquid given to the liquid cells around
	// refined mesh: computed blocks - particles coming from the finer level (added in applyFlux),
	// ghost blocks of a coarser level - particles the computed cells sent into the ghost cell (Mesh moves them)
	double* exchangeCorr;
	bool anyReleased = false;
	double* UP, * VP, * WP;                     // particle velocity
	double* PP;
	double* dPPx, * dPPy, * dPPz;
	double* UF, * VF, * WF;                     // fluid velocity seen by the particles
	double* Gx, * Gy, * Gz;                     // force of the particles of the cell acting on the fluid
	double* GxSmooth, * GySmooth, * GzSmooth;   // force spread to the liquid cells around
	double* dpX, * dpY, * dpZ;
	double dParticle, rhoParticle, crossectionParticle, volumeParticle, massParticle;
	double Cmax;

	BlockFreeSurfaceEulerParticlesFVM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par);

	~BlockFreeSurfaceEulerParticlesFVM();

	void readInitConcentration(int id, string simulPath);
	void collision();
	void cellForce(int ind, double* acc, double& epsF);
	void computeGradP();
	void computeNewParticleVelocityField();
	double cd(double Uslip);
	void interpolateFluidVelocity();
	void normalizeForce();
	void smoothForce();
	void upwind();
	// flux of concentration through face p (of the 6 directions) leaving cell ind
	double faceFlux(int ind, int p);
	double upwindFlux(int ind, int p, BlockFreeSurfaceEulerParticlesFVM* bs, int loc, double Un);
	void computeFlux();
	void applyFlux();
	// the particles follow the liquid: cells that lost the liquid release their particles (returns true if any),
	// the liquid cells around collect them in a separate pass
	bool releaseDryCells();
	void collectReleased();
	void computeParticlePressure();
	// ghost cells: the concentration is the number of particles in the cell, so it is divided (summed) by 8
	void prolongCell(int ind, const Block* src, int sInd);
	void restrictCell(int ind, Block* const* src, const int* sInd);
	void computeParticlePressureGradient();

private:
	void liquidGradient(const double* (*field)(Block*), double scale, double* outX, double* outY, double* outZ);
};

#endif /* BLOCKFREESURFACEEULERPARTICLES_H */
