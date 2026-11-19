/*
 * File:   BlockFreeSurfaceEulerParticlesLBM.h
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

#ifndef BLOCKFREESURFACEEULERPARTICLESLBM_H
#define BLOCKFREESURFACEEULERPARTICLESLBM_H

class BlockFreeSurfaceEulerParticlesLBM : public BlockFreeSurface {
public:

	double* concentration, * gc, * g;           // number of particles in the cell and its populations
	double* released;                           // particles of a cell without liquid given to the liquid cells around
	// refined mesh: change of the sum of the populations of a computed cell balancing the exchange with the finer
	// level (Mesh::particle*), added in applyExchangeCorr
	double* exchangeCorr;
	bool anyReleased = false;
	double* UP, * VP, * WP;                     // particle velocity
	double* gUx, * gUy, * gUz, * gVx, * gVy, * gVz, * gWx, * gWy, * gWz;
	double* UF, * VF, * WF;                     // fluid velocity seen by the particles
	double* Gx, * Gy, * Gz;                     // force of the particles of the cell acting on the fluid
	double* GxSmooth, * GySmooth, * GzSmooth;   // force spread to the liquid cells around
	double* dpX, * dpY, * dpZ;
	double* dcX, * dcY, * dcZ;
	double* gVOFpX, * gVOFpY, * gVOFpZ;
	double dParticle, rhoParticle, crossectionParticle, volumeParticle, massParticle;
	double taus, diffusLB;
	double tauPlus;                             // TRT relaxation time of the symmetric part
	double tauFront;                            // extra tau- at steep fronts (frontDiffusivity)
	static constexpr double TRT_TAU_PLUS = 1.0;
	double concentrationLimit, interPartFCoef, concentrationSum, concentrationMax;
	double* tauDir;                             // tau- of the pairs of directions
	double tauOverPack;
	int* positiveDirs;
	double theta;
	double* user_def;

	BlockFreeSurfaceEulerParticlesLBM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par);

	~BlockFreeSurfaceEulerParticlesLBM();

	void collision();
	void cellForce(int ind, double* acc, double& epsF);
	void equilibriumScalar(double* geq, double c, double Up, double Vp, double Wp);
	void computeInitDistributionScalar();
	virtual void readInitConcentration(int id, string simulPath);
	void calculConcentration();
	void collisionScalar();
	void propagationScalar();
	virtual void propagateScalarCell(int ind);
	// streaming of the listed cells only (ghost layers of a refined mesh)
	void propagateScalarCells(const vector<int>& cells);
	// adds exchangeCorr to the populations (refined mesh), after the propagation of the scalar
	virtual void applyExchangeCorr();
	// ghost cells: the populations hold numbers of particles, so they are divided (summed) by 8
	void prolongCell(int ind, const Block* src, int sInd);
	void restrictCell(int ind, Block* const* src, const int* sInd);
	// the particles follow the liquid: cells that lost the liquid release their particles (returns true if any),
	// the liquid cells around collect them in a separate pass
	virtual bool releaseDryCells();
	virtual void collectReleased();
	void computeGradP();
	void computeGradU();
	void computeGradC();
	void computeGradVOFp();
	// gradient of the particle velocity and the fluid velocity seen by the particles (before the particle velocity)
	virtual void prepareParticleStep();
	virtual void computeNewParticleVelocityField(double* distance, double* gradDistance);
	double cd(double Uslip);
	void interpolateFluidVelocity();
	void normalizeForce();
	void smoothForce();
	void setTau();
	double computeTau(double cFor, double cMid, double cBack);
	double limitFraction(double x);

protected:
	// next to a wall the particles do not move towards it faster than the fluid
	void wallSlipLimit(int ind, const double* dist, const double* gradDist);

private:
	void fieldGradient(const double* (*field)(Block*), double scale, bool liquidOnly, double* outX, double* outY, double* outZ);
};

#endif /* BLOCKFREESURFACEEULERPARTICLESLBM_H */
