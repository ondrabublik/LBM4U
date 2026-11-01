/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

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

	double* concentration, * gc, * g;
	double* UP, * VP, * WP;
	double* gUx, * gUy, * gUz, * gVx, * gVy, * gVz, * gWx, * gWy, * gWz;
	double* UF, * VF, * WF;
	double* Gx, * Gy, * Gz;
	double* GxSmooth, * GySmooth, * GzSmooth;
	double* dpX, * dpY, * dpZ;
	double* dcX, * dcY, * dcZ;
	double* gVOFpX, * gVOFpY, * gVOFpZ;
	double dParticle, rhoParticle, crossectionParticle, volumeParticle, massParticle;
	double taus, diffusLB;
	double* deltaValue = nullptr;
	int interpLevel;
	double concentrationLimit, interPartFCoef, concentrationSum, concentrationMax;
	int couplingWayNumber;
	double* tauDir;
	double tauOverPack;
	int *positiveDirs;
	double theta;
	double* user_def;

	BlockFreeSurfaceEulerParticlesLBM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par);

	~BlockFreeSurfaceEulerParticlesLBM();

	void collision();
	void cellForce(int ind, double* acc, double& epsF);
	void equilibriumScalar(double* geq, double c, double Up, double Vp, double Wp);
	void computeInitDistributionScalar();
	void readInitConcentration(int id, string simulPath);
	void calculConcentration();
	void collisionScalar();
	void propagationScalar();
	void computeGradP();
	void computeGradU();
	void computeGradC();
	void computeGradVOFp();
	void computeNewParticleVelocityField(double* distance, double* gradDistance);
	double cd(double Uslip);
	void interpolateFluidVelocity();
	void smoothForce();
	void setTau();
	double computeTau(double cFor, double cMid, double cBack);
	double delta(double x);
	double limitFraction(double x);
};

#endif /* BLOCKFREESURFACEEULERPARTICLESLBM_H */

