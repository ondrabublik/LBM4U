/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

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

    double *concentration, *fluxConcentration;
	double *UP, *VP, *WP;
	double *PP;
	double* dPPx, * dPPy, * dPPz;
	double *UF, *VF, *WF;
	double *Gx, *Gy, *Gz;
	double *GxSmooth, *GySmooth, *GzSmooth;
	double *dpX, *dpY, *dpZ;
	double *fcP;
	double dParticle, rhoParticle, crossectionParticle, volumeParticle, massParticle;
	double *deltaValue;
	int interpLevel;
	double Cmax;
	int *exLW, *eyLW, *ezLW;
	int *oppLW;

    BlockFreeSurfaceEulerParticlesFVM(bool isFirst, Block **blocks, int nBlock, int nBlockMax, int *blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int *indAll, int *indIn, int *indBound, bool isMRTModel, int viscosityModel, int nf, int *ex, int *ey, int *ez, double *w, int* opp, int ni, int nj, int nk, unordered_map<string, string> par);

    ~BlockFreeSurfaceEulerParticlesFVM();

	void readInitConcentration(int id, string simulPath);
	void collision();
	void computeGradP();
	void computeNewParticleVelocityField();
	double cd(double Uslip);
	void interpolateFluidVelocity();
	void smoothForce();
	double delta(double x);
	void upwind();
	void computeParticlePressure();
	void computeParticlePressureGradient();
};

#endif /* BLOCKFREESURFACEEULERPARTICLES_H */

