/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   lagrangeParticles.h
 * Author: obublik
 *
 * Created on 9. ledna 2019, 14:27
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "math.h"
#include <windows.h>
#include <string>
#include <unordered_map>

#include "Block.h"

#ifndef LAGRANGEPARTICLES_H
#define LAGRANGEPARTICLES_H

using namespace std;

class LagrangeParticles {
public:

    double Ch, Cu, Ct, Ca, Crho, Cvis, Cf;
	double vLB;
    double *G, *G0;
    Block **blocks;
    int nBlock, nBlockMax;
    int *blockMap;
    int nBi, nBj, nBk;
    int ni, nj, nk;
    int nParticles;
	double* partD, * partRho, * partMass, * partVolume, *partCrossection;
	double* partX, * partY, * partZ;
    double *veloX, *veloY, *veloZ;
	double* veloFX, * veloFY, * veloFZ;
	double* veloX0, * veloY0, * veloZ0;
    double *pX, *pY, *pZ;
	double* Gx, * Gy, * Gz;
	double* FinterX, * FinterY, * FinterZ;
	double* FwallX, * FwallY, * FwallZ;
	double* volumeSolid, * areaSolid;
	int* nInter;
	double* partWallNormalX, *partWallNormalY, *partWallNormalZ, *partWallDistance;
	bool* inGas;
	int interpLevel;
    
    LagrangeParticles(Block **blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, unordered_map<string, string> par);
    
    ~LagrangeParticles();

    void readInitParticlesDistribution(string simulPath, double Ch, double Cu);
    void saveParticlesPositionsToVTK(string simulPath, int iter);
    void interpolateFluidVelocity();
    void applyForce();
	void boxParticles();
    void newParticlePosition(double dt);
	void particleParticleInteraction();
	void particleWallInteraction();
	double gammaPP(double L, double avgD, double uRel);
	double gammaWall(double wallDistance, double D);
	double delta(double x);
	double cd(double d, double uSlip);
	void getParticlesInfo();
	void particleBarrierForce(int iter);
	double wallDistance(int p, double x, double y, double z);
};

#endif /* LAGRANGEPARTICLES_H */

