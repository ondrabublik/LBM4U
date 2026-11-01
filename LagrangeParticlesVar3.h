
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

#ifndef LAGRANGEPARTICLESVAR3_H
#define LAGRANGEPARTICLESVAR3_H

using namespace std;

class LagrangeParticlesVar3 {
public:

	double Ch, Cu, Ct, Cf, CRho, Cm, Cvis;
	double vLB;
	double* G;
	Block** blocks;
	int nBlock, nBlockMax;
	int* blockMap;
	int nBi, nBj, nBk;
	int ni, nj, nk, nf;
	int nParticles;
	double* partCrossection, *partMass;
	double* partD, * partRho;
	double* partX, * partY, * partZ;
	double* veloX, * veloY, * veloZ;
	double* veloX0, * veloY0, * veloZ0;
	double* Gx, * Gy, * Gz; // force density
	double* FinterX, * FinterY, * FinterZ;
	int nVeloInterp, nForceInterp;
	double wForceNorm;

	LagrangeParticlesVar3(Block** blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, int nf, unordered_map<string, string> par);

	~LagrangeParticlesVar3();

	void readInitParticlesDistribution(string simulPath, double Ch, double Cu);
	void saveParticlesPositionsToVTK(string simulPath, int iter);
	void getParticleForce();
	void particleParticleInteraction();
	void applyForce();
	void newParticlePosition();
	double delta(double x);
	double gammaPP(double x, double avgD);
	double gammaWall(double x);
	double cd(double Re);
};

#endif /* LAGRANGEPARTICLESVAR3_H */
