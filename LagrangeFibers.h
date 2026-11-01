/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 /*
  * File:   lagrangeFibers.h
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
#include <string>
#include <unordered_map>

#include "Block.h"

#ifndef LAGRANGEFIBERS_H
#define LAGRANGEFIBERS_H

using namespace std;

class LagrangeFibers {
public:

	double Ch, Cu, Ct, Ca, Crho, Cvis, Cf;
	double vLB;
	double* G;
	Block** blocks;
	int nBlock, nBlockMax;
	int* blockMap;
	int nBi, nBj, nBk;
	int ni, nj, nk;
	int nFibers;
	int nChain;
	double rChain;
	double* partD, * partRho, * partMass, * partVolume, * partCrossection;
	double* fibreCoords, *fibreVelo, *fibreFluidVelo;
	double* fibreFluidForce;
	double* Fwall;
	double* volumeSolid, * areaSolid;
	double* partWallNormal;
	bool* inGas;
	int interpLevel;

	LagrangeFibers(Block** blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, unordered_map<string, string> par);

	~LagrangeFibers();

	void readInitFibersDistribution(string simulPath, double Ch, double Cu);
	void saveFibersPositionsToVTK(string simulPath, int iter);
	void interpolateFluidVelocity();
	void applyForce();
	void newFiberPosition(double dt);
	void fiberWallInteraction();
	double gammaWall(double wallDistance, double D);
	double delta(double x);
	double cd(double d, double uSlip);
	double wallDistance(int p, double x, double y, double z);
};

#endif /* LAGRANGEFIBERS_H */


