
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 /*
  * File:   BlockFreeSurfaceLagrangeFibers.h
  * Author: obublik
  *
  * Created on 9. ledna 2019, 13:27
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

#ifndef BLOCKFREESURFACELAGRANGEFIBERS_H
#define BLOCKFREESURFACELAGRANGEFIBERS_H

class BlockFreeSurfaceLagrangeFibers : public BlockFreeSurface {
public:

	double* particleForceX, * particleForceY, * particleForceZ;
	int nParticleBox;
	int* particleBox;
	double* boxVolume, * boxArea;
	double* concentration;
	double* particleFraction;
	double cLES;
	int couplingWayNumber;

	BlockFreeSurfaceLagrangeFibers(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int* indAll, int* indIn, int* indBound, bool isMRTModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, unordered_map<string, string> par);

	~BlockFreeSurfaceLagrangeFibers();

	void collision();
	void computeConcentrationAndFluidFraction();
};

#endif /* BLOCKFREESURFACELAGRANGEFIBERS_H */

