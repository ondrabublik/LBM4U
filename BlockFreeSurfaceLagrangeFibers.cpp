
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 /*
  * File:   Block.cpp
  * Author: obublik
  *
  * Created on 2. �nora 2018, 12:27
  */

#include "BlockFreeSurfaceLagrangeFibers.h"

using namespace std;

BlockFreeSurfaceLagrangeFibers::BlockFreeSurfaceLagrangeFibers(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurface(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, coarseLevel, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {

	particleForceX = new double[ntot]();
	particleForceY = new double[ntot]();
	particleForceZ = new double[ntot]();

	nParticleBox = getParInt("maxParticlesInCell", par, false, 200);
	particleBox = new int[ntot * (nParticleBox + 1)]();

	concentration = new double[ntot]();
	particleFraction = new double[ntot]();
	boxVolume = new double[ntot]();
	boxArea = new double[ntot]();

	cLES = getParDouble("cLES", par, false, 0);
	couplingWayNumber = getParInt("couplingWayNumber", par, false, 2);

	if (isFirst) {
		cout << "number of coupling ways = " << couplingWayNumber << endl;
	}
}

BlockFreeSurfaceLagrangeFibers::~BlockFreeSurfaceLagrangeFibers() {
	delete[] particleForceX;
	delete[] particleForceY;
	delete[] particleForceZ;
	delete[] particleBox;
	delete[] concentration;
	delete[] particleFraction;
	delete[] boxVolume;
	delete[] boxArea;
}

void BlockFreeSurfaceLagrangeFibers::cellForce(int ind, double* acc, double& epsF) {
	// fiber force is a force density, the exact difference method needs an acceleration
	double rho = Rho[ind] > 0 ? Rho[ind] : 1.0;
	epsF = (couplingWayNumber == 4) ? 1 - particleFraction[ind] : 1.0;
	acc[0] = G[0] * epsF + particleForceX[ind] * vof[ind] / rho;
	acc[1] = G[1] * epsF + particleForceY[ind] * vof[ind] / rho;
	acc[2] = G[2] * epsF + particleForceZ[ind] * vof[ind] / rho;
}

double BlockFreeSurfaceLagrangeFibers::localTau(int ind) {
	if (cLES == 0) {
		return BlockFreeSurface::localTau(ind);
	}
	double tauLoc = tau + 3 * cLES * computeShearRate(ind);
	return min(max(tauLoc, tauMin), tauMax);
}

void BlockFreeSurfaceLagrangeFibers::collision() {
	BlockFreeSurface::collision();
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		particleForceX[ind] = 0;
		particleForceY[ind] = 0;
		particleForceZ[ind] = 0;
	}
}

void BlockFreeSurfaceLagrangeFibers::computeConcentrationAndFluidFraction() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		int indBlockShift = 0;
		int indShift = 0;
		concentration[ind] = particleBox[ind * (nParticleBox + 1)];
		particleFraction[ind] = boxVolume[ind];
		int sum = 1;
		for (int dir = 1; dir < nf; dir++) {
			indBlockShift = nbBlock(ind, dir);
			indShift = nbLocal(ind, dir);
			if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
				if (!(blocks[indBlockShift]->boundary[indShift])) {
					concentration[ind] += ((BlockFreeSurfaceLagrangeFibers*)blocks[indBlockShift])->particleBox[indShift * (nParticleBox + 1)];
					particleFraction[ind] += ((BlockFreeSurfaceLagrangeFibers*)blocks[indBlockShift])->boxVolume[indShift];
					sum++;
				}
			}
		}
		concentration[ind] /= sum;
		particleFraction[ind] /= sum;
	}
}