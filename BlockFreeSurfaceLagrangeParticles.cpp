
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 /*
  * File:   Block.cpp
  * Author: obublik
  *
  * Created on 2. února 2018, 12:27
  */

#include "BlockFreeSurfaceLagrangeParticles.h"

using namespace std;

BlockFreeSurfaceLagrangeParticles::BlockFreeSurfaceLagrangeParticles(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurface(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {

	particleForceX = new double[ntot]();
	particleForceY = new double[ntot]();
	particleForceZ = new double[ntot]();

	dpX = new double[ntot]();
	dpY = new double[ntot]();
	dpZ = new double[ntot]();

	nParticleBox = getParInt("maxParticlesInCell", par, false, 500);
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

BlockFreeSurfaceLagrangeParticles::~BlockFreeSurfaceLagrangeParticles() {
	delete[] particleForceX;
	delete[] particleForceY;
	delete[] particleForceZ;
	delete[] particleBox;
	delete[] dpX;
	delete[] dpY;
	delete[] dpZ;
	delete[] concentration;
	delete[] particleFraction;
	delete[] boxVolume;
	delete[] boxArea;
}

void BlockFreeSurfaceLagrangeParticles::cellForce(int ind, double* acc, double& epsF) {
	// particle force is a force density, the exact difference method needs an acceleration
	double rho = Rho[ind] > 0 ? Rho[ind] : 1.0;
	epsF = (couplingWayNumber == 4) ? 1 - particleFraction[ind] : 1.0;
	acc[0] = G[0] * epsF + particleForceX[ind] * vof[ind] / rho;
	acc[1] = G[1] * epsF + particleForceY[ind] * vof[ind] / rho;
	acc[2] = G[2] * epsF + particleForceZ[ind] * vof[ind] / rho;
}

double BlockFreeSurfaceLagrangeParticles::localTau(int ind) {
	if (cLES == 0) {
		return BlockFreeSurface::localTau(ind);
	}
	double tauLoc = tau + 3 * cLES * concentration[ind] * computeShearRate(ind);
	return min(max(tauLoc, tauMin), tauMax);
}

void BlockFreeSurfaceLagrangeParticles::collision() {
	BlockFreeSurface::collision();
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		particleForceX[ind] = 0;
		particleForceY[ind] = 0;
		particleForceZ[ind] = 0;
		particleBox[ind * (nParticleBox + 1)] = 0;
		boxVolume[ind] = 0;
		boxArea[ind] = 0;
	}
}


void BlockFreeSurfaceLagrangeParticles::computeGradP() {
	int indBlockShift = 0;
	int indShift = 0;
	double* auxRho = new double[nf];

	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indShift = ind + (ni * nj) * ez[dir] + ni * ey[dir] + ex[dir];
				bool dirExist = false;
				if (!boundary[indShift] && (flagMap[indShift] == 1 || flagMap[indShift] == 3)) {
					auxRho[dir] = Rho[indShift];
					dirExist = true;
				}
				if (!dirExist) {
					auxRho[dir] = Rho[ind];
				}
			}

			double px = 0;
			double py = 0;
			double pz = 0;
			for (int dir = 1; dir < nf; dir++) {
				px += w[dir] * ex[dir] * auxRho[dir];
				py += w[dir] * ey[dir] * auxRho[dir];
				pz += w[dir] * ez[dir] * auxRho[dir];
			}
			dpX[ind] = 1.5 * px / 3;
			dpY[ind] = 1.5 * py / 3;
			dpZ[ind] = 1.5 * pz / 3;
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indBlockShift = nbBlock(ind, dir);
				int indShift = nbLocal(ind, dir);
				bool dirExist = false;
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
					if (!(blocks[indBlockShift]->boundary[indShift]) && (fm == 1 || fm == 3)) {
						auxRho[dir] = ((BlockFreeSurface*)blocks[indBlockShift])->Rho[indShift];
						dirExist = true;
					}
				}
				if (!dirExist) {
					auxRho[dir] = Rho[ind];
				}
			}

			double px = 0;
			double py = 0;
			double pz = 0;
			for (int dir = 1; dir < nf; dir++) {
				px += w[dir] * ex[dir] * auxRho[dir];
				py += w[dir] * ey[dir] * auxRho[dir];
				pz += w[dir] * ez[dir] * auxRho[dir];
			}
			dpX[ind] = 1.5 * px / 3;
			dpY[ind] = 1.5 * py / 3;
			dpZ[ind] = 1.5 * pz / 3;
		}
	}
}

void BlockFreeSurfaceLagrangeParticles::computeConcentrationAndFluidFraction() {
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
					concentration[ind] += ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->particleBox[indShift * (nParticleBox + 1)];
					particleFraction[ind] += ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->boxVolume[indShift];
					sum++;
				}
			}
		}
		concentration[ind] /= sum;
		particleFraction[ind] /= sum;
	}
}