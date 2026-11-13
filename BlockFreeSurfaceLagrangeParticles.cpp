/*
 * File:   BlockFreeSurfaceLagrangeParticles.cpp
 * Author: obublik
 *
 * Created on 2. února 2018, 12:27
 *
 * Free surface block coupled with the Lagrangian particles (LagrangeParticles).
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

	nParticleBox = max(getParInt("maxParticlesInCell", par, false, 500), 1);
	particleBox = new int[(size_t)ntot * (nParticleBox + 1)]();

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
	acc[0] = fluidGravity * G[0] * epsF + particleForceX[ind] * vof[ind] / rho;
	acc[1] = fluidGravity * G[1] * epsF + particleForceY[ind] * vof[ind] / rho;
	acc[2] = fluidGravity * G[2] * epsF + particleForceZ[ind] * vof[ind] / rho;
}

double BlockFreeSurfaceLagrangeParticles::localTau(int ind) {
	if (cLES == 0) {
		return BlockFreeSurface::localTau(ind);
	}
	double tauLoc = tau + 3 * cLES * concentration[ind] * computeShearRate(ind);
	return min(max(tauLoc, tauMin), tauMax);
}

// the forces of the particles act until LagrangeParticles::applyForce replaces them (all steps of a finer level)
void BlockFreeSurfaceLagrangeParticles::collision() {
	BlockFreeSurface::collision();
}

// dp = grad(p) / 2 with p = rho / 3 in the liquid cells (interpolated to the particles for the output)
void BlockFreeSurfaceLagrangeParticles::computeGradP() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] != 1 && flagMap[ind] != 3) {
			dpX[ind] = 0;
			dpY[ind] = 0;
			dpZ[ind] = 0;
			continue;
		}
		double px = 0, py = 0, pz = 0;
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			double rho = Rho[ind];
			if (bl > -1 && bl < nBlockMax && blocks[bl] && !blocks[bl]->boundary[loc]) {
				int fm = ((BlockFreeSurface*)blocks[bl])->flagMap[loc];
				if (fm == 1 || fm == 3) {
					rho = blocks[bl]->Rho[loc];
				}
			}
			px += w[dir] * ex[dir] * rho;
			py += w[dir] * ey[dir] * rho;
			pz += w[dir] * ez[dir] * rho;
		}
		dpX[ind] = 0.5 * px;
		dpY[ind] = 0.5 * py;
		dpZ[ind] = 0.5 * pz;
	}
}

// number of particles and their volume fraction averaged over the cell and its neighbours
void BlockFreeSurfaceLagrangeParticles::computeConcentrationAndFluidFraction() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		concentration[ind] = particleBox[ind * (nParticleBox + 1)];
		particleFraction[ind] = boxVolume[ind];
		int sum = 1;
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			if (bl > -1 && bl < nBlockMax && blocks[bl] && !blocks[bl]->boundary[loc]) {
				BlockFreeSurfaceLagrangeParticles* b = (BlockFreeSurfaceLagrangeParticles*)blocks[bl];
				concentration[ind] += b->particleBox[loc * (nParticleBox + 1)];
				particleFraction[ind] += b->boxVolume[loc];
				sum++;
			}
		}
		concentration[ind] /= sum;
		particleFraction[ind] /= sum;
	}
}
