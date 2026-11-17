/*
 * File:   BlockFreeSurfaceLagrangeNoSphereParticles.cpp
 *
 * Free surface block for the non-spherical Lagrangian particles (problemID 22).
 */

#include "BlockFreeSurfaceLagrangeNoSphereParticles.h"

using namespace std;

BlockFreeSurfaceLagrangeNoSphereParticles::BlockFreeSurfaceLagrangeNoSphereParticles(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurfaceLagrangeParticles(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {
	gradU = new double[(size_t)9 * ntot]();
}

BlockFreeSurfaceLagrangeNoSphereParticles::~BlockFreeSurfaceLagrangeNoSphereParticles() {
	delete[] gradU;
}

void BlockFreeSurfaceLagrangeNoSphereParticles::computeParticleFields() {
	computeGradP();
	computeGradU();
}

// see BlockFreeSurface::velocityGradient
void BlockFreeSurfaceLagrangeNoSphereParticles::computeGradU() {
	velocityGradient(gradU);
}
