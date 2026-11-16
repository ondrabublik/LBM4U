/*
 * File:   BlockFreeSurfaceLagrangeNoSphereParticles.h
 *
 * Free surface block coupled with non-spherical (ellipsoidal) Lagrangian particles (problemID 22,
 * LagrangeNoSphereParticles). Same coupling as BlockFreeSurfaceLagrangeParticles, in addition the
 * gradient of the fluid velocity, which drives the rotation of the particles (Jeffery).
 */

#ifndef BLOCKFREESURFACELAGRANGENOSPHEREPARTICLES_H
#define BLOCKFREESURFACELAGRANGENOSPHEREPARTICLES_H

#include "BlockFreeSurfaceLagrangeParticles.h"

class BlockFreeSurfaceLagrangeNoSphereParticles : public BlockFreeSurfaceLagrangeParticles {
public:
	// velocity gradient du_a/dx_b of the liquid cells, gradU[9 * ind + 3 * a + b], lattice units of the level
	double* gradU;

	BlockFreeSurfaceLagrangeNoSphereParticles(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par);
	~BlockFreeSurfaceLagrangeNoSphereParticles();

	void computeParticleFields() override;
	void computeGradU();
};

#endif /* BLOCKFREESURFACELAGRANGENOSPHEREPARTICLES_H */
