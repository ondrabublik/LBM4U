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

// du_a/dx_b = 3 sum_i w_i e_ib u_a(x + e_i) (second order for a smooth field). The value behind a wall
// (bounce-back, the wall is half way to the neighbour) is the mirror image 2 u_wall - u, so the wall shear
// rate is exact for a linear profile; behind the free surface or an outlet the cell value is repeated
// (zero normal gradient, a stress free surface).
void BlockFreeSurfaceLagrangeNoSphereParticles::computeGradU() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		double* g = &gradU[9 * ind];
		for (int a = 0; a < 9; a++) {
			g[a] = 0;
		}
		if (flagMap[ind] != 1 && flagMap[ind] != 3) {
			continue;
		}
		double u0[3] = { U[ind], V[ind], W[ind] };
		for (int dir = 1; dir < nf; dir++) {
			int bl, loc;
			int ct = neighbourType(ind, dir, bl, loc);
			double un[3] = { u0[0], u0[1], u0[2] };
			if (ct == CT_FLUID && bl >= 0) {
				int fm = ((BlockFreeSurface*)blocks[bl])->flagMap[loc];
				if (fm == 1 || fm == 3) {
					un[0] = blocks[bl]->U[loc];
					un[1] = blocks[bl]->V[loc];
					un[2] = blocks[bl]->W[loc];
				}
			}
			else if (ct == CT_WALL || ct == CT_MOVING || ct == CT_INLET) {
				const double* uw = ct == CT_MOVING ? velWall : (ct == CT_INLET ? velInlet : nullptr);
				for (int a = 0; a < 3; a++) {
					un[a] = 2 * (uw ? uw[a] : 0.0) - u0[a];
				}
			}
			double e[3] = { (double)ex[dir], (double)ey[dir], (double)ez[dir] };
			for (int a = 0; a < 3; a++) {
				for (int b = 0; b < 3; b++) {
					g[3 * a + b] += 3 * w[dir] * e[b] * un[a];
				}
			}
		}
	}
}
