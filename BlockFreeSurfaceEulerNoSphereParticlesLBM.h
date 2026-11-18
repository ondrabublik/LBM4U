/*
 * File:   BlockFreeSurfaceEulerNoSphereParticlesLBM.h
 *
 * Eulerian non-spherical particles transported by LBM (problemID 42): spheroids described by the number of
 * particles in the cell c, transported by the scalar populations g as problemID 4, and by the orientation
 * tensor A = <p p> of the cell (SpheroidOrientation, as problemID 33).
 *
 * Drag: K = <M>^-1 g(Re) with the mobility averaged over the orientations of the cell. The particle velocity
 * solves (m I + K) v_new = m (v - (v . grad) v) + K u + (m - rho V) g (problemID 4 with a tensor, implicit),
 * next to the walls it is limited as in problemID 4, the liquid gets c K (v - u). A sphere (r = 1) gives
 * problemID 4.
 *
 * Orientation transport: the net flux of particles through each link (population out minus population in) carries
 * the orientation of the upwind cell (as the flux of c A in problemID 33):
 *   S A_new = S A + sum_links max(in - out, 0) (A_neighbour - A),   S the sum of the populations after the streaming.
 * The measure (sum of the populations) A is conserved as the number of particles, A stays a convex combination
 * (trace 1, positive semidefinite) and a uniform A is kept exactly; A diffuses only by the particle diffusion and
 * the upwind error. Then A follows the Folgar-Tucker equation in every liquid cell.
 *
 * Refined mesh: the net flux exchanged with the other level (exchangeCorr of problemID 4) carries the orientation
 * of its upwind cell, collected in exchangeOrient and added with exchangeCorr.
 *
 * Parameters: those of problemID 4 and particleAspectRatio, orientationDiffusion, orientationClosure,
 * init/orientation<tag>.txt (SpheroidOrientation).
 */

#ifndef BLOCKFREESURFACEEULERNOSPHEREPARTICLESLBM_H
#define BLOCKFREESURFACEEULERNOSPHEREPARTICLESLBM_H

#include "BlockFreeSurfaceEulerParticlesLBM.h"
#include "SpheroidOrientation.h"

class BlockFreeSurfaceEulerNoSphereParticlesLBM : public BlockFreeSurfaceEulerParticlesLBM, public SpheroidOrientation {
public:
	// S A after the propagation (S the sum of the populations), the neighbours read the old orient until applyExchangeCorr
	double* orientNew;
	// (sum of the populations) A exchanged with the other level, added to the cell with exchangeCorr
	double* exchangeOrient;
	// released * A of a cell without liquid given to the liquid cells around
	double* releasedOrient;

	BlockFreeSurfaceEulerNoSphereParticlesLBM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par);
	~BlockFreeSurfaceEulerNoSphereParticlesLBM();

	void readInitConcentration(int id, string simulPath) override;
	void prepareParticleStep() override;
	void computeNewParticleVelocityField(double* distance, double* gradDistance) override;
	void propagateScalarCell(int ind) override;
	// new A of the liquid cells (transport and exchange with the other level), then the Folgar-Tucker step
	void applyExchangeCorr() override;
	bool releaseDryCells() override;
	void collectReleased() override;
	void prolongCell(int ind, const Block* src, int sInd) override;
	void restrictCell(int ind, Block* const* src, const int* sInd) override;
};

#endif /* BLOCKFREESURFACEEULERNOSPHEREPARTICLESLBM_H */
