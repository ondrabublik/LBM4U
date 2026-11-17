/*
 * File:   BlockFreeSurfaceEulerNoSphereParticlesFVM.h
 *
 * Eulerian non-spherical particles (problemID 33): spheroids (rotationally symmetric ellipsoids, aspect ratio
 * r = a / b, a along the symmetry axis p) described by the number of particles in the cell c (as problemID 3) and
 * the second order orientation tensor A = <p p> of the cell (symmetric, trace 1; A = I/3 random orientations,
 * A = p p all particles aligned).
 *
 * Drag: K = <M>^-1 g(Re) with the mobility <M> averaged over the orientations (SpheroidOrientation). The particle
 * velocity solves (m I + K)(v_new) = m v + K u + (m - rho V) g (implicit, Schiller-Naumann g), the liquid gets
 * c K (v - u). A sphere (r = 1) gives problemID 3.
 *
 * Orientation: q = c A is transported with the same (limited) flux as c, the flux of q through a face is the
 * flux of c times A of the upwind cell (A stays a weighted mean: trace 1, positive semidefinite). Then
 * A follows the Folgar-Tucker equation in every liquid cell (SpheroidOrientation, RK4).
 *
 * Parameters: particleDiameter (volume equivalent sphere), particleDensity and the orientation parameters of
 * SpheroidOrientation (particleAspectRatio, orientationDiffusion, orientationClosure, init/orientation<tag>.txt).
 */

#ifndef BLOCKFREESURFACEEULERNOSPHEREPARTICLESFVM_H
#define BLOCKFREESURFACEEULERNOSPHEREPARTICLESFVM_H

#include "BlockFreeSurfaceEulerParticlesFVM.h"
#include "SpheroidOrientation.h"

class BlockFreeSurfaceEulerNoSphereParticlesFVM : public BlockFreeSurfaceEulerParticlesFVM, public SpheroidOrientation {
public:
	BlockFreeSurfaceEulerNoSphereParticlesFVM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par);
	~BlockFreeSurfaceEulerNoSphereParticlesFVM();

	void readInitConcentration(int id, string simulPath) override;
	void prepareParticleStep() override;
	void computeNewParticleVelocityField() override;
	void applyFlux() override;
	void forwardExchange(int ind, BlockFreeSurfaceEulerParticlesFVM* dst, int dloc) override;
	bool releaseDryCells() override;
	void collectReleased() override;
	void prolongCell(int ind, const Block* src, int sInd) override;
	void restrictCell(int ind, Block* const* src, const int* sInd) override;

protected:
	double* fluxOrient;     // flux of c A leaving the cell
	double* exchangeOrient; // c A coming from the finer level / collected by a ghost (as exchangeCorr)
	double* releasedOrient; // c A of a cell without liquid given to the liquid cells around (as released)
	void faceTransport(int ind, double fl, BlockFreeSurfaceEulerParticlesFVM* nb, int loc, bool toCoarser) override;
};

#endif /* BLOCKFREESURFACEEULERNOSPHEREPARTICLESFVM_H */
