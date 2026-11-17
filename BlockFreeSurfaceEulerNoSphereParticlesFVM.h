/*
 * File:   BlockFreeSurfaceEulerNoSphereParticlesFVM.h
 *
 * Eulerian non-spherical particles (problemID 33): spheroids (rotationally symmetric ellipsoids, aspect ratio
 * r = a / b, a along the symmetry axis p) described by the number of particles in the cell c (as problemID 3) and
 * the second order orientation tensor A = <p p> of the cell (symmetric, trace 1; A = I/3 random orientations,
 * A = p p all particles aligned).
 *
 * Drag: the mobility of a spheroid, M(p) = I / KT_perp + (1 / KT_par - 1 / KT_perp) p p (Stokes resistances
 * KT of LagrangeNoSphereParticles::resistance, times the viscosity), is linear in p p, so its average over the
 * orientations of the cell is exact:  M = I / KT_perp + (1 / KT_par - 1 / KT_perp) A.  The particle velocity
 * solves (m I + K)(v_new) = m v + K u + (m - rho V) g with K = M^-1 g(Re) (implicit, Schiller-Naumann g),
 * the liquid gets c K (v - u). A sphere (r = 1) gives problemID 3.
 *
 * Orientation: q = c A is transported with the same (limited) flux as c, the flux of q through a face is the
 * flux of c times A of the upwind cell (A stays a weighted mean: trace 1, positive semidefinite). Then
 * A follows the Folgar-Tucker equation in every liquid cell (RK4)
 *   dA/dt = W A - A W + lambda (D A + A D - 2 A4:D) + 2 C_I gamma (I - 3 A),   lambda = (r^2 - 1) / (r^2 + 1)
 * W, D vorticity and strain rate tensors of the liquid, gamma = sqrt(2 D:D), C_I the rotary diffusion
 * coefficient (Folgar & Tucker 1984), A4 the fourth order tensor from a closure (hybrid or quadratic).
 *
 * Parameters: particleDiameter (volume equivalent sphere), particleDensity, particleAspectRatio (r > 1 prolate,
 * r < 1 oblate, default 1), orientationDiffusion (C_I, default 0), orientationClosure (0 hybrid, 1 quadratic).
 * init/orientation<tag>.txt (optional): "A11 A22 A33 A12 A13 A23" per cell in the order of the block files,
 * default A = I/3.
 */

#ifndef BLOCKFREESURFACEEULERNOSPHEREPARTICLESFVM_H
#define BLOCKFREESURFACEEULERNOSPHEREPARTICLESFVM_H

#include "BlockFreeSurfaceEulerParticlesFVM.h"

class BlockFreeSurfaceEulerNoSphereParticlesFVM : public BlockFreeSurfaceEulerParticlesFVM {
public:
	// orientation tensor, orient[6 * ind + k], k = xx, yy, zz, xy, xz, yz
	double* orient;
	// velocity gradient du_a/dx_b, gradU[9 * ind + 3 * a + b] (lattice units of the level)
	double* gradU;
	double aspectRatio, lambdaJ, rotaryDiffusion;
	int closure; // 0 hybrid, 1 quadratic
	double KTpar, KTperp; // Stokes resistance per unit viscosity along and across the axis (lattice units)

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

	// drag tensor K for the orientation tensor A per unit viscosity, without the Schiller-Naumann correction
	void dragTensor(const double A[6], double K[9]) const;
	// one step dt of the Folgar-Tucker equation, grad = du_a/dx_b
	void orientationStep(double A[6], const double grad[9], double dt) const;
	// dA/dt of the Folgar-Tucker equation
	void orientationRate(const double A[6], const double grad[9], double dA[6]) const;

	// --- symmetric 3x3 tensors stored as xx, yy, zz, xy, xz, yz ---
	// eigenvalues (descending) and eigenvectors (columns of vec, vec[3 * row + col]), Jacobi rotations
	static void eigenSym(const double A[6], double val[3], double vec[9]);
	// trace 1 and positive semidefinite (negative eigenvalues are cut)
	static void projectOrientation(double A[6]);

protected:
	double* fluxOrient;     // flux of c A leaving the cell
	double* exchangeOrient; // c A coming from the finer level / collected by a ghost (as exchangeCorr)
	double* releasedOrient; // c A of a cell without liquid given to the liquid cells around (as released)
	void faceTransport(int ind, double fl, BlockFreeSurfaceEulerParticlesFVM* nb, int loc, bool toCoarser) override;
};

#endif /* BLOCKFREESURFACEEULERNOSPHEREPARTICLESFVM_H */
