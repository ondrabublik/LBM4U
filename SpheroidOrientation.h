/*
 * File:   SpheroidOrientation.h
 *
 * Orientation of the Eulerian non-spherical particles (problemID 33 FVM, problemID 42 LBM), shared by both
 * models: spheroids (rotationally symmetric ellipsoids, aspect ratio r = a / b, a along the symmetry axis p)
 * described by the second order orientation tensor of the cell A = <p p> (symmetric, trace 1; A = I/3 random
 * orientations, A = p p all particles aligned).
 *
 * Drag: the mobility of a spheroid, M(p) = I / KT_perp + (1 / KT_par - 1 / KT_perp) p p (Stokes resistances
 * KT of LagrangeNoSphereParticles::resistance, times the viscosity), is linear in p p, so its average over the
 * orientations of the cell is exact:  M = I / KT_perp + (1 / KT_par - 1 / KT_perp) A,  K = M^-1.
 *
 * Orientation: A follows the Folgar-Tucker equation in every liquid cell (RK4)
 *   dA/dt = W A - A W + lambda (D A + A D - 2 A4:D) + 2 C_I gamma (I - 3 A),   lambda = (r^2 - 1) / (r^2 + 1)
 * W, D vorticity and strain rate tensors of the liquid, gamma = sqrt(2 D:D), C_I the rotary diffusion
 * coefficient (Folgar & Tucker 1984), A4 the fourth order tensor from a closure (hybrid or quadratic).
 *
 * Parameters: particleAspectRatio (r > 1 prolate, r < 1 oblate, default 1), orientationDiffusion (C_I,
 * default 0), orientationClosure (0 hybrid, 1 quadratic). init/orientation<tag>.txt (optional):
 * "A11 A22 A33 A12 A13 A23" per cell in the order of the block files, default A = I/3.
 */

#ifndef SPHEROIDORIENTATION_H
#define SPHEROIDORIENTATION_H

#include <string>
#include "Block.h"

class SpheroidOrientation {
public:
	// orientation tensor, orient[6 * ind + k], k = xx, yy, zz, xy, xz, yz
	double* orient = nullptr;
	// velocity gradient du_a/dx_b, gradU[9 * ind + 3 * a + b] (lattice units of the level)
	double* gradU = nullptr;
	double aspectRatio = 1, lambdaJ = 0, rotaryDiffusion = 0;
	int closure = 0;                 // 0 hybrid, 1 quadratic
	double KTpar = 0, KTperp = 0;    // Stokes resistance per unit viscosity along and across the axis (lattice units)

	static const double ISO[6];      // random orientation I/3

	SpheroidOrientation() {}
	~SpheroidOrientation();
	SpheroidOrientation(const SpheroidOrientation&) = delete;
	SpheroidOrientation& operator=(const SpheroidOrientation&) = delete;

	// parameters of the particles (dParticle: diameter of the volume equivalent sphere, lattice units),
	// A = I/3 in all ntot cells
	void initOrientation(const ParamMap& par, double dParticle, int ntot, const std::string& model);
	// init/orientation<tag>.txt of the block, nothing without the file
	void readOrientation(const std::string& fileName, int nAll, const int* indAll);

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
	// inverse of a regular 3x3 matrix M[3 * row + col]
	static void inverse3(const double M[9], double inv[9]);
	// Schiller-Naumann correction of the Stokes drag, Cd = 0.44 above Re = 1000 (as the spherical models)
	static double dragCorrection(double Re);
};

#endif /* SPHEROIDORIENTATION_H */
