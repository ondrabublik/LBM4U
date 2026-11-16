/*
 * File:   LagrangeNoSphereParticles.h
 *
 * Non-spherical Lagrangian point particles (problemID 22): rigid ellipsoids with the semi-axes a1, a2, a3
 * along their body axes, an orientation (quaternion) and an angular velocity.
 *
 * Translation: the drag is a tensor, K = R diag(k1, k2, k3) R^T (R: body -> lab), with
 *   k_i = f_i mu KT_i g(Re),  KT_i = 16 pi a1 a2 a3 / (chi0 + a_i^2 alpha_i)    (Stokes, Oberbeck 1876)
 * where chi0 and alpha_i are the elliptic integrals of the ellipsoid (Carlson R_F, R_D), g(Re) is the
 * Schiller-Naumann correction with the volume equivalent diameter and f_i are optional user factors.
 * A sphere gives KT = 6 pi a, the model is then identical to problemID 2. The velocity is updated
 * implicitly (3x3 system), the reaction K (v - u) goes to the liquid as in problemID 2.
 *
 * Rotation: Euler equations in the body frame with the hydrodynamic torque of Jeffery (1922),
 *   T_1 = mu KR_1 (Omega_1 + lambda_1 E_23 - omega_1),  lambda_1 = (a2^2 - a3^2) / (a2^2 + a3^2)   (cyclic)
 *   KR_1 = 16 pi a1 a2 a3 (a2^2 + a3^2) / (3 (a2^2 alpha_2 + a3^2 alpha_3))
 * Omega is half of the fluid vorticity and E the strain rate at the particle (both in the body frame).
 * The update is implicit, so particles with a small moment of inertia follow the Jeffery orbits.
 * Limits: creeping flow torque (no inertial torque at finite Re), the torque is not returned to the liquid.
 *
 * init/particles.txt: the number of particles, then "x y z a1 a2 a3 density qw qx qy qz [f1 f2 f3]" (SI units,
 * the quaternion rotates the body axes to the lab frame, it is normalised).
 */

#ifndef LAGRANGENOSPHEREPARTICLES_H
#define LAGRANGENOSPHEREPARTICLES_H

#include <array>
#include "LagrangeParticles.h"

class LagrangeNoSphereParticles : public LagrangeParticles {
public:
	// per particle, lattice units of the level 0: semi-axes, user drag factors, quaternion (w, x, y, z),
	// angular velocity in the body frame, resistance coefficients (translation KT, rotation KR, without mu),
	// principal moments of inertia, fluid velocity gradient at the particle (du_a/dx_b at 3 a + b)
	vector<double> semiAxes, dragFactor, quat, omegaBody, KT, KR, inertia, gradF;

	LagrangeNoSphereParticles(Block** blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, const ParamMap& par);

	void readInitParticlesDistribution(string simulPath, double Ch, double Cu) override;
	void interpolateFluidVelocity() override;
	void newParticlePosition(double dt) override;

	// adds a particle (lattice units), used by readInitParticlesDistribution and the tests
	void setParticles(const vector<array<double, 14> >& p);

	// --- geometry and mechanics of an ellipsoid (lattice or any consistent units) ---
	static double carlsonRF(double x, double y, double z);
	static double carlsonRD(double x, double y, double z);
	// translational (KT) and rotational (KR) Stokes resistance per unit viscosity along the body axes
	static void resistance(const double a[3], double KT[3], double KR[3]);
	// rotation matrix body -> lab, R[3 * row + col], its columns are the body axes
	static void rotationMatrix(const double q[4], double R[9]);
	// one step of the angular velocity (body frame) and of the orientation, grad = du_a/dx_b (lab frame)
	static void rotationStep(double q[4], double wB[3], const double a[3], const double KR[3], const double I[3], double nu, const double grad[9], double dt);
	// Jeffery angular velocity (body frame) of a torque free ellipsoid
	static void jefferyOmega(const double q[4], const double a[3], const double grad[9], double wJ[3]);
	// lab angular velocity of the particle p
	void omegaLab(int p, double w[3]) const;

protected:
	double contactRadius(int p) const override;
	void writeExtraVTK(ofstream& f) override;
};

#endif /* LAGRANGENOSPHEREPARTICLES_H */
