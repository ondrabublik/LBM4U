/*
 * File:   LagrangeNoSphereParticles.cpp
 *
 * Non-spherical (ellipsoidal) Lagrangian particles, problemID 22. See LagrangeNoSphereParticles.h.
 */

#include "LagrangeNoSphereParticles.h"
#include "BlockFreeSurfaceLagrangeNoSphereParticles.h"
#include <cmath>
#include <stdexcept>

using namespace std;

namespace {
	const double PI = 3.14159265358979323846;

	// Schiller-Naumann correction of the Stokes drag, Cd = 0.44 above Re = 1000 (as LagrangeParticles::cd)
	double dragCorrection(double Re) {
		return Re < 1000 ? 1 + 0.15 * pow(Re, 0.687) : 0.44 * Re / 24;
	}

	void quatMultiply(const double a[4], const double b[4], double out[4]) {
		double r[4] = {
			a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
			a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
			a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
			a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0] };
		for (int i = 0; i < 4; i++) {
			out[i] = r[i];
		}
	}

	// solves A x = b for a symmetric positive definite 3x3 matrix (Cramer's rule)
	void solve3(const double A[9], const double b[3], double x[3]) {
		double c00 = A[4] * A[8] - A[5] * A[7];
		double c01 = A[5] * A[6] - A[3] * A[8];
		double c02 = A[3] * A[7] - A[4] * A[6];
		double det = A[0] * c00 + A[1] * c01 + A[2] * c02;
		double inv[9] = {
			c00, A[2] * A[7] - A[1] * A[8], A[1] * A[5] - A[2] * A[4],
			c01, A[0] * A[8] - A[2] * A[6], A[2] * A[3] - A[0] * A[5],
			c02, A[1] * A[6] - A[0] * A[7], A[0] * A[4] - A[1] * A[3] };
		for (int i = 0; i < 3; i++) {
			x[i] = (inv[3 * i] * b[0] + inv[3 * i + 1] * b[1] + inv[3 * i + 2] * b[2]) / det;
		}
	}
}

LagrangeNoSphereParticles::LagrangeNoSphereParticles(Block** blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, const ParamMap& par)
	: LagrangeParticles(blocks, nBlock, nBlockMax, nBi, nBj, nBk, ni, nj, nk, par) {
}

// ------------------------------------------------------------------------------------------------
// ellipsoid geometry
// ------------------------------------------------------------------------------------------------

// Carlson's symmetric elliptic integral R_F (duplication algorithm, Numerical Recipes 3rd ed.)
double LagrangeNoSphereParticles::carlsonRF(double x, double y, double z) {
	const double ERRTOL = 0.0025, C1 = 1.0 / 24, C2 = 0.1, C3 = 3.0 / 44, C4 = 1.0 / 14;
	double ave, delx, dely, delz;
	for (int it = 0; it < 100; it++) {
		double sx = sqrt(x), sy = sqrt(y), sz = sqrt(z);
		double alamb = sx * (sy + sz) + sy * sz;
		x = 0.25 * (x + alamb);
		y = 0.25 * (y + alamb);
		z = 0.25 * (z + alamb);
		ave = (x + y + z) / 3;
		delx = (ave - x) / ave;
		dely = (ave - y) / ave;
		delz = (ave - z) / ave;
		if (max(max(fabs(delx), fabs(dely)), fabs(delz)) <= ERRTOL) {
			break;
		}
	}
	double e2 = delx * dely - delz * delz;
	double e3 = delx * dely * delz;
	return (1 + (C1 * e2 - C2 - C3 * e3) * e2 + C4 * e3) / sqrt(ave);
}

// Carlson's symmetric elliptic integral R_D (Numerical Recipes 3rd ed.)
double LagrangeNoSphereParticles::carlsonRD(double x, double y, double z) {
	const double ERRTOL = 0.0015, C1 = 3.0 / 14, C2 = 1.0 / 6, C3 = 9.0 / 22, C4 = 3.0 / 26, C5 = 0.25 * C3, C6 = 1.5 * C4;
	double sum = 0, fac = 1, ave, delx, dely, delz;
	for (int it = 0; it < 100; it++) {
		double sx = sqrt(x), sy = sqrt(y), sz = sqrt(z);
		double alamb = sx * (sy + sz) + sy * sz;
		sum += fac / (sz * (z + alamb));
		fac *= 0.25;
		x = 0.25 * (x + alamb);
		y = 0.25 * (y + alamb);
		z = 0.25 * (z + alamb);
		ave = 0.2 * (x + y + 3 * z);
		delx = (ave - x) / ave;
		dely = (ave - y) / ave;
		delz = (ave - z) / ave;
		if (max(max(fabs(delx), fabs(dely)), fabs(delz)) <= ERRTOL) {
			break;
		}
	}
	double ea = delx * dely, eb = delz * delz, ec = ea - eb, ed = ea - 6 * eb, ee = ed + ec + ec;
	return 3 * sum + fac * (1 + ed * (-C1 + C5 * ed - C6 * delz * ee) + delz * (C2 * ee + delz * (-C3 * ec + delz * C4 * ea))) / (ave * sqrt(ave));
}

// chi0 = abc int_0^inf dl / D,  alpha_i = abc int_0^inf dl / ((a_i^2 + l) D),  D = sqrt((a^2 + l)(b^2 + l)(c^2 + l))
// int dl / D = 2 R_F(a^2, b^2, c^2),  int dl / ((a_i^2 + l) D) = 2/3 R_D(a_j^2, a_k^2, a_i^2)
void LagrangeNoSphereParticles::resistance(const double a[3], double KT[3], double KR[3]) {
	double abc = a[0] * a[1] * a[2];
	double s[3] = { a[0] * a[0], a[1] * a[1], a[2] * a[2] };
	double chi0 = abc * 2 * carlsonRF(s[0], s[1], s[2]);
	double alpha[3];
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3, k = (i + 2) % 3;
		alpha[i] = abc * 2.0 / 3.0 * carlsonRD(s[j], s[k], s[i]);
	}
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3, k = (i + 2) % 3;
		KT[i] = 16 * PI * abc / (chi0 + s[i] * alpha[i]);
		KR[i] = 16 * PI * abc * (s[j] + s[k]) / (3 * (s[j] * alpha[j] + s[k] * alpha[k]));
	}
}

void LagrangeNoSphereParticles::rotationMatrix(const double q[4], double R[9]) {
	double w = q[0], x = q[1], y = q[2], z = q[3];
	R[0] = 1 - 2 * (y * y + z * z); R[1] = 2 * (x * y - w * z);     R[2] = 2 * (x * z + w * y);
	R[3] = 2 * (x * y + w * z);     R[4] = 1 - 2 * (x * x + z * z); R[5] = 2 * (y * z - w * x);
	R[6] = 2 * (x * z - w * y);     R[7] = 2 * (y * z + w * x);     R[8] = 1 - 2 * (x * x + y * y);
}

// Omega' + lambda E' in the body frame (zero torque of Jeffery's creeping flow solution)
void LagrangeNoSphereParticles::jefferyOmega(const double q[4], const double a[3], const double g[9], double wJ[3]) {
	double R[9];
	rotationMatrix(q, R);
	// half vorticity and strain rate in the lab frame
	double Om[3] = { 0.5 * (g[7] - g[5]), 0.5 * (g[2] - g[6]), 0.5 * (g[3] - g[1]) };
	double E[9];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			E[3 * i + j] = 0.5 * (g[3 * i + j] + g[3 * j + i]);
		}
	}
	// body frame: Om' = R^T Om, E' = R^T E R
	double OmB[3], EB[9];
	for (int i = 0; i < 3; i++) {
		OmB[i] = R[i] * Om[0] + R[3 + i] * Om[1] + R[6 + i] * Om[2];
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			double s = 0;
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					s += R[3 * k + i] * E[3 * k + l] * R[3 * l + j];
				}
			}
			EB[3 * i + j] = s;
		}
	}
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3, k = (i + 2) % 3;
		double sj = a[j] * a[j], sk = a[k] * a[k];
		wJ[i] = OmB[i] + (sj - sk) / (sj + sk) * EB[3 * j + k];
	}
}

// I dw/dt + w x (I w) = nu KR (wJ - w) in the body frame, implicit in the torque, then q <- q exp(w dt / 2)
void LagrangeNoSphereParticles::rotationStep(double q[4], double wB[3], const double a[3], const double KR[3], const double I[3], double nu, const double grad[9], double dt) {
	double wJ[3] = { 0, 0, 0 };
	if (nu > 0) {
		jefferyOmega(q, a, grad, wJ);
	}
	double gyro[3] = {
		(I[2] - I[1]) * wB[1] * wB[2],
		(I[0] - I[2]) * wB[2] * wB[0],
		(I[1] - I[0]) * wB[0] * wB[1] };
	for (int i = 0; i < 3; i++) {
		double c = nu * KR[i];
		double den = I[i] / dt + c;
		if (den > 0) {
			wB[i] = (I[i] / dt * wB[i] - gyro[i] + c * wJ[i]) / den;
		}
	}
	double w = sqrt(wB[0] * wB[0] + wB[1] * wB[1] + wB[2] * wB[2]);
	double angle = w * dt;
	if (angle > 0) {
		double s = sin(angle / 2) / w;
		double dq[4] = { cos(angle / 2), wB[0] * s, wB[1] * s, wB[2] * s };
		quatMultiply(q, dq, q);
	}
	double n = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
	for (int i = 0; i < 4; i++) {
		q[i] /= n;
	}
}

void LagrangeNoSphereParticles::omegaLab(int p, double w[3]) const {
	double R[9];
	rotationMatrix(&quat[4 * p], R);
	const double* wb = &omegaBody[3 * p];
	for (int i = 0; i < 3; i++) {
		w[i] = R[3 * i] * wb[0] + R[3 * i + 1] * wb[1] + R[3 * i + 2] * wb[2];
	}
}

// ------------------------------------------------------------------------------------------------
// particles
// ------------------------------------------------------------------------------------------------

// p: x y z a1 a2 a3 rho qw qx qy qz f1 f2 f3 in lattice units (rho relative to the fluid)
void LagrangeNoSphereParticles::setParticles(const vector<array<double, 14> >& p) {
	int n = (int)p.size();
	resize(n);
	for (vector<double>* v : { &semiAxes, &dragFactor, &omegaBody, &KT, &KR, &inertia }) {
		v->assign(3 * (size_t)n, 0.0);
	}
	quat.assign(4 * (size_t)n, 0.0);
	gradF.assign(9 * (size_t)n, 0.0);
	dMax = 0;
	for (int i = 0; i < n; i++) {
		const array<double, 14>& r = p[i];
		partX[i] = r[0];
		partY[i] = r[1];
		partZ[i] = r[2];
		double* a = &semiAxes[3 * i];
		for (int d = 0; d < 3; d++) {
			a[d] = r[3 + d];
			dragFactor[3 * i + d] = r[11 + d];
		}
		partRho[i] = r[6];
		double qn = sqrt(r[7] * r[7] + r[8] * r[8] + r[9] * r[9] + r[10] * r[10]);
		double* q = &quat[4 * i];
		if (qn > 0) {
			for (int d = 0; d < 4; d++) {
				q[d] = r[7 + d] / qn;
			}
		}
		else {
			q[0] = 1;
		}
		partVolume[i] = 4.0 / 3.0 * PI * a[0] * a[1] * a[2];
		partD[i] = 2 * cbrt(a[0] * a[1] * a[2]); // volume equivalent diameter
		partCrossection[i] = PI / 4 * partD[i] * partD[i];
		partMass[i] = partVolume[i] * partRho[i];
		double* I = &inertia[3 * i];
		I[0] = partMass[i] * (a[1] * a[1] + a[2] * a[2]) / 5;
		I[1] = partMass[i] * (a[2] * a[2] + a[0] * a[0]) / 5;
		I[2] = partMass[i] * (a[0] * a[0] + a[1] * a[1]) / 5;
		resistance(a, &KT[3 * i], &KR[3 * i]);
		dMax = max(dMax, partD[i]);
	}
}

// init/particles.txt: the number of particles, then "x y z a1 a2 a3 density qw qx qy qz [f1 f2 f3]" (SI units)
void LagrangeNoSphereParticles::readInitParticlesDistribution(string simulPath, double Ch, double Cu) {
	string fileName = simulPath + "init/particles.txt";
	ifstream myfile(fileName);
	if (!myfile.is_open()) {
		throw runtime_error("file " + fileName + " with the particles not found");
	}
	string line;
	int n = 0;
	if (!getline(myfile, line) || !(istringstream(line) >> n) || n < 0) {
		throw runtime_error("Error in " + fileName + ": the first line has to hold the number of particles");
	}
	vector<array<double, 14> > p(n);
	for (int i = 0; i < n; i++) {
		array<double, 14>& r = p[i];
		if (!getline(myfile, line)) {
			throw runtime_error("Error in " + fileName + ": expected " + to_string(n) + " particles");
		}
		istringstream iss(line);
		for (int k = 0; k < 11; k++) {
			if (!(iss >> r[k])) {
				throw runtime_error("Error in " + fileName + ": particle " + to_string(i + 1) + " needs x y z a1 a2 a3 density qw qx qy qz [f1 f2 f3]");
			}
		}
		for (int k = 11; k < 14; k++) {
			if (!(iss >> r[k])) {
				r[k] = 1;
			}
		}
		for (int k = 0; k < 6; k++) {
			r[k] /= Ch;
		}
		r[6] /= Crho;
		if (!(r[3] > 0 && r[4] > 0 && r[5] > 0 && r[6] > 0 && r[11] > 0 && r[12] > 0 && r[13] > 0)) {
			throw runtime_error("Error in " + fileName + ": particle " + to_string(i + 1) + " needs positive semi-axes, density and drag factors");
		}
	}
	setParticles(p);
	cout << nParticles << " non-spherical Lagrangian particles (ellipsoids) were read" << endl;
}

// fluid velocity as for the spheres and in addition the velocity gradient (per cell of the level 0)
void LagrangeNoSphereParticles::interpolateFluidVelocity() {
	LagrangeParticles::interpolateFluidVelocity();
#pragma omp parallel for
	for (int p = 0; p < nParticles; p++) {
		double* g = &gradF[9 * p];
		for (int a = 0; a < 9; a++) {
			g[a] = 0;
		}
		Cell c;
		if (inGas[p] || !locate(partX[p], partY[p], partZ[p], c)) {
			continue;
		}
		double sum = 0;
		forLiquidCells(c, [&](BlockFreeSurfaceLagrangeParticles* b, int ind, double weight) {
			const double* gc = &((BlockFreeSurfaceLagrangeNoSphereParticles*)b)->gradU[9 * ind];
			for (int a = 0; a < 9; a++) {
				g[a] += weight * gc[a];
			}
			sum += weight;
		});
		if (sum > 0) {
			double s = (double)(1 << c.level) / sum; // gradient per cell of the level 0
			for (int a = 0; a < 9; a++) {
				g[a] *= s;
			}
		}
	}
}

void LagrangeNoSphereParticles::newParticlePosition(double dtP) {
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nParticles; i++) {
		double x0 = partX[i], y0 = partY[i], z0 = partZ[i];
		double L0 = wallDistance(i, x0, y0, z0);
		double n0[3] = { partWallNormalX[i], partWallNormalY[i], partWallNormalZ[i] };
		double* q = &quat[4 * i];
		double* wB = &omegaBody[3 * i];
		if (!inGas[i]) {
			// the fluid does not flow into the wall
			if (L0 >= 0 && L0 < 1) {
				double norm = n0[0] * veloFX[i] + n0[1] * veloFY[i] + n0[2] * veloFZ[i];
				if (norm < 0) {
					veloFX[i] -= norm * n0[0];
					veloFY[i] -= norm * n0[1];
					veloFZ[i] -= norm * n0[2];
				}
			}
			double v[3] = { veloX[i], veloY[i], veloZ[i] };
			double u[3] = { veloFX[i], veloFY[i], veloFZ[i] };
			double slip = sqrt((v[0] - u[0]) * (v[0] - u[0]) + (v[1] - u[1]) * (v[1] - u[1]) + (v[2] - u[2]) * (v[2] - u[2]));
			double corr = dragCorrection(partD[i] * slip / vLB);
			// drag tensor K = R diag(k) R^T
			double R[9], K[9];
			rotationMatrix(q, R);
			double k[3];
			for (int d = 0; d < 3; d++) {
				k[d] = dragFactor[3 * i + d] * vLB * KT[3 * i + d] * corr;
			}
			for (int r = 0; r < 3; r++) {
				for (int s = 0; s < 3; s++) {
					K[3 * r + s] = R[3 * r] * k[0] * R[3 * s] + R[3 * r + 1] * k[1] * R[3 * s + 1] + R[3 * r + 2] * k[2] * R[3 * s + 2];
				}
			}
			double m = partMass[i];
			double buoy = m - partVolume[i];
			double Fint[3] = { FinterX[i], FinterY[i], FinterZ[i] };
			double A[9], rhs[3];
			for (int r = 0; r < 3; r++) {
				rhs[r] = m / dtP * v[r] + G[r] * buoy + Fint[r];
				for (int s = 0; s < 3; s++) {
					A[3 * r + s] = K[3 * r + s] + (r == s ? m / dtP : 0.0);
					rhs[r] += K[3 * r + s] * u[s];
				}
			}
			double vn[3];
			solve3(A, rhs, vn);
			veloX[i] = vn[0];
			veloY[i] = vn[1];
			veloZ[i] = vn[2];
			// reaction of the drag: impulse given to the fluid, summed over the sub-steps of the LBM step
			double w[3] = { vn[0] - u[0], vn[1] - u[1], vn[2] - u[2] };
			Gx[i] += dtP * (K[0] * w[0] + K[1] * w[1] + K[2] * w[2]);
			Gy[i] += dtP * (K[3] * w[0] + K[4] * w[1] + K[5] * w[2]);
			Gz[i] += dtP * (K[6] * w[0] + K[7] * w[1] + K[8] * w[2]);
			rotationStep(q, wB, &semiAxes[3 * i], &KR[3 * i], &inertia[3 * i], vLB, &gradF[9 * i], dtP);
		}
		else {
			veloX[i] += dtP * G[0];
			veloY[i] += dtP * G[1];
			veloZ[i] += dtP * G[2];
			rotationStep(q, wB, &semiAxes[3 * i], &KR[3 * i], &inertia[3 * i], 0.0, &gradF[9 * i], dtP); // torque free
		}
		FinterX[i] = 0;
		FinterY[i] = 0;
		FinterZ[i] = 0;
		nInter[i] = 1;

		partX[i] += dtP * veloX[i];
		partY[i] += dtP * veloY[i];
		partZ[i] += dtP * veloZ[i];

		wallContact(i, x0, y0, z0, L0, n0);

		if (grid) {
			double size[3] = { (double)nBi * ni, (double)nBj * nj, (double)nBk * nk };
			double* pos[3] = { &partX[i], &partY[i], &partZ[i] };
			for (int d = 0; d < 3; d++) {
				if (grid->periodic[d] && std::isfinite(*pos[d])) {
					*pos[d] -= size[d] * floor(*pos[d] / size[d]);
				}
			}
		}
	}
}

// support function of the ellipsoid along the wall normal: sqrt(n^T R diag(a^2) R^T n)
double LagrangeNoSphereParticles::contactRadius(int p) const {
	double R[9];
	rotationMatrix(&quat[4 * p], R);
	double n[3] = { partWallNormalX[p], partWallNormalY[p], partWallNormalZ[p] };
	double s = 0;
	for (int d = 0; d < 3; d++) {
		double nb = R[d] * n[0] + R[3 + d] * n[1] + R[6 + d] * n[2]; // body component
		s += semiAxes[3 * p + d] * semiAxes[3 * p + d] * nb * nb;
	}
	double nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
	return nn > 0 ? sqrt(s / nn) : partD[p] / 2;
}

// shape tensor R diag(a) R^T (cells of the level 0, for the Tensor Glyph filter of ParaView), body axes and
// the angular velocity [rad/s]
void LagrangeNoSphereParticles::writeExtraVTK(ofstream& f) {
	vector<double> R(9 * (size_t)nParticles);
	for (int p = 0; p < nParticles; p++) {
		rotationMatrix(&quat[4 * p], &R[9 * p]);
	}
	f << "TENSORS shape float\n";
	for (int p = 0; p < nParticles; p++) {
		const double* r = &R[9 * p];
		const double* a = &semiAxes[3 * p];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				f << r[3 * i] * a[0] * r[3 * j] + r[3 * i + 1] * a[1] * r[3 * j + 1] + r[3 * i + 2] * a[2] * r[3 * j + 2] << (j < 2 ? " " : "\n");
			}
		}
	}
	f << "\n";
	for (int ax = 0; ax < 3; ax++) {
		f << "VECTORS axis" << ax + 1 << " float\n";
		for (int p = 0; p < nParticles; p++) {
			const double* r = &R[9 * p];
			f << r[ax] << " " << r[3 + ax] << " " << r[6 + ax] << "\n";
		}
		f << "\n";
	}
	f << "VECTORS angularVelocity float\n";
	for (int p = 0; p < nParticles; p++) {
		double w[3];
		omegaLab(p, w);
		f << w[0] / Ct << " " << w[1] / Ct << " " << w[2] / Ct << "\n";
	}
	f << "\n";
}
