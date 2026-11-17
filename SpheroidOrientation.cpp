/*
 * File:   SpheroidOrientation.cpp
 *
 * Orientation tensor of the Eulerian non-spherical particles (problemID 33 and 42). See SpheroidOrientation.h.
 */

#include "SpheroidOrientation.h"
#include "LagrangeNoSphereParticles.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

using namespace std;

namespace {
	// symmetric tensor xx, yy, zz, xy, xz, yz <-> full matrix M[3 * row + col]
	void toFull(const double s[6], double M[9]) {
		M[0] = s[0]; M[1] = s[3]; M[2] = s[4];
		M[3] = s[3]; M[4] = s[1]; M[5] = s[5];
		M[6] = s[4]; M[7] = s[5]; M[8] = s[2];
	}

	void toSym(const double M[9], double s[6]) {
		s[0] = M[0];
		s[1] = M[4];
		s[2] = M[8];
		s[3] = (M[1] + M[3]) / 2;
		s[4] = (M[2] + M[6]) / 2;
		s[5] = (M[5] + M[7]) / 2;
	}

	void mul(const double A[9], const double B[9], double C[9]) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				C[3 * i + j] = A[3 * i] * B[j] + A[3 * i + 1] * B[3 + j] + A[3 * i + 2] * B[6 + j];
			}
		}
	}

	double det3(const double M[9]) {
		return M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) + M[2] * (M[3] * M[7] - M[4] * M[6]);
	}
}

const double SpheroidOrientation::ISO[6] = { 1.0 / 3, 1.0 / 3, 1.0 / 3, 0, 0, 0 };

SpheroidOrientation::~SpheroidOrientation() {
	delete[] orient;
	delete[] gradU;
}

void SpheroidOrientation::initOrientation(const ParamMap& par, double dParticle, int ntot, const string& model) {
	orient = new double[(size_t)6 * ntot];
	for (int i = 0; i < ntot; i++) {
		for (int k = 0; k < 6; k++) {
			orient[6 * i + k] = ISO[k];
		}
	}
	gradU = new double[(size_t)9 * ntot]();

	aspectRatio = Block::getParDouble("particleAspectRatio", par, false, 1.0);
	rotaryDiffusion = Block::getParDouble("orientationDiffusion", par, false, 0.0);
	closure = Block::getParInt("orientationClosure", par, false, 0);
	if (!(aspectRatio > 0) || rotaryDiffusion < 0 || closure < 0 || closure > 1) {
		throw invalid_argument(model + ": particleAspectRatio > 0, orientationDiffusion >= 0, orientationClosure 0 or 1 required");
	}
	lambdaJ = (aspectRatio * aspectRatio - 1) / (aspectRatio * aspectRatio + 1);
	// semi-axes of the spheroid with the volume of the sphere of diameter dParticle
	double b = dParticle / 2 / cbrt(aspectRatio);
	double a[3] = { aspectRatio * b, b, b };
	double KT[3], KR[3];
	LagrangeNoSphereParticles::resistance(a, KT, KR);
	KTpar = KT[0];
	KTperp = KT[1];
}

void SpheroidOrientation::readOrientation(const string& fileName, int nAll, const int* indAll) {
	ifstream myfile(fileName);
	if (!myfile.is_open()) {
		return;
	}
	string line;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		double A[6];
		if (!getline(myfile, line)) {
			throw runtime_error("Error in " + fileName + ": expected " + to_string(nAll) + " lines");
		}
		istringstream is(line);
		for (int k = 0; k < 6; k++) {
			if (!(is >> A[k])) {
				throw runtime_error("Error in " + fileName + ": expected 6 values A11 A22 A33 A12 A13 A23 on line " + to_string(i + 1));
			}
		}
		projectOrientation(A);
		for (int k = 0; k < 6; k++) {
			orient[6 * ind + k] = A[k];
		}
	}
}

// K = <M>^-1 per unit viscosity, <M> = I / KT_perp + (1 / KT_par - 1 / KT_perp) A
void SpheroidOrientation::dragTensor(const double A[6], double K[9]) const {
	double M[9];
	toFull(A, M);
	double dm = 1 / KTpar - 1 / KTperp;
	for (int i = 0; i < 9; i++) {
		M[i] *= dm;
	}
	M[0] += 1 / KTperp;
	M[4] += 1 / KTperp;
	M[8] += 1 / KTperp;
	inverse3(M, K);
}

void SpheroidOrientation::inverse3(const double M[9], double inv[9]) {
	double d = det3(M);
	inv[0] = (M[4] * M[8] - M[5] * M[7]) / d;
	inv[1] = (M[2] * M[7] - M[1] * M[8]) / d;
	inv[2] = (M[1] * M[5] - M[2] * M[4]) / d;
	inv[3] = (M[5] * M[6] - M[3] * M[8]) / d;
	inv[4] = (M[0] * M[8] - M[2] * M[6]) / d;
	inv[5] = (M[2] * M[3] - M[0] * M[5]) / d;
	inv[6] = (M[3] * M[7] - M[4] * M[6]) / d;
	inv[7] = (M[1] * M[6] - M[0] * M[7]) / d;
	inv[8] = (M[0] * M[4] - M[1] * M[3]) / d;
}

double SpheroidOrientation::dragCorrection(double Re) {
	return Re < 1000 ? 1 + 0.15 * pow(Re, 0.687) : 0.44 * Re / 24;
}

// ------------------------------------------------------------------------------------------------
// orientation dynamics
// ------------------------------------------------------------------------------------------------

void SpheroidOrientation::orientationRate(const double As[6], const double L[9], double dA[6]) const {
	double A[9], D[9], W[9];
	toFull(As, A);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			D[3 * i + j] = (L[3 * i + j] + L[3 * j + i]) / 2;
			W[3 * i + j] = (L[3 * i + j] - L[3 * j + i]) / 2;
		}
	}
	double WA[9], AW[9], DA[9], AD[9];
	mul(W, A, WA);
	mul(A, W, AW);
	mul(D, A, DA);
	mul(A, D, AD);
	double AdD = 0, DdD = 0;
	for (int k = 0; k < 9; k++) {
		AdD += A[k] * D[k];
		DdD += D[k] * D[k];
	}
	double trD = D[0] + D[4] + D[8];
	// closure A4:D, hybrid: f quadratic + (1 - f) linear, f = 1 - 27 det(A) (f = 1 for aligned, 0 for random particles)
	double f = 1;
	if (closure == 0) {
		f = min(max(1 - 27 * det3(A), 0.0), 1.0);
	}
	double out[9];
	double gamma = sqrt(2 * DdD);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			int k = 3 * i + j;
			double delta = i == j ? 1.0 : 0.0;
			double quad = A[k] * AdD;
			double lin = -(delta * trD + 2 * D[k]) / 35 + (A[k] * trD + 2 * (AD[k] + DA[k]) + delta * AdD) / 7;
			double A4D = f * quad + (1 - f) * lin;
			out[k] = WA[k] - AW[k] + lambdaJ * (DA[k] + AD[k] - 2 * A4D) + 2 * rotaryDiffusion * gamma * (delta - 3 * A[k]);
		}
	}
	toSym(out, dA);
}

// classical Runge-Kutta, sub-steps when the velocity gradient is large compared with 1 / dt
void SpheroidOrientation::orientationStep(double A[6], const double grad[9], double dt) const {
	double g2 = 0;
	for (int k = 0; k < 9; k++) {
		g2 += grad[k] * grad[k];
	}
	if (g2 == 0) {
		return;
	}
	double rate = sqrt(g2) * (1 + fabs(lambdaJ) + 6 * rotaryDiffusion);
	int nSub = max(1, (int)ceil(rate * dt / 0.05));
	double h = dt / nSub;
	for (int s = 0; s < nSub; s++) {
		double k1[6], k2[6], k3[6], k4[6], t[6];
		orientationRate(A, grad, k1);
		for (int k = 0; k < 6; k++) t[k] = A[k] + h / 2 * k1[k];
		orientationRate(t, grad, k2);
		for (int k = 0; k < 6; k++) t[k] = A[k] + h / 2 * k2[k];
		orientationRate(t, grad, k3);
		for (int k = 0; k < 6; k++) t[k] = A[k] + h * k3[k];
		orientationRate(t, grad, k4);
		for (int k = 0; k < 6; k++) {
			A[k] += h / 6 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]);
		}
	}
	projectOrientation(A);
}

// ------------------------------------------------------------------------------------------------
// symmetric tensors
// ------------------------------------------------------------------------------------------------

void SpheroidOrientation::eigenSym(const double As[6], double val[3], double vec[9]) {
	double a[9];
	toFull(As, a);
	for (int k = 0; k < 9; k++) {
		vec[k] = (k % 4 == 0) ? 1.0 : 0.0;
	}
	double scale = fabs(a[0]) + fabs(a[4]) + fabs(a[8]) + fabs(a[1]) + fabs(a[2]) + fabs(a[5]);
	for (int sweep = 0; sweep < 50; sweep++) {
		double off = fabs(a[1]) + fabs(a[2]) + fabs(a[5]);
		if (off <= 1e-17 * scale || off == 0) {
			break;
		}
		static const int P[3] = { 0, 0, 1 }, Q[3] = { 1, 2, 2 };
		for (int r = 0; r < 3; r++) {
			int p = P[r], q = Q[r];
			double apq = a[3 * p + q];
			if (apq == 0) {
				continue;
			}
			double theta = (a[3 * q + q] - a[3 * p + p]) / (2 * apq);
			double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1));
			double c = 1 / sqrt(t * t + 1), s = t * c;
			for (int k = 0; k < 3; k++) { // A J
				double akp = a[3 * k + p], akq = a[3 * k + q];
				a[3 * k + p] = c * akp - s * akq;
				a[3 * k + q] = s * akp + c * akq;
			}
			for (int k = 0; k < 3; k++) { // J^T A
				double apk = a[3 * p + k], aqk = a[3 * q + k];
				a[3 * p + k] = c * apk - s * aqk;
				a[3 * q + k] = s * apk + c * aqk;
			}
			for (int k = 0; k < 3; k++) { // V J
				double vkp = vec[3 * k + p], vkq = vec[3 * k + q];
				vec[3 * k + p] = c * vkp - s * vkq;
				vec[3 * k + q] = s * vkp + c * vkq;
			}
		}
	}
	int order[3] = { 0, 1, 2 };
	double d[3] = { a[0], a[4], a[8] };
	for (int i = 0; i < 2; i++) {
		for (int j = i + 1; j < 3; j++) {
			if (d[order[j]] > d[order[i]]) {
				swap(order[i], order[j]);
			}
		}
	}
	double v2[9];
	for (int c = 0; c < 3; c++) {
		val[c] = d[order[c]];
		// sign: the largest component of the eigenvector is positive
		int big = 0;
		for (int r = 1; r < 3; r++) {
			if (fabs(vec[3 * r + order[c]]) > fabs(vec[3 * big + order[c]])) {
				big = r;
			}
		}
		double sg = vec[3 * big + order[c]] < 0 ? -1.0 : 1.0;
		for (int r = 0; r < 3; r++) {
			v2[3 * r + c] = sg * vec[3 * r + order[c]];
		}
	}
	for (int k = 0; k < 9; k++) {
		vec[k] = v2[k];
	}
}

void SpheroidOrientation::projectOrientation(double A[6]) {
	double tr = A[0] + A[1] + A[2];
	if (!(tr > 0)) {
		for (int k = 0; k < 6; k++) {
			A[k] = ISO[k];
		}
		return;
	}
	for (int k = 0; k < 6; k++) {
		A[k] /= tr;
	}
	// positive semidefinite: all principal minors >= 0
	double M[9];
	toFull(A, M);
	const double tol = -1e-14; // round-off of an aligned (singular) tensor is not corrected
	if (A[0] >= tol && A[1] >= tol && A[2] >= tol && A[0] * A[1] - A[3] * A[3] >= tol && A[0] * A[2] - A[4] * A[4] >= tol
		&& A[1] * A[2] - A[5] * A[5] >= tol && det3(M) >= tol) {
		return;
	}
	double val[3], vec[9];
	eigenSym(A, val, vec);
	double sum = 0;
	for (int c = 0; c < 3; c++) {
		val[c] = max(val[c], 0.0);
		sum += val[c];
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			double m = 0;
			for (int c = 0; c < 3; c++) {
				m += vec[3 * i + c] * val[c] * vec[3 * j + c];
			}
			M[3 * i + j] = m / sum;
		}
	}
	toSym(M, A);
}
