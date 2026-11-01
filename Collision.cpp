/*
 * File:   Collision.cpp
 *
 * BGK, TRT and MRT collision operators.
 */

#include "Collision.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

using namespace std;

bool invertMatrix(const vector<double>& A, vector<double>& inv, int n) {
	vector<double> a(A);
	inv.assign(n * n, 0.0);
	for (int i = 0; i < n; i++) {
		inv[n * i + i] = 1.0;
	}
	for (int c = 0; c < n; c++) {
		int piv = c;
		for (int r = c + 1; r < n; r++) {
			if (fabs(a[n * r + c]) > fabs(a[n * piv + c])) {
				piv = r;
			}
		}
		if (fabs(a[n * piv + c]) < 1e-14) {
			return false;
		}
		if (piv != c) {
			for (int k = 0; k < n; k++) {
				swap(a[n * c + k], a[n * piv + k]);
				swap(inv[n * c + k], inv[n * piv + k]);
			}
		}
		double d = a[n * c + c];
		for (int k = 0; k < n; k++) {
			a[n * c + k] /= d;
			inv[n * c + k] /= d;
		}
		for (int r = 0; r < n; r++) {
			if (r != c) {
				double fct = a[n * r + c];
				if (fct != 0.0) {
					for (int k = 0; k < n; k++) {
						a[n * r + k] -= fct * a[n * c + k];
						inv[n * r + k] -= fct * inv[n * c + k];
					}
				}
			}
		}
	}
	return true;
}

CollisionOperator::CollisionOperator() {
	model = COLLISION_BGK;
	nf = 0;
	opp = nullptr;
	magic = 3.0 / 16.0;
	tauMinusMax = 1e30;
	mrtMagicSq = true;
	cachedTau = -1;
}

void CollisionOperator::init(int model, int nf, const int* ex, const int* ey, const int* ez, const int* opp, double magic, bool mrtMagicSq, double tau0, double tauMinusMax) {
	this->tauMinusMax = tauMinusMax;
	this->model = model;
	this->nf = nf;
	this->opp = opp;
	this->magic = magic;
	this->mrtMagicSq = mrtMagicSq;
	if (model == COLLISION_MRT) {
		buildMatrices(ex, ey, ez);
		buildT(tau0, cachedT);
		cachedTau = tau0;
	}
}

void CollisionOperator::buildMatrices(const int* ex, const int* ey, const int* ez) {
	M.assign(nf * nf, 0.0);
	S0.assign(nf, 0.0);
	stressMoments.clear();
	fluxMoments.clear();
	if (nf == 19) {
		// d'Humieres et al. (2002), velocity ordering of Mesh::model(319)
		static const double m19[19 * 19] = {
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			-30, -11, -11, -11, -11, -11, -11, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
			12, -4, -4, -4, -4, -4, -4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0,
			0, -4, 4, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0,
			0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1,
			0, 0, 0, -4, 4, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1,
			0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1,
			0, 0, 0, 0, 0, -4, 4, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1,
			0, 2, 2, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -2, -2, -2, -2,
			0, -4, -4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, -2, -2, -2, -2,
			0, 0, 0, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0,
			0, 0, 0, -2, -2, 2, 2, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, -1, -1, 1, 1, 0, 0, 0, 0, 1, -1, 1, -1,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, -1, -1, 1, 1
		};
		M.assign(m19, m19 + 19 * 19);
		double s[19] = { 0, 1.19, 1.4, 0, 1.2, 0, 1.2, 0, 1.2, 0, 1.4, 0, 1.4, 0, 0, 0, 1.98, 1.98, 1.98 };
		S0.assign(s, s + 19);
		stressMoments = { 9, 11, 13, 14, 15 };
		fluxMoments = { 4, 6, 8, 16, 17, 18 }; // all odd non-conserved moments
	}
	else if (nf == 15) {
		static const double m15[15 * 15] = {
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			-2, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1,
			16, -4, -4, -4, -4, -4, -4, 1, 1, 1, 1, 1, 1, 1, 1,
			0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1,
			0, -4, 4, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1,
			0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1,
			0, 0, 0, -4, 4, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1,
			0, 0, 0, 0, 0, 1, -1, 1, 1, 1, 1, -1, -1, -1, -1,
			0, 0, 0, 0, 0, -4, 4, 1, 1, 1, 1, -1, -1, -1, -1,
			0, 2, 2, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1,
			0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, -1, -1, 1, 1,
			0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1,
			0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 1, -1
		};
		M.assign(m15, m15 + 15 * 15);
		double s[15] = { 0, 1.19, 1.4, 0, 1.2, 0, 1.2, 0, 1.2, 0, 0, 0, 0, 0, 1.98 };
		S0.assign(s, s + 15);
		stressMoments = { 9, 10, 11, 12, 13 };
		fluxMoments = { 4, 6, 8, 14 }; // all odd non-conserved moments
	}
	else if (nf == 9) {
		// Lallemand & Luo (2000), rows rho, e, eps, jx, qx, jy, qy, pxx, pxy
		for (int p = 0; p < 9; p++) {
			int c2 = ex[p] * ex[p] + ey[p] * ey[p];
			M[9 * 0 + p] = 1;
			M[9 * 1 + p] = -4 + 3 * c2;
			M[9 * 2 + p] = 4 - 10.5 * c2 + 4.5 * c2 * c2;
			M[9 * 3 + p] = ex[p];
			M[9 * 4 + p] = (-5 + 3 * c2) * ex[p];
			M[9 * 5 + p] = ey[p];
			M[9 * 6 + p] = (-5 + 3 * c2) * ey[p];
			M[9 * 7 + p] = ex[p] * ex[p] - ey[p] * ey[p];
			M[9 * 8 + p] = ex[p] * ey[p];
		}
		double s[9] = { 0, 1.19, 1.4, 0, 1.2, 0, 1.2, 0, 0 };
		S0.assign(s, s + 9);
		stressMoments = { 7, 8 };
		fluxMoments = { 4, 6 };
	}
	else {
		throw invalid_argument("MRT collision is not defined for this velocity model");
	}
	if (!invertMatrix(M, iM, nf)) {
		throw runtime_error("MRT matrix is singular");
	}
}

void CollisionOperator::diagonal(double tau, double* S) const {
	double sNu = 1.0 / tau;
	for (int i = 0; i < nf; i++) {
		S[i] = S0[i];
	}
	for (int m : stressMoments) {
		S[m] = sNu;
	}
	if (mrtMagicSq) {
		double sq = max(8.0 * (2.0 - sNu) / (8.0 - sNu), 1.0 / tauMinusMax);
		for (int m : fluxMoments) {
			S[m] = sq;
		}
	}
}

void CollisionOperator::buildT(double tau, vector<double>& T) const {
	vector<double> S(nf);
	diagonal(tau, S.data());
	T.assign(nf * nf, 0.0);
	for (int i = 0; i < nf; i++) {
		for (int j = 0; j < nf; j++) {
			double sum = 0;
			for (int k = 0; k < nf; k++) {
				sum += iM[nf * i + k] * S[k] * M[nf * k + j];
			}
			T[nf * i + j] = sum;
		}
	}
}

double CollisionOperator::tauMinus(double tau) const {
	return min(0.5 + magic / (tau - 0.5), tauMinusMax);
}

void CollisionOperator::relax(const double* fneq, double tau, double* delta) const {
	switch (model) {
	case COLLISION_BGK: {
		double omega = 1.0 / tau;
		for (int p = 0; p < nf; p++) {
			delta[p] = omega * fneq[p];
		}
		break;
	}
	case COLLISION_TRT: {
		double omegaP = 1.0 / tau;
		double omegaM = 1.0 / tauMinus(tau);
		for (int p = 0; p < nf; p++) {
			double sym = 0.5 * (fneq[p] + fneq[opp[p]]);
			double asym = 0.5 * (fneq[p] - fneq[opp[p]]);
			delta[p] = omegaP * sym + omegaM * asym;
		}
		break;
	}
	case COLLISION_MRT: {
		if (tau == cachedTau) {
			for (int p = 0; p < nf; p++) {
				double sum = 0;
				const double* row = &cachedT[nf * p];
				for (int q = 0; q < nf; q++) {
					sum += row[q] * fneq[q];
				}
				delta[p] = sum;
			}
		}
		else {
			double S[27], m[27];
			diagonal(tau, S);
			for (int k = 0; k < nf; k++) {
				double sum = 0;
				const double* row = &M[nf * k];
				for (int q = 0; q < nf; q++) {
					sum += row[q] * fneq[q];
				}
				m[k] = S[k] * sum;
			}
			for (int p = 0; p < nf; p++) {
				double sum = 0;
				const double* row = &iM[nf * p];
				for (int k = 0; k < nf; k++) {
					sum += row[k] * m[k];
				}
				delta[p] = sum;
			}
		}
		break;
	}
	}
}
