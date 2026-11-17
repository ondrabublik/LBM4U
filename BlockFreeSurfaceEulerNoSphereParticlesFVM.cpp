/*
 * File:   BlockFreeSurfaceEulerNoSphereParticlesFVM.cpp
 *
 * Eulerian non-spherical particles with the orientation tensor, problemID 33.
 * See BlockFreeSurfaceEulerNoSphereParticlesFVM.h.
 */

#include "BlockFreeSurfaceEulerNoSphereParticlesFVM.h"
#include <cmath>
#include <stdexcept>

using namespace std;

BlockFreeSurfaceEulerNoSphereParticlesFVM::BlockFreeSurfaceEulerNoSphereParticlesFVM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurfaceEulerParticlesFVM(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {
	extraTransport = true;
	initOrientation(par, dParticle, ntot, "problemID 33");
	fluxOrient = new double[(size_t)6 * ntot]();
	exchangeOrient = new double[(size_t)6 * ntot]();
	releasedOrient = new double[(size_t)6 * ntot]();
}

BlockFreeSurfaceEulerNoSphereParticlesFVM::~BlockFreeSurfaceEulerNoSphereParticlesFVM() {
	delete[] fluxOrient;
	delete[] exchangeOrient;
	delete[] releasedOrient;
}

// the concentration as problemID 3, then init/orientation<tag>.txt (A11 A22 A33 A12 A13 A23 per cell), A = I/3 without the file
void BlockFreeSurfaceEulerNoSphereParticlesFVM::readInitConcentration(int id, string simulPath) {
	BlockFreeSurfaceEulerParticlesFVM::readInitConcentration(id, simulPath);
	readOrientation(simulPath + "init/orientation" + fileTag + ".txt", nAll, indAll);
}

void BlockFreeSurfaceEulerNoSphereParticlesFVM::prepareParticleStep() {
	BlockFreeSurfaceEulerParticlesFVM::prepareParticleStep();
	velocityGradient(gradU);
}

// implicit step of the particle velocity with the drag tensor of the cell (problemID 3 with a tensor)
void BlockFreeSurfaceEulerNoSphereParticlesFVM::computeNewParticleVelocityField() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			UP[ind] = 0;
			VP[ind] = 0;
			WP[ind] = 0;
			continue;
		}
		double v[3] = { UP[ind], VP[ind], WP[ind] };
		double u[3] = { UF[ind], VF[ind], WF[ind] };
		double slip = sqrt((v[0] - u[0]) * (v[0] - u[0]) + (v[1] - u[1]) * (v[1] - u[1]) + (v[2] - u[2]) * (v[2] - u[2]));
		double K[9];
		dragTensor(&orient[6 * ind], K);
		double s = vLB * dragCorrection(dParticle * slip / vLB);
		for (int k = 0; k < 9; k++) {
			K[k] *= s;
		}
		double mass = massParticle;
		double sys[9], inv[9], rhs[3];
		for (int k = 0; k < 9; k++) {
			sys[k] = K[k];
		}
		sys[0] += mass;
		sys[4] += mass;
		sys[8] += mass;
		for (int a = 0; a < 3; a++) {
			rhs[a] = mass * v[a] + K[3 * a] * u[0] + K[3 * a + 1] * u[1] + K[3 * a + 2] * u[2] + G[a] * (mass - volumeParticle);
		}
		inverse3(sys, inv);
		for (int a = 0; a < 3; a++) {
			v[a] = inv[3 * a] * rhs[0] + inv[3 * a + 1] * rhs[1] + inv[3 * a + 2] * rhs[2];
		}
		UP[ind] = v[0];
		VP[ind] = v[1];
		WP[ind] = v[2];

		double c = max(concentration[ind], 0.0);
		double w3[3] = { v[0] - u[0], v[1] - u[1], v[2] - u[2] };
		Gx[ind] = c * (K[0] * w3[0] + K[1] * w3[1] + K[2] * w3[2]);
		Gy[ind] = c * (K[3] * w3[0] + K[4] * w3[1] + K[5] * w3[2]);
		Gz[ind] = c * (K[6] * w3[0] + K[7] * w3[1] + K[8] * w3[2]);
	}
}

// the flux of c A is the flux of c times A of the upwind cell
void BlockFreeSurfaceEulerNoSphereParticlesFVM::faceTransport(int ind, double fl, BlockFreeSurfaceEulerParticlesFVM* nb, int loc, bool toCoarser) {
	BlockFreeSurfaceEulerNoSphereParticlesFVM* n = (BlockFreeSurfaceEulerNoSphereParticlesFVM*)nb;
	const double* Aup = fl > 0 ? &orient[6 * ind] : &n->orient[6 * loc];
	for (int k = 0; k < 6; k++) {
		fluxOrient[6 * ind + k] += fl * Aup[k];
	}
	if (toCoarser) {
		for (int k = 0; k < 6; k++) {
			double* acc = &n->exchangeOrient[6 * loc + k];
#pragma omp atomic
			*acc += fl * Aup[k];
		}
	}
}

void BlockFreeSurfaceEulerNoSphereParticlesFVM::forwardExchange(int ind, BlockFreeSurfaceEulerParticlesFVM* dst, int dloc) {
	BlockFreeSurfaceEulerNoSphereParticlesFVM* d = (BlockFreeSurfaceEulerNoSphereParticlesFVM*)dst;
	for (int k = 0; k < 6; k++) {
		double q = exchangeOrient[6 * ind + k];
		if (q != 0) {
			exchangeOrient[6 * ind + k] = 0;
			double* acc = &d->exchangeOrient[6 * dloc + k];
#pragma omp atomic
			*acc += q;
		}
	}
	BlockFreeSurfaceEulerParticlesFVM::forwardExchange(ind, dst, dloc);
}

// c A after the fluxes, A = (c A) / trace(c A); then the Folgar-Tucker step in the liquid cells
void BlockFreeSurfaceEulerNoSphereParticlesFVM::applyFlux() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		double* A = &orient[6 * ind];
		double* fq = &fluxOrient[6 * ind];
		double* eq = &exchangeOrient[6 * ind];
		if (fluxConcentration[ind] != 0 || exchangeCorr[ind] != 0) {
			double c = concentration[ind];
			double q[6];
			for (int k = 0; k < 6; k++) {
				q[k] = c * A[k] + eq[k] - fq[k];
			}
			if (q[0] + q[1] + q[2] > 0) {
				for (int k = 0; k < 6; k++) {
					A[k] = q[k];
				}
				projectOrientation(A);
			}
		}
		for (int k = 0; k < 6; k++) {
			fq[k] = 0;
			eq[k] = 0;
		}
	}
	BlockFreeSurfaceEulerParticlesFVM::applyFlux();

	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		orientationStep(&orient[6 * ind], &gradU[9 * ind], 1.0);
	}
}

bool BlockFreeSurfaceEulerNoSphereParticlesFVM::releaseDryCells() {
	bool any = BlockFreeSurfaceEulerParticlesFVM::releaseDryCells();
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		for (int k = 0; k < 6; k++) {
			releasedOrient[6 * ind + k] = released[ind] * orient[6 * ind + k];
		}
	}
	return any;
}

void BlockFreeSurfaceEulerNoSphereParticlesFVM::collectReleased() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		double add = 0, addQ[6] = { 0, 0, 0, 0, 0, 0 };
		stencil27(ind, [&](BlockFreeSurface* b, int loc, double wgt) {
			BlockFreeSurfaceEulerNoSphereParticlesFVM* e = (BlockFreeSurfaceEulerNoSphereParticlesFVM*)b;
			add += wgt * e->released[loc];
			for (int k = 0; k < 6; k++) {
				addQ[k] += wgt * e->releasedOrient[6 * loc + k];
			}
		});
		if (add == 0) {
			continue;
		}
		double* A = &orient[6 * ind];
		double c = concentration[ind];
		for (int k = 0; k < 6; k++) {
			A[k] = c * A[k] + addQ[k];
		}
		projectOrientation(A);
		concentration[ind] += add;
	}
}

void BlockFreeSurfaceEulerNoSphereParticlesFVM::prolongCell(int ind, const Block* src, int sInd) {
	BlockFreeSurfaceEulerParticlesFVM::prolongCell(ind, src, sInd);
	const BlockFreeSurfaceEulerNoSphereParticlesFVM* s = (const BlockFreeSurfaceEulerNoSphereParticlesFVM*)src;
	for (int k = 0; k < 6; k++) {
		orient[6 * ind + k] = s->orient[6 * sInd + k];
	}
}

// A of the coarse cell: mean of the 8 fine cells weighted by their concentration (plain mean without particles)
void BlockFreeSurfaceEulerNoSphereParticlesFVM::restrictCell(int ind, Block* const* src, const int* sInd) {
	BlockFreeSurfaceEulerParticlesFVM::restrictCell(ind, src, sInd);
	double q[6] = { 0, 0, 0, 0, 0, 0 }, m[6] = { 0, 0, 0, 0, 0, 0 };
	double c = 0;
	int n = 0;
	for (int o = 0; o < 8; o++) {
		const BlockFreeSurfaceEulerNoSphereParticlesFVM* s = (const BlockFreeSurfaceEulerNoSphereParticlesFVM*)src[o];
		if (!s || s->boundary[sInd[o]]) {
			continue;
		}
		const double* A = &s->orient[6 * sInd[o]];
		double cs = max(s->concentration[sInd[o]], 0.0);
		for (int k = 0; k < 6; k++) {
			q[k] += cs * A[k];
			m[k] += A[k];
		}
		c += cs;
		n++;
	}
	if (n == 0) {
		return;
	}
	double* A = &orient[6 * ind];
	for (int k = 0; k < 6; k++) {
		A[k] = c > 0 ? q[k] : m[k];
	}
	projectOrientation(A);
}
