/*
 * File:   BlockFreeSurfaceEulerNoSphereParticlesLBM.cpp
 *
 * Eulerian non-spherical particles transported by LBM, problemID 42.
 * See BlockFreeSurfaceEulerNoSphereParticlesLBM.h.
 */

#include "BlockFreeSurfaceEulerNoSphereParticlesLBM.h"
#include <cmath>

using namespace std;

BlockFreeSurfaceEulerNoSphereParticlesLBM::BlockFreeSurfaceEulerNoSphereParticlesLBM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurfaceEulerParticlesLBM(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {
	initOrientation(par, dParticle, ntot, "problemID 42");
	orientNew = new double[(size_t)6 * ntot]();
	exchangeOrient = new double[(size_t)6 * ntot]();
	releasedOrient = new double[(size_t)6 * ntot]();
}

BlockFreeSurfaceEulerNoSphereParticlesLBM::~BlockFreeSurfaceEulerNoSphereParticlesLBM() {
	delete[] orientNew;
	delete[] exchangeOrient;
	delete[] releasedOrient;
}

// the concentration as problemID 4, then init/orientation<tag>.txt (A11 A22 A33 A12 A13 A23 per cell), A = I/3 without the file
void BlockFreeSurfaceEulerNoSphereParticlesLBM::readInitConcentration(int id, string simulPath) {
	BlockFreeSurfaceEulerParticlesLBM::readInitConcentration(id, simulPath);
	readOrientation(simulPath + "init/orientation" + fileTag + ".txt", nAll, indAll);
}

void BlockFreeSurfaceEulerNoSphereParticlesLBM::prepareParticleStep() {
	BlockFreeSurfaceEulerParticlesLBM::prepareParticleStep();
	velocityGradient(gradU);
}

// explicit convection and implicit drag with the drag tensor of the cell (problemID 4 with a tensor)
void BlockFreeSurfaceEulerNoSphereParticlesLBM::computeNewParticleVelocityField(double* dist, double* gradDist) {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] != 1 && flagMap[ind] != 3) {
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
		const double conv[3] = {
			v[0] * gUx[ind] + v[1] * gUy[ind] + v[2] * gUz[ind],
			v[0] * gVx[ind] + v[1] * gVy[ind] + v[2] * gVz[ind],
			v[0] * gWx[ind] + v[1] * gWy[ind] + v[2] * gWz[ind] };
		double sys[9], inv[9], rhs[3];
		for (int k = 0; k < 9; k++) {
			sys[k] = K[k];
		}
		sys[0] += mass;
		sys[4] += mass;
		sys[8] += mass;
		for (int a = 0; a < 3; a++) {
			rhs[a] = mass * (v[a] - conv[a]) + K[3 * a] * u[0] + K[3 * a + 1] * u[1] + K[3 * a + 2] * u[2] + G[a] * (mass - volumeParticle);
		}
		inverse3(sys, inv);
		for (int a = 0; a < 3; a++) {
			v[a] = inv[3 * a] * rhs[0] + inv[3 * a + 1] * rhs[1] + inv[3 * a + 2] * rhs[2];
		}
		UP[ind] = v[0];
		VP[ind] = v[1];
		WP[ind] = v[2];
		wallSlipLimit(ind, dist, gradDist);

		double c = max(concentration[ind], 0.0);
		double w3[3] = { UP[ind] - u[0], VP[ind] - u[1], WP[ind] - u[2] };
		Gx[ind] = c * (K[0] * w3[0] + K[1] * w3[1] + K[2] * w3[2]);
		Gy[ind] = c * (K[3] * w3[0] + K[4] * w3[1] + K[5] * w3[2]);
		Gz[ind] = c * (K[6] * w3[0] + K[7] * w3[1] + K[8] * w3[2]);
	}
}

// Streaming of problemID 4 and the transport of the orientation by the net fluxes of the links (as the upwind A of
// problemID 33): through the link to the neighbour y the cell sends out = gc_p(x) and receives in = gc_opp(y), the net
// flux carries A of the upwind cell. With the sum of the populations S after the streaming
//   S A_new = S A + sum_links max(in - out, 0) (A_y - A)       (orientNew holds S A_new)
// which conserves sum(S A) exactly and is a convex combination while the cell does not send out more than it holds.
// Assigning every population the A of its source cell instead would mix A of the neighbours in every step (the
// gross exchange w_i c is much larger than the net flux), A would diffuse ~100 times faster than the particles.
// Bounce-back (wall, gas): in = out, no transport. If the weight of the cell's own A is negative (negative
// populations at sharp fronts) the convex weights are renormalised (the measure is then not exactly conserved);
// with S <= 0 the cell keeps its A.
void BlockFreeSurfaceEulerNoSphereParticlesLBM::propagateScalarCell(int ind) {
	if (flagMap[ind] != 1 && flagMap[ind] != 3) {
		return;
	}
	double dq[6] = { 0, 0, 0, 0, 0, 0 };
	double S = 0, M = 0;
	const double* own = &orient[6 * ind];
	for (int p = 0; p < nf; p++) {
		int bl = nbBlock(ind, p);
		int loc = nbLocal(ind, p);
		bool open = p > 0 && bl > -1 && bl < nBlockMax && blocks[bl] && !blocks[bl]->boundary[loc];
		if (open) {
			int fm = ((BlockFreeSurface*)blocks[bl])->flagMap[loc];
			open = fm == 1 || fm == 3;
		}
		double gp;
		if (open) {
			BlockFreeSurfaceEulerNoSphereParticlesLBM* nb = (BlockFreeSurfaceEulerNoSphereParticlesLBM*)blocks[bl];
			gp = nb->gc[nf * loc + opp[p]];
			double m = gp - gc[nf * ind + p]; // net flux into the cell
			if (m > 0) {
				const double* A = &nb->orient[6 * loc];
				M += m;
				for (int k = 0; k < 6; k++) {
					dq[k] += m * (A[k] - own[k]);
				}
			}
		}
		else {
			gp = gc[nf * ind + p];
		}
		g[nf * ind + opp[p]] = gp;
		S += gp;
	}
	double* out = &orientNew[6 * ind];
	if (S <= 0) {
		for (int k = 0; k < 6; k++) {
			out[k] = S * own[k];
		}
		return;
	}
	double f = S - M >= 0 ? 1.0 : S / M; // weight of the own A is S - f M >= 0
	for (int k = 0; k < 6; k++) {
		out[k] = S * own[k] + f * dq[k];
	}
}

void BlockFreeSurfaceEulerNoSphereParticlesLBM::applyExchangeCorr() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		double* eq = &exchangeOrient[6 * ind];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			double q[6];
			bool exchange = !boundary[ind];
			for (int k = 0; k < 6; k++) {
				q[k] = orientNew[6 * ind + k] + (exchange ? eq[k] : 0.0);
			}
			// trace q = sum of the populations (with the exchange); without particles the cell keeps its A
			if (q[0] + q[1] + q[2] > 0) {
				double* A = &orient[6 * ind];
				for (int k = 0; k < 6; k++) {
					A[k] = q[k];
				}
				projectOrientation(A);
			}
		}
		for (int k = 0; k < 6; k++) {
			eq[k] = 0;
		}
	}
	BlockFreeSurfaceEulerParticlesLBM::applyExchangeCorr();

	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		orientationStep(&orient[6 * ind], &gradU[9 * ind], 1.0);
	}
}

bool BlockFreeSurfaceEulerNoSphereParticlesLBM::releaseDryCells() {
	bool any = BlockFreeSurfaceEulerParticlesLBM::releaseDryCells();
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		for (int k = 0; k < 6; k++) {
			releasedOrient[6 * ind + k] = released[ind] * orient[6 * ind + k];
		}
	}
	return any;
}

// A of a liquid cell collecting the released particles: mean weighted by the sums of the populations
void BlockFreeSurfaceEulerNoSphereParticlesLBM::collectReleased() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		double add = 0, addQ[6] = { 0, 0, 0, 0, 0, 0 };
		stencil27(ind, [&](BlockFreeSurface* b, int loc, double wgt) {
			BlockFreeSurfaceEulerNoSphereParticlesLBM* e = (BlockFreeSurfaceEulerNoSphereParticlesLBM*)b;
			add += wgt * e->released[loc];
			for (int k = 0; k < 6; k++) {
				addQ[k] += wgt * e->releasedOrient[6 * loc + k];
			}
		});
		if (add == 0) {
			continue;
		}
		double S = 0;
		for (int p = 0; p < nf; p++) {
			S += gc[nf * ind + p];
		}
		double* A = &orient[6 * ind];
		double q[6];
		for (int k = 0; k < 6; k++) {
			q[k] = S * A[k] + addQ[k];
		}
		if (q[0] + q[1] + q[2] > 0) {
			for (int k = 0; k < 6; k++) {
				A[k] = q[k];
			}
			projectOrientation(A);
		}
	}
	BlockFreeSurfaceEulerParticlesLBM::collectReleased();
}

void BlockFreeSurfaceEulerNoSphereParticlesLBM::prolongCell(int ind, const Block* src, int sInd) {
	BlockFreeSurfaceEulerParticlesLBM::prolongCell(ind, src, sInd);
	const BlockFreeSurfaceEulerNoSphereParticlesLBM* s = (const BlockFreeSurfaceEulerNoSphereParticlesLBM*)src;
	for (int k = 0; k < 6; k++) {
		orient[6 * ind + k] = s->orient[6 * sInd + k];
	}
}

// A of the coarse cell: mean of the 8 fine cells weighted by their concentration (plain mean without particles)
void BlockFreeSurfaceEulerNoSphereParticlesLBM::restrictCell(int ind, Block* const* src, const int* sInd) {
	BlockFreeSurfaceEulerParticlesLBM::restrictCell(ind, src, sInd);
	double q[6] = { 0, 0, 0, 0, 0, 0 }, m[6] = { 0, 0, 0, 0, 0, 0 };
	double c = 0;
	int n = 0;
	for (int o = 0; o < 8; o++) {
		const BlockFreeSurfaceEulerNoSphereParticlesLBM* s = (const BlockFreeSurfaceEulerNoSphereParticlesLBM*)src[o];
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
