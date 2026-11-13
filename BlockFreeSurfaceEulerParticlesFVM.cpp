/*
 * File:   BlockFreeSurfaceEulerParticlesFVM.cpp
 * Author: obublik
 *
 * Created on 2. února 2018, 12:27
 *
 * Free surface LBM with Eulerian particles: the number of particles in a cell (concentration) is transported
 * by a finite volume scheme (upwind / Lax-Wendroff with the van Leer limiter), the particle velocity field
 * follows the equation of motion of a single particle (drag, gravity, buoyancy).
 */

#include "BlockFreeSurfaceEulerParticlesFVM.h"
#include <stdexcept>

using namespace std;

namespace {
	const int exLW[6] = { 1, 0, 0, -1,  0,  0 };
	const int eyLW[6] = { 0, 1, 0,  0, -1,  0 };
	const int ezLW[6] = { 0, 0, 1,  0,  0, -1 };
	const int oppLW[6] = { 3, 4, 5, 0, 1, 2 };
}

BlockFreeSurfaceEulerParticlesFVM::BlockFreeSurfaceEulerParticlesFVM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurface(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {

	concentration = new double[ntot]();
	fluxConcentration = new double[ntot]();
	released = new double[ntot]();
	exchangeCorr = new double[ntot]();

	dpX = new double[ntot]();
	dpY = new double[ntot]();
	dpZ = new double[ntot]();

	UP = new double[ntot]();
	VP = new double[ntot]();
	WP = new double[ntot]();

	PP = new double[ntot](); // particle pressure
	dPPx = new double[ntot]();
	dPPy = new double[ntot]();
	dPPz = new double[ntot]();

	UF = new double[ntot]();
	VF = new double[ntot]();
	WF = new double[ntot]();

	Gx = new double[ntot]();
	Gy = new double[ntot]();
	Gz = new double[ntot]();

	GxSmooth = new double[ntot]();
	GySmooth = new double[ntot]();
	GzSmooth = new double[ntot]();

	dParticle = getParDouble("particleDiameter", par, true, 0) / Ch;
	rhoParticle = getParDouble("particleDensity", par, true, 0) / Crho;

	crossectionParticle = PI * dParticle * dParticle / 4;
	volumeParticle = PI * dParticle * dParticle * dParticle / 6;
	massParticle = rhoParticle * volumeParticle;
	Cmax = 1 / (dParticle * dParticle * dParticle);
}

BlockFreeSurfaceEulerParticlesFVM::~BlockFreeSurfaceEulerParticlesFVM() {
	delete[] concentration;
	delete[] fluxConcentration;
	delete[] released;
	delete[] exchangeCorr;
	delete[] dpX;
	delete[] dpY;
	delete[] dpZ;
	delete[] UP;
	delete[] VP;
	delete[] WP;
	delete[] PP;
	delete[] dPPx;
	delete[] dPPy;
	delete[] dPPz;
	delete[] UF;
	delete[] VF;
	delete[] WF;
	delete[] Gx;
	delete[] Gy;
	delete[] Gz;
	delete[] GxSmooth;
	delete[] GySmooth;
	delete[] GzSmooth;
}

// init/concentration<block>.txt: number of particles in each cell (order of the block files), zero if the file is missing
void BlockFreeSurfaceEulerParticlesFVM::readInitConcentration(int id, string simulPath) {
	string fileName = simulPath + "init/concentration" + fileTag + ".txt";
	ifstream myfile(fileName);
	if (myfile.is_open()) {
		string line;
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			double c = 0;
			if (!getline(myfile, line) || !(istringstream(line) >> c)) {
				throw runtime_error("Error in " + fileName + ": expected " + to_string(nAll) + " values");
			}
			concentration[ind] = c;
		}
	}

	for (int i = 0; i < ntot; i++) {
		UP[i] = U[i];
		VP[i] = V[i];
		WP[i] = W[i];
	}
}

void BlockFreeSurfaceEulerParticlesFVM::cellForce(int ind, double* acc, double& epsF) {
	// smoothed drag is a force density, the exact difference method needs an acceleration
	double rho = Rho[ind] > 0 ? Rho[ind] : 1.0;
	epsF = 1.0;
	acc[0] = fluidGravity * G[0] + GxSmooth[ind] * vof[ind] / rho;
	acc[1] = fluidGravity * G[1] + GySmooth[ind] * vof[ind] / rho;
	acc[2] = fluidGravity * G[2] + GzSmooth[ind] * vof[ind] / rho;
}

void BlockFreeSurfaceEulerParticlesFVM::collision() {
	BlockFreeSurface::collision();
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		Gx[ind] = 0;
		Gy[ind] = 0;
		Gz[ind] = 0;
		GxSmooth[ind] = 0;
		GySmooth[ind] = 0;
		GzSmooth[ind] = 0;
	}
}

// Flux of the concentration leaving cell ind through the face in the direction p (of the 6 directions).
// It is evaluated from the side of the upwind cell, so both cells of the face get the same value with
// opposite signs (conservation) and the limiter uses the upwind ratio of the gradients.
double BlockFreeSurfaceEulerParticlesFVM::faceFlux(int ind, int p) {
	int bl, loc;
	local2globalShift(indBlock, ind, bl, loc, exLW[p], eyLW[p], ezLW[p]);
	if (bl < 0 || bl >= nBlockMax || !blocks[bl] || blocks[bl]->boundary[loc]) {
		return 0;
	}
	BlockFreeSurfaceEulerParticlesFVM* bs = (BlockFreeSurfaceEulerParticlesFVM*)blocks[bl];
	if (bs->flagMap[loc] != 1 && bs->flagMap[loc] != 3) {
		return 0;
	}
	if (ghost == GHOST_NONE && bs->ghost == GHOST_FROM_FINER) {
		return 0; // the flux through a face of a finer level is computed by its cells (Mesh::forwardFluxFVM)
	}
	double Un = (bs->UP[loc] + UP[ind]) / 2 * exLW[p] + (bs->VP[loc] + VP[ind]) / 2 * eyLW[p] + (bs->WP[loc] + WP[ind]) / 2 * ezLW[p];
	if (Un < 0) {
		return -bs->upwindFlux(loc, oppLW[p], this, ind, -Un);
	}
	return upwindFlux(ind, p, bs, loc, Un);
}

// flux from cell ind (upwind) to the cell loc of the block bs in the direction p, Un >= 0
double BlockFreeSurfaceEulerParticlesFVM::upwindFlux(int ind, int p, BlockFreeSurfaceEulerParticlesFVM* bs, int loc, double Un) {
	double cA = concentration[ind];
	double cB = bs->concentration[loc];
	// low order (upwind) and high order (Lax-Wendroff) flux
	double fluxLow = Un * cA;
	double fluxHigh = Un * (cA + cB) / 2 - Un * Un * (cB - cA) / 2;

	// van Leer limiter, first order next to a wall or to the gas
	double limiter = 0;
	int blO, locO;
	int pop = oppLW[p];
	local2globalShift(indBlock, ind, blO, locO, exLW[pop], eyLW[pop], ezLW[pop]);
	if (blO > -1 && blO < nBlockMax && blocks[blO] && !blocks[blO]->boundary[locO]) {
		BlockFreeSurfaceEulerParticlesFVM* bo = (BlockFreeSurfaceEulerParticlesFVM*)blocks[blO];
		if ((bo->flagMap[locO] == 1 || bo->flagMap[locO] == 3) && cB != cA) {
			double r = (cA - bo->concentration[locO]) / (cB - cA);
			limiter = (r + fabs(r)) / (1 + fabs(r));
		}
	}
	return fluxLow - limiter * (fluxLow - fluxHigh);
}

void BlockFreeSurfaceEulerParticlesFVM::computeFlux() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		fluxConcentration[ind] = 0;
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		for (int p = 0; p < 6; p++) {
			double fl = faceFlux(ind, p);
			fluxConcentration[ind] += fl;
			if (fl != 0 && nbIndex[ind] >= 0) {
				int bl, loc;
				local2globalShift(indBlock, ind, bl, loc, exLW[p], eyLW[p], ezLW[p]);
				if (blocks[bl]->ghost == GHOST_FROM_COARSER) {
					// particles sent to (taken from) the coarser level, Mesh::forwardFluxFVM gives them to the coarser cell
					double* acc = &((BlockFreeSurfaceEulerParticlesFVM*)blocks[bl])->exchangeCorr[loc];
#pragma omp atomic
					*acc += fl;
				}
			}
		}
	}
}

void BlockFreeSurfaceEulerParticlesFVM::applyFlux() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		concentration[ind] += exchangeCorr[ind] - fluxConcentration[ind];
		exchangeCorr[ind] = 0;
	}
}

void BlockFreeSurfaceEulerParticlesFVM::prolongCell(int ind, const Block* src, int sInd) {
	BlockFreeSurface::prolongCell(ind, src, sInd);
	const BlockFreeSurfaceEulerParticlesFVM* s = (const BlockFreeSurfaceEulerParticlesFVM*)src;
	concentration[ind] = s->concentration[sInd] / 8;
	UP[ind] = s->UP[sInd];
	VP[ind] = s->VP[sInd];
	WP[ind] = s->WP[sInd];
	PP[ind] = s->PP[sInd];
}

void BlockFreeSurfaceEulerParticlesFVM::restrictCell(int ind, Block* const* src, const int* sInd) {
	BlockFreeSurface::restrictCell(ind, src, sInd);
	double c = 0, u = 0, v = 0, w3 = 0, pp = 0;
	int n = 0;
	for (int o = 0; o < 8; o++) {
		const BlockFreeSurfaceEulerParticlesFVM* s = (const BlockFreeSurfaceEulerParticlesFVM*)src[o];
		if (!s || s->boundary[sInd[o]]) {
			continue;
		}
		int si = sInd[o];
		c += s->concentration[si];
		u += s->UP[si];
		v += s->VP[si];
		w3 += s->WP[si];
		pp += s->PP[si];
		n++;
	}
	concentration[ind] = c;
	UP[ind] = n ? u / n : 0;
	VP[ind] = n ? v / n : 0;
	WP[ind] = n ? w3 / n : 0;
	PP[ind] = n ? pp / n : 0;
}

bool BlockFreeSurfaceEulerParticlesFVM::releaseDryCells() {
	anyReleased = false;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		released[ind] = 0;
		if (concentration[ind] == 0 || boundary[ind] || flagMap[ind] == 1 || flagMap[ind] == 3) {
			continue;
		}
		double s = 0;
		liquidStencil(ind, [&s](BlockFreeSurface*, int, double wgt) { s += wgt; });
		if (s > 0) { // otherwise the particles wait until the liquid returns
			released[ind] = concentration[ind] / s;
			concentration[ind] = 0;
			anyReleased = true;
		}
	}
	return anyReleased;
}

void BlockFreeSurfaceEulerParticlesFVM::collectReleased() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		double add = 0;
		stencil27(ind, [&add](BlockFreeSurface* b, int loc, double wgt) {
			add += wgt * ((BlockFreeSurfaceEulerParticlesFVM*)b)->released[loc];
		});
		concentration[ind] += add;
	}
}

void BlockFreeSurfaceEulerParticlesFVM::upwind() {
	computeFlux();
	applyFlux();
}

// Implicit step of the equation of motion of a single particle. The velocity does not depend on the
// concentration, so it is defined in all liquid cells (the fronts of the concentration need it);
// the concentration scales only the force acting on the fluid.
void BlockFreeSurfaceEulerParticlesFVM::computeNewParticleVelocityField() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			UP[ind] = 0;
			VP[ind] = 0;
			WP[ind] = 0;
			continue;
		}
		double Uslip = sqrt((UP[ind] - UF[ind]) * (UP[ind] - UF[ind]) + (VP[ind] - VF[ind]) * (VP[ind] - VF[ind]) + (WP[ind] - WF[ind]) * (WP[ind] - WF[ind]));
		double coef = 0.5 * cd(Uslip) * crossectionParticle;
		double mass = massParticle;
		double volume = volumeParticle;

		UP[ind] = (mass * UP[ind] + coef * UF[ind] + G[0] * (mass - volume)) / (mass + coef);
		VP[ind] = (mass * VP[ind] + coef * VF[ind] + G[1] * (mass - volume)) / (mass + coef);
		WP[ind] = (mass * WP[ind] + coef * WF[ind] + G[2] * (mass - volume)) / (mass + coef);

		double c = max(concentration[ind], 0.0);
		Gx[ind] = c * coef * (UP[ind] - UF[ind]);
		Gy[ind] = c * coef * (VP[ind] - VF[ind]);
		Gz[ind] = c * coef * (WP[ind] - WF[ind]);
	}
}

double BlockFreeSurfaceEulerParticlesFVM::cd(double uSlip) { // return Cd*uSlip
	double CdUslip = 0.44 * uSlip;
	double Re = dParticle * uSlip / vLB;
	if (Re < 1000) {
		CdUslip = 24.0 / dParticle * vLB * (1 + 0.15 * pow(Re, 0.687));
	}
	return CdUslip;
}

// gradient of a cell field over the liquid cells: grad = 3 * sum(w * e * value), out = scale * grad
void BlockFreeSurfaceEulerParticlesFVM::liquidGradient(const double* (*field)(Block*), double scale, double* outX, double* outY, double* outZ) {
	const double* own = field(this);
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		outX[ind] = 0;
		outY[ind] = 0;
		outZ[ind] = 0;
		if (flagMap[ind] != 1 && flagMap[ind] != 3) {
			continue;
		}
		double gx = 0, gy = 0, gz = 0;
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			double val = own[ind];
			if (bl > -1 && bl < nBlockMax && blocks[bl] && !blocks[bl]->boundary[loc]) {
				int fm = ((BlockFreeSurface*)blocks[bl])->flagMap[loc];
				if (fm == 1 || fm == 3) {
					val = field(blocks[bl])[loc];
				}
			}
			gx += w[dir] * ex[dir] * val;
			gy += w[dir] * ey[dir] * val;
			gz += w[dir] * ez[dir] * val;
		}
		outX[ind] = 3 * scale * gx;
		outY[ind] = 3 * scale * gy;
		outZ[ind] = 3 * scale * gz;
	}
}

// dp = grad(p) / 2 with p = rho / 3 (kept for the pressure gradient force 2 * dp * volume)
void BlockFreeSurfaceEulerParticlesFVM::computeGradP() {
	liquidGradient([](Block* b) -> const double* { return b->Rho; }, 1.0 / 6, dpX, dpY, dpZ);
}

// fluid velocity seen by the particles of the cell: velocity averaged over the liquid cells around
void BlockFreeSurfaceEulerParticlesFVM::interpolateFluidVelocity() {
	for (int c = 0; c < nAll; c++) {
		int ind = indAll[c];
		UF[ind] = 0;
		VF[ind] = 0;
		WF[ind] = 0;
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		double sum = 0, u = 0, v = 0, w3 = 0;
		liquidStencil(ind, [&](BlockFreeSurface* b, int loc, double wgt) {
			u += wgt * b->U[loc];
			v += wgt * b->V[loc];
			w3 += wgt * b->W[loc];
			sum += wgt;
		});
		UF[ind] = u / sum;
		VF[ind] = v / sum;
		WF[ind] = w3 / sum;
	}
}

// divides the force of the cell by the part of the kernel that reaches the liquid (see smoothForce)
void BlockFreeSurfaceEulerParticlesFVM::normalizeForce() {
	for (int c = 0; c < nAll; c++) {
		int ind = indAll[c];
		if (Gx[ind] == 0 && Gy[ind] == 0 && Gz[ind] == 0) {
			continue;
		}
		double s = max(liquidKernelSum(ind), FORCE_KERNEL_MIN);
		Gx[ind] /= s;
		Gy[ind] /= s;
		Gz[ind] /= s;
	}
}

// Spreads the force of the particles to the liquid cells around. Every cell gathers the normalized forces
// of its neighbours (it writes only its own cells, so the blocks can run in parallel) and the sum of
// GSmooth * vof over the cells equals the sum of the forces G.
void BlockFreeSurfaceEulerParticlesFVM::smoothForce() {
	for (int c = 0; c < nAll; c++) {
		int ind = indAll[c];
		GxSmooth[ind] = 0;
		GySmooth[ind] = 0;
		GzSmooth[ind] = 0;
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		double gx = 0, gy = 0, gz = 0;
		liquidStencil(ind, [&](BlockFreeSurface* b, int loc, double wgt) {
			BlockFreeSurfaceEulerParticlesFVM* e = (BlockFreeSurfaceEulerParticlesFVM*)b;
			gx += wgt * e->Gx[loc];
			gy += wgt * e->Gy[loc];
			gz += wgt * e->Gz[loc];
		});
		GxSmooth[ind] = gx;
		GySmooth[ind] = gy;
		GzSmooth[ind] = gz;
	}
}

void BlockFreeSurfaceEulerParticlesFVM::computeParticlePressure() {
	// the particle pressure (packing limit) is not modelled
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		PP[ind] = 0;
	}
}

void BlockFreeSurfaceEulerParticlesFVM::computeParticlePressureGradient() {
	liquidGradient([](Block* b) -> const double* { return ((BlockFreeSurfaceEulerParticlesFVM*)b)->PP; }, 1.0 / 6, dPPx, dPPy, dPPz);
}
