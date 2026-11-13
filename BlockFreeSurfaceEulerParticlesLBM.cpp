/*
 * File:   BlockFreeSurfaceEulerParticlesLBM.cpp
 * Author: obublik
 *
 * Created on 2. února 2018, 12:27
 *
 * Free surface LBM with Eulerian particles: the number of particles in a cell (concentration) is transported
 * by a second set of populations (advection-diffusion LBM), the particle velocity field follows the equation
 * of motion of a single particle (convection, drag, gravity, buoyancy).
 */

#include "BlockFreeSurfaceEulerParticlesLBM.h"
#include <stdexcept>

using namespace std;

BlockFreeSurfaceEulerParticlesLBM::BlockFreeSurfaceEulerParticlesLBM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurface(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {

	concentration = new double[ntot]();
	gc = new double[ntot * nf]();
	g = new double[ntot * nf]();
	released = new double[ntot]();

	dpX = new double[ntot]();
	dpY = new double[ntot]();
	dpZ = new double[ntot]();

	dcX = new double[ntot]();
	dcY = new double[ntot]();
	dcZ = new double[ntot]();

	gVOFpX = new double[ntot]();
	gVOFpY = new double[ntot]();
	gVOFpZ = new double[ntot]();

	UP = new double[ntot]();
	VP = new double[ntot]();
	WP = new double[ntot]();

	gUx = new double[ntot]();
	gUy = new double[ntot]();
	gUz = new double[ntot]();

	gVx = new double[ntot]();
	gVy = new double[ntot]();
	gVz = new double[ntot]();

	gWx = new double[ntot]();
	gWy = new double[ntot]();
	gWz = new double[ntot]();

	UF = new double[ntot]();
	VF = new double[ntot]();
	WF = new double[ntot]();

	Gx = new double[ntot]();
	Gy = new double[ntot]();
	Gz = new double[ntot]();

	GxSmooth = new double[ntot]();
	GySmooth = new double[ntot]();
	GzSmooth = new double[ntot]();

	tauDir = new double[ntot * nf]();
	user_def = new double[ntot]();

	dParticle = getParDouble("particleDiameter", par, true, 0) / Ch;
	rhoParticle = getParDouble("particleDensity", par, true, 0) / Crho;
	diffusLB = getParDouble("particleDiffusivity", par, true, 0) / Cvis;
	interPartFCoef = getParDouble("particleInterForces", par, false, 0) / Cf;
	// The equilibrium has the zeroth moment theta * c and the second moment c / 3, so the lattice
	// "speed of sound" of the scalar is 1 / (3 theta) and the diffusivity is (taus - 1/2) / (3 theta).
	// The rest population c * (theta - 1 + w0) has to stay positive (theta = 1 / 1.4 is too small for D3Q15).
	theta = max(1.0 / 1.4, 1.0 - w[0] + (1.0 / 1.4 - 2.0 / 3.0));
	taus = 1.0 / 2 + 3 * theta * diffusLB;
	double diffusLBoverPack = getParDouble("overpackDiffusivity", par, false, 0);
	tauOverPack = 3 * diffusLBoverPack;

	positiveDirs = new int[(nf - 1) / 2];
	vector<int> aux(nf, 0);
	int s = 0;
	for (int p = 1; p < nf; p++) {
		if (aux[p] != 1) {
			positiveDirs[s] = p;
			aux[p] = 1;
			aux[opp[p]] = 1;
			s++;
		}
	}

	for (int i = 0; i < ntot * nf; i++) {
		tauDir[i] = taus;
	}

	crossectionParticle = PI * dParticle * dParticle / 4;
	volumeParticle = PI * dParticle * dParticle * dParticle / 6;
	massParticle = rhoParticle * volumeParticle;
	concentrationLimit = 0.74 / volumeParticle; // close packing of spheres
	concentrationSum = 0;
	concentrationMax = 0;

	if (isFirst) {
		cout << "taus = " << taus << endl;
		cout << "maximal allowed value of concetration = " << concentrationLimit << endl;
		if (taus < 0.505) {
			cout << "warning: taus is close to 0.5, the transport of the particles may oscillate (increase particleDiffusivity)" << endl;
		}
	}
}

BlockFreeSurfaceEulerParticlesLBM::~BlockFreeSurfaceEulerParticlesLBM() {
	delete[] gc;
	delete[] g;
	delete[] released;
	delete[] dpX;
	delete[] dpY;
	delete[] dpZ;
	delete[] dcX;
	delete[] dcY;
	delete[] dcZ;
	delete[] gVOFpX;
	delete[] gVOFpY;
	delete[] gVOFpZ;
	delete[] UP;
	delete[] VP;
	delete[] WP;
	delete[] UF;
	delete[] VF;
	delete[] WF;
	delete[] Gx;
	delete[] Gy;
	delete[] Gz;
	delete[] GxSmooth;
	delete[] GySmooth;
	delete[] GzSmooth;
	delete[] concentration;
	delete[] gUx;
	delete[] gUy;
	delete[] gUz;
	delete[] gVx;
	delete[] gVy;
	delete[] gVz;
	delete[] gWx;
	delete[] gWy;
	delete[] gWz;
	delete[] tauDir;
	delete[] user_def;
	delete[] positiveDirs;
}

// init/concentration<block>.txt: number of particles in each cell (order of the block files), zero if the file is missing
void BlockFreeSurfaceEulerParticlesLBM::readInitConcentration(int id, string simulPath) {
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

void BlockFreeSurfaceEulerParticlesLBM::equilibriumScalar(double* geq, double c, double Up, double Vp, double Wp) {
	double sumg = 0.0;
	for (int p = 1; p < nf; p++) {
		double x = ex[p] * Up + ey[p] * Vp + ez[p] * Wp;
		geq[p] = w[p] * c * (1.0 + 3.0 * x * theta);
		sumg += geq[p];
	}
	geq[0] = c * theta - sumg;
}

void BlockFreeSurfaceEulerParticlesLBM::computeInitDistributionScalar() {
	double geq[27];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			equilibriumScalar(geq, concentration[ind], UP[ind], VP[ind], WP[ind]);
			for (int p = 0; p < nf; p++) {
				gc[nf * ind + p] = geq[p];
				g[nf * ind + p] = geq[p];
			}
		}
	}
}

// streaming of the scalar populations, walls and the gas reflect them (no flux of particles)
void BlockFreeSurfaceEulerParticlesLBM::propagationScalar() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] != 1 && flagMap[ind] != 3) {
			continue;
		}
		for (int p = 0; p < nf; p++) {
			int bl = nbBlock(ind, p);
			int loc = nbLocal(ind, p);
			bool open = bl > -1 && bl < nBlockMax && blocks[bl] && !blocks[bl]->boundary[loc];
			if (open) {
				int fm = ((BlockFreeSurface*)blocks[bl])->flagMap[loc];
				open = fm == 1 || fm == 3;
			}
			if (open) {
				g[nf * ind + opp[p]] = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[bl])->gc[nf * loc + opp[p]];
			}
			else {
				g[nf * ind + opp[p]] = gc[nf * ind + p];
			}
		}
	}
}

bool BlockFreeSurfaceEulerParticlesLBM::releaseDryCells() {
	anyReleased = false;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		released[ind] = 0;
		if (boundary[ind] || flagMap[ind] == 1 || flagMap[ind] == 3) {
			continue;
		}
		double sum = 0;
		for (int p = 0; p < nf; p++) {
			sum += gc[nf * ind + p];
		}
		if (sum == 0) {
			continue;
		}
		double s = 0;
		liquidStencil(ind, [&s](BlockFreeSurface*, int, double wgt) { s += wgt; });
		if (s > 0) { // otherwise the particles wait until the liquid returns
			released[ind] = sum / s;
			for (int p = 0; p < nf; p++) {
				gc[nf * ind + p] = 0;
				g[nf * ind + p] = 0;
			}
			anyReleased = true;
		}
	}
	return anyReleased;
}

// the collected particles are added as the equilibrium at rest
void BlockFreeSurfaceEulerParticlesLBM::collectReleased() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		double add = 0;
		stencil27(ind, [&add](BlockFreeSurface* b, int loc, double wgt) {
			add += wgt * ((BlockFreeSurfaceEulerParticlesLBM*)b)->released[loc];
		});
		if (add == 0) {
			continue;
		}
		double c = add / theta;
		concentration[ind] += c;
		gc[nf * ind] += c * (theta - 1.0 + w[0]);
		for (int p = 1; p < nf; p++) {
			gc[nf * ind + p] += c * w[p];
		}
	}
}

// The concentration is defined so that the collision with the relaxation times of the directions
// conserves the sum of the populations: sum((g - geq(c)) / tau) = 0.
void BlockFreeSurfaceEulerParticlesLBM::calculConcentration() {
	concentrationSum = 0.0;
	concentrationMax = 0.0;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		concentration[ind] = 0.0;
		user_def[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			const double* td = &tauDir[nf * ind];
			// sum(geq / tau) / c, the velocity terms cancel because opposite directions share tau
			double auxsum = (theta - 1.0 + w[0]) / td[0];
			double sum = g[nf * ind] / td[0];
			for (int p = 1; p < nf; p++) {
				auxsum += w[p] / td[p];
				sum += g[nf * ind + p] / td[p];
				user_def[ind] = max(user_def[ind], td[p]);
			}
			concentration[ind] = sum / auxsum;
			concentrationSum += concentration[ind];
			concentrationMax = max(concentrationMax, concentration[ind]);
		}
	}
}

void BlockFreeSurfaceEulerParticlesLBM::collisionScalar() {
	double geq[27];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			equilibriumScalar(geq, concentration[ind], UP[ind], VP[ind], WP[ind]);
			for (int p = 0; p < nf; p++) {
				int pind = nf * ind + p;
				gc[pind] = g[pind] - 1.0 / tauDir[nf * ind + p] * (g[pind] - geq[p]);
			}
		}
	}
}

// Explicit convection and implicit drag of the equation of motion of a single particle. The velocity
// does not depend on the concentration, so it is defined in all liquid cells; the concentration scales
// only the force acting on the fluid.
void BlockFreeSurfaceEulerParticlesLBM::computeNewParticleVelocityField(double* dist, double* gradDist) {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] != 1 && flagMap[ind] != 3) {
			UP[ind] = 0;
			VP[ind] = 0;
			WP[ind] = 0;
			continue;
		}
		double mass = massParticle;
		double volume = volumeParticle;
		double Uslip = sqrt((UP[ind] - UF[ind]) * (UP[ind] - UF[ind]) + (VP[ind] - VF[ind]) * (VP[ind] - VF[ind]) + (WP[ind] - WF[ind]) * (WP[ind] - WF[ind]));
		double coef = 0.5 * cd(Uslip) * crossectionParticle;

		double up = UP[ind], vp = VP[ind], wp = WP[ind];
		UP[ind] = (mass * (up - up * gUx[ind] - vp * gUy[ind] - wp * gUz[ind]) + coef * UF[ind] + G[0] * (mass - volume)) / (mass + coef);
		VP[ind] = (mass * (vp - up * gVx[ind] - vp * gVy[ind] - wp * gVz[ind]) + coef * VF[ind] + G[1] * (mass - volume)) / (mass + coef);
		WP[ind] = (mass * (wp - up * gWx[ind] - vp * gWy[ind] - wp * gWz[ind]) + coef * WF[ind] + G[2] * (mass - volume)) / (mass + coef);

		if (dist[ind] < 1) {
			// next to a wall the particles do not move towards it faster than the fluid, along the wall they are free
			const double* n = &gradDist[3 * ind];
			double nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
			double slip = (UP[ind] - UF[ind]) * n[0] + (VP[ind] - VF[ind]) * n[1] + (WP[ind] - WF[ind]) * n[2];
			if (nn > 0 && slip < 0) {
				UP[ind] -= slip * n[0] / nn;
				VP[ind] -= slip * n[1] / nn;
				WP[ind] -= slip * n[2] / nn;
			}
		}

		double c = max(concentration[ind], 0.0);
		Gx[ind] = c * coef * (UP[ind] - UF[ind]);
		Gy[ind] = c * coef * (VP[ind] - VF[ind]);
		Gz[ind] = c * coef * (WP[ind] - WF[ind]);
	}
}

double BlockFreeSurfaceEulerParticlesLBM::cd(double uSlip) { // return Cd*uSlip
	double CdUslip = 0.44 * uSlip;
	double Re = dParticle * uSlip / vLB;
	if (Re < 1000) {
		CdUslip = 24.0 / dParticle * vLB * (1 + 0.15 * pow(Re, 0.687));
	}
	return CdUslip;
}

// gradient of a cell field: grad = 3 * sum(w * e * value), out = scale * grad. Missing neighbours (walls,
// with liquidOnly also the gas) take the value of the cell; with liquidOnly the gradient is zero outside of the liquid.
void BlockFreeSurfaceEulerParticlesLBM::fieldGradient(const double* (*field)(Block*), double scale, bool liquidOnly, double* outX, double* outY, double* outZ) {
	const double* own = field(this);
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		outX[ind] = 0;
		outY[ind] = 0;
		outZ[ind] = 0;
		if (liquidOnly && flagMap[ind] != 1 && flagMap[ind] != 3) {
			continue;
		}
		double gx = 0, gy = 0, gz = 0;
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			double val = own[ind];
			if (bl > -1 && bl < nBlockMax && blocks[bl] && !blocks[bl]->boundary[loc]) {
				int fm = ((BlockFreeSurface*)blocks[bl])->flagMap[loc];
				if (!liquidOnly || fm == 1 || fm == 3) {
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

#define LBM_FIELD(member) [](Block* b) -> const double* { return ((BlockFreeSurfaceEulerParticlesLBM*)b)->member; }

// dp = grad(p) / 2 with p = rho / 3 (kept for the pressure gradient force 2 * dp * volume)
void BlockFreeSurfaceEulerParticlesLBM::computeGradP() {
	fieldGradient(LBM_FIELD(Rho), 1.0 / 6, true, dpX, dpY, dpZ);
}

// gradient of the particle velocity for the convective term
void BlockFreeSurfaceEulerParticlesLBM::computeGradU() {
	fieldGradient(LBM_FIELD(UP), 1.0, true, gUx, gUy, gUz);
	fieldGradient(LBM_FIELD(VP), 1.0, true, gVx, gVy, gVz);
	fieldGradient(LBM_FIELD(WP), 1.0, true, gWx, gWy, gWz);
}

// dc = grad(c) / 6
void BlockFreeSurfaceEulerParticlesLBM::computeGradC() {
	fieldGradient(LBM_FIELD(concentration), 1.0 / 6, true, dcX, dcY, dcZ);
}

// sum(w * e * c) / (cLimit * sum(w)): gradient of the volume fraction of the particles relative to the packing limit
void BlockFreeSurfaceEulerParticlesLBM::computeGradVOFp() {
	fieldGradient(LBM_FIELD(concentration), 1.0 / (3 * concentrationLimit * (1.0 - w[0])), false, gVOFpX, gVOFpY, gVOFpZ);
}

void BlockFreeSurfaceEulerParticlesLBM::cellForce(int ind, double* acc, double& epsF) {
	// smoothed drag is a force density, the exact difference method needs an acceleration
	double rho = Rho[ind] > 0 ? Rho[ind] : 1.0;
	epsF = limitFraction(1 - concentration[ind] * volumeParticle);
	acc[0] = G[0] * epsF + GxSmooth[ind] * vof[ind] / rho;
	acc[1] = G[1] * epsF + GySmooth[ind] * vof[ind] / rho;
	acc[2] = G[2] * epsF + GzSmooth[ind] * vof[ind] / rho;
}

void BlockFreeSurfaceEulerParticlesLBM::collision() {
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

double BlockFreeSurfaceEulerParticlesLBM::limitFraction(double a) {
	if (a > 1) {
		a = 1;
	}
	if (a < 0.3) {
		a = 0.3;
	}

	return a;
}

// fluid velocity seen by the particles of the cell: velocity averaged over the liquid cells around
void BlockFreeSurfaceEulerParticlesLBM::interpolateFluidVelocity() {
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
void BlockFreeSurfaceEulerParticlesLBM::normalizeForce() {
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
void BlockFreeSurfaceEulerParticlesLBM::smoothForce() {
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
			BlockFreeSurfaceEulerParticlesLBM* e = (BlockFreeSurfaceEulerParticlesLBM*)b;
			gx += wgt * e->Gx[loc];
			gy += wgt * e->Gy[loc];
			gz += wgt * e->Gz[loc];
		});
		GxSmooth[ind] = gx;
		GySmooth[ind] = gy;
		GzSmooth[ind] = gz;
	}
}

// relaxation times of the pairs of opposite directions (artificial diffusion at the fronts and above the packing limit)
void BlockFreeSurfaceEulerParticlesLBM::setTau() {
	int nPos = (nf - 1) / 2;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] != 1 && flagMap[ind] != 3) {
			continue;
		}
		double cm = concentration[ind];
		for (int p = 0; p < nPos; p++) {
			int pFor = positiveDirs[p];
			int pBack = opp[pFor];
			double cf = cm;
			double cb = cm;
			int bl = nbBlock(ind, pFor);
			int loc = nbLocal(ind, pFor);
			if (bl > -1 && bl < nBlockMax && blocks[bl] && !blocks[bl]->boundary[loc]) {
				cf = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[bl])->concentration[loc];
			}
			bl = nbBlock(ind, pBack);
			loc = nbLocal(ind, pBack);
			if (bl > -1 && bl < nBlockMax && blocks[bl] && !blocks[bl]->boundary[loc]) {
				cb = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[bl])->concentration[loc];
			}
			double tauLoc = computeTau(cf, cm, cb);
			tauDir[nf * ind + pFor] = tauLoc;
			tauDir[nf * ind + pBack] = tauLoc;
		}
	}
}

double BlockFreeSurfaceEulerParticlesLBM::computeTau(double cf, double cm, double cb) {
	double d2 = cf - 2 * cm + cb;
	double overpack = max((cm - concentrationLimit) / concentrationLimit, 0.0);
	return taus + abs(d2) / 1000 + overpack * tauOverPack;
}
