/*
 * File:   LagrangeParticles.cpp
 * Author: obublik
 */

#include "LagrangeParticles.h"
#include "Block.h"
#include "BlockFreeSurfaceLagrangeParticles.h"
#include <cmath>
#include <stdexcept>

using namespace std;

typedef BlockFreeSurfaceLagrangeParticles LagBlock;

namespace {
	const double PI = 3.14159265358979323846;
}

LagrangeParticles::LagrangeParticles(Block** blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, const ParamMap& par) {
	this->blocks = blocks;
	this->nBlock = nBlock;
	this->nBlockMax = nBlockMax;
	this->nBi = nBi;
	this->nBj = nBj;
	this->nBk = nBk;
	this->ni = ni;
	this->nj = nj;
	this->nk = nk;

	Cu = Block::getParDouble("uMax", par, true, 0) / 0.1;
	Ch = Block::getParDouble("H", par, true, 0) / (this->ni);
	Ct = Ch / Cu;
	Ca = Ch / (Ct * Ct);
	Crho = Block::getParDouble("density", par, true, 0);
	Cf = Crho * Ch * Ch * Cu * Cu;
	G[0] = G0[0] = Block::getParDouble("gx", par, false, 0) / Ca;
	G[1] = G0[1] = Block::getParDouble("gy", par, false, 0) / Ca;
	G[2] = G0[2] = Block::getParDouble("gz", par, false, 0) / Ca;

	Cvis = Ch * Ch / Ct;
	vLB = Block::getParDouble("vKin", par, true, 0) / Cvis;
	restitution = min(max(Block::getParDouble("particleRestitution", par, false, 1.0), 0.0), 1.0);

	interpLevel = 2;
}

void LagrangeParticles::resize(int n) {
	nParticles = n;
	for (vector<double>* v : { &partD, &partRho, &partMass, &partVolume, &partCrossection, &partX, &partY, &partZ,
		&veloX, &veloY, &veloZ, &veloFX, &veloFY, &veloFZ, &pX, &pY, &pZ, &Gx, &Gy, &Gz, &FinterX, &FinterY, &FinterZ,
		&volumeSolid, &areaSolid, &partWallNormalX, &partWallNormalY, &partWallNormalZ }) {
		v->assign(n, 0.0);
	}
	nInter.assign(n, 1);
	inGas.assign(n, 0);
}

// kinetic energy of the particles [J]
void LagrangeParticles::getParticlesInfo() {
	double kinEnergy = 0;
	for (int i = 0; i < nParticles; i++) {
		kinEnergy += 0.5 * partMass[i] * (veloX[i] * veloX[i] + veloY[i] * veloY[i] + veloZ[i] * veloZ[i]);
	}
	cout << "Kinetic energy of the particles: " << kinEnergy * Crho * Ch * Ch * Ch * Cu * Cu << " J" << endl;
}

bool LagrangeParticles::locate(double x, double y, double z, Cell& c) const {
	if (!(std::isfinite(x) && std::isfinite(y) && std::isfinite(z))) {
		return false;
	}
	// refined mesh: the finest computed block holding the point (the coordinates are in the cells of the level 0)
	int maxL = grid ? grid->maxLevel : 0;
	for (int L = maxL; L >= 0; L--) {
		double s = (double)(1 << L);
		double fb[3] = { floor(x * s / ni), floor(y * s / nj), floor(z * s / nk) };
		bool big = false;
		for (double v : fb) {
			big = big || fabs(v) > 1e9;
		}
		if (big) {
			return false;
		}
		int iB = (int)fb[0], jB = (int)fb[1], kB = (int)fb[2];
		if (grid) {
			c.bl = grid->find(L, iB, jB, kB); // wraps the periodic directions
		}
		else {
			bool in = iB >= 0 && iB < nBi && jB >= 0 && jB < nBj && kB >= 0 && kB < nBk;
			c.bl = in ? (nBi * nBj) * kB + nBi * jB + iB : -1;
		}
		if (c.bl < 0 || c.bl >= nBlockMax || !blocks[c.bl] || blocks[c.bl]->ghost != GHOST_NONE) {
			continue;
		}
		c.level = L;
		c.x = x * s - fb[0] * ni;
		c.y = y * s - fb[1] * nj;
		c.z = z * s - fb[2] * nk;
		c.i = min(max((int)c.x, 0), ni - 1);
		c.j = min(max((int)c.y, 0), nj - 1);
		c.k = min(max((int)c.z, 0), nk - 1);
		return true;
	}
	return false;
}

template <class F> void LagrangeParticles::liquidCells(const Cell& c, F fn) {
	Block* b0 = blocks[c.bl];
	int indC = (ni * nj) * c.k + ni * c.j + c.i;
	for (int k = -interpLevel; k <= interpLevel; k++) {
		double wk = delta(c.k + k + 0.5 - c.z);
		if (wk == 0) {
			continue;
		}
		for (int j = -interpLevel; j <= interpLevel; j++) {
			double wj = wk * delta(c.j + j + 0.5 - c.y);
			if (wj == 0) {
				continue;
			}
			for (int i = -interpLevel; i <= interpLevel; i++) {
				double weight = wj * delta(c.i + i + 0.5 - c.x);
				if (weight == 0) {
					continue;
				}
				int ip = c.i + i, jp = c.j + j, kp = c.k + k;
				int bl = c.bl, ind;
				if (ip > -1 && ip < ni && jp > -1 && jp < nj && kp > -1 && kp < nk) {
					ind = (ni * nj) * kp + ni * jp + ip;
				}
				else {
					b0->local2globalShift(c.bl, indC, bl, ind, i, j, k);
					if (bl < 0 || bl >= nBlockMax || !blocks[bl]) {
						continue;
					}
				}
				LagBlock* b = (LagBlock*)blocks[bl];
				// ghost blocks (copies of another level near a level interface) are not used
				if (b->ghost == GHOST_NONE && !b->boundary[ind] && (b->flagMap[ind] == 1 || b->flagMap[ind] == 3)) {
					fn(b, ind, weight);
				}
			}
		}
	}
}

void LagrangeParticles::forLiquidCells(const Cell& c, const function<void(LagBlock*, int, double)>& fn) {
	liquidCells(c, fn);
}

// one step of the equation of motion (implicit drag) and of the position, dtP is a fraction of the LBM step
void LagrangeParticles::newParticlePosition(double dtP) {
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nParticles; i++) {
		double x0 = partX[i], y0 = partY[i], z0 = partZ[i];
		double L0 = wallDistance(i, x0, y0, z0);
		double n0[3] = { partWallNormalX[i], partWallNormalY[i], partWallNormalZ[i] };
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

			double Uslip = sqrt((veloX[i] - veloFX[i]) * (veloX[i] - veloFX[i]) + (veloY[i] - veloFY[i]) * (veloY[i] - veloFY[i]) + (veloZ[i] - veloFZ[i]) * (veloZ[i] - veloFZ[i]));
			double coef = 0.5 * cd(partD[i], Uslip) * partCrossection[i];

			veloX[i] = (veloX[i] * partMass[i] / dtP + coef * veloFX[i] + G[0] * (partMass[i] - partVolume[i]) + FinterX[i]) / (partMass[i] / dtP + coef);
			veloY[i] = (veloY[i] * partMass[i] / dtP + coef * veloFY[i] + G[1] * (partMass[i] - partVolume[i]) + FinterY[i]) / (partMass[i] / dtP + coef);
			veloZ[i] = (veloZ[i] * partMass[i] / dtP + coef * veloFZ[i] + G[2] * (partMass[i] - partVolume[i]) + FinterZ[i]) / (partMass[i] / dtP + coef);

			// reaction of the drag: impulse given to the fluid, summed over the sub-steps of the LBM step
			Gx[i] += dtP * coef * (veloX[i] - veloFX[i]);
			Gy[i] += dtP * coef * (veloY[i] - veloFY[i]);
			Gz[i] += dtP * coef * (veloZ[i] - veloFZ[i]);
		}
		else {
			veloX[i] += dtP * G[0];
			veloY[i] += dtP * G[1];
			veloZ[i] += dtP * G[2];
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

// The particle (a sphere of the radius d/2) must not penetrate the wall: it is moved back along the wall
// normal to the distance d/2 and the velocity component towards the wall is reflected. The tangential
// motion is kept, so the particle can slide along the wall.
void LagrangeParticles::wallContact(int p, double x0, double y0, double z0, double L0, const double* n0) {
	double L = wallDistance(p, partX[p], partY[p], partZ[p]); // also sets the wall normal
	double r = contactRadius(p);
	bool contact = false;
	for (int it = 0; it < 4 && L >= 0 && L < r; it++) {
		// the distance function is not exact near the edges, a few corrections are enough
		contact = true;
		partX[p] += (r - L) * partWallNormalX[p];
		partY[p] += (r - L) * partWallNormalY[p];
		partZ[p] += (r - L) * partWallNormalZ[p];
		L = wallDistance(p, partX[p], partY[p], partZ[p]);
		r = contactRadius(p);
	}
	if (L >= 0) {
		if (contact) {
			double nx = partWallNormalX[p], ny = partWallNormalY[p], nz = partWallNormalZ[p];
			double vn = nx * veloX[p] + ny * veloY[p] + nz * veloZ[p];
			if (vn < 0) {
				veloX[p] -= (1 + restitution) * vn * nx;
				veloY[p] -= (1 + restitution) * vn * ny;
				veloZ[p] -= (1 + restitution) * vn * nz;
			}
		}
		return;
	}
	// the step ended inside of a wall or outside of the domain: back to the previous position
	partX[p] = x0;
	partY[p] = y0;
	partZ[p] = z0;
	partWallNormalX[p] = n0[0];
	partWallNormalY[p] = n0[1];
	partWallNormalZ[p] = n0[2];
	if (L0 < 0) {
		// the particle was placed outside of the fluid, it stays there
		veloX[p] = veloY[p] = veloZ[p] = 0;
		return;
	}
	double vn = n0[0] * veloX[p] + n0[1] * veloY[p] + n0[2] * veloZ[p];
	if (vn < 0) {
		veloX[p] -= (1 + restitution) * vn * n0[0];
		veloY[p] -= (1 + restitution) * vn * n0[1];
		veloZ[p] -= (1 + restitution) * vn * n0[2];
	}
	else if (vn == 0) {
		veloX[p] = veloY[p] = veloZ[p] = 0; // no normal known, stop the particle
	}
}

// fluid velocity (and the other fields) at the particles, starts the LBM step of the particles
void LagrangeParticles::interpolateFluidVelocity() {
#pragma omp parallel for
	for (int p = 0; p < nParticles; p++) {
		veloFX[p] = 0;
		veloFY[p] = 0;
		veloFZ[p] = 0;
		pX[p] = 0;
		pY[p] = 0;
		pZ[p] = 0;
		volumeSolid[p] = 0;
		areaSolid[p] = 0;
		Gx[p] = 0;
		Gy[p] = 0;
		Gz[p] = 0;

		Cell c;
		if (!locate(partX[p], partY[p], partZ[p], c)) {
			inGas[p] = 1; // outside of the domain: no coupling with the fluid
			continue;
		}
		LagBlock* b0 = (LagBlock*)blocks[c.bl];
		int ind0 = (ni * nj) * c.k + ni * c.j + c.i;
		if (b0->boundary[ind0]) {
			continue;
		}
		inGas[p] = !(b0->flagMap[ind0] == 1 || b0->flagMap[ind0] == 3);
		if (inGas[p]) {
			continue;
		}
		double sum = 0, u = 0, v = 0, w = 0, px = 0, py = 0, pz = 0, vol = 0, area = 0;
		liquidCells(c, [&](LagBlock* b, int ind, double weight) {
			u += weight * b->U[ind];
			v += weight * b->V[ind];
			w += weight * b->W[ind];
			px += weight * b->dpX[ind];
			py += weight * b->dpY[ind];
			pz += weight * b->dpZ[ind];
			vol += weight * b->boxVolume[ind];
			area += weight * b->boxArea[ind];
			sum += weight;
		});
		if (sum > 0) {
			veloFX[p] = u / sum;
			veloFY[p] = v / sum;
			veloFZ[p] = w / sum;
			pX[p] = px / sum;
			pY[p] = py / sum;
			pZ[p] = pz / sum;
			volumeSolid[p] = vol / sum;
			areaSolid[p] = area / sum;
		}
	}
}

double LagrangeParticles::wallDistance(int p, double xin, double yin, double zin) {
	double L = -1;
	Cell c;
	if (!locate(xin, yin, zin, c)) {
		return L;
	}
	Block* b = blocks[c.bl];
	int ind = (ni * nj) * c.k + ni * c.j + c.i;
	if (b->boundary[ind]) {
		return L;
	}
	// trilinear interpolation of the wall distance in the cell
	const double* coeff = &b->wallDistancePolynome[8 * ind];
	double x = c.x - c.i;
	double y = c.y - c.j;
	double z = c.z - c.k;
	L = coeff[0] + coeff[1] * x + coeff[2] * y + coeff[3] * z + coeff[4] * x * y + coeff[5] * y * z + coeff[6] * z * x + coeff[7] * x * y * z;

	double nx = coeff[1] + coeff[4] * y + coeff[6] * z + coeff[7] * y * z;
	double ny = coeff[2] + coeff[4] * x + coeff[5] * z + coeff[7] * x * z;
	double nz = coeff[3] + coeff[5] * y + coeff[6] * x + coeff[7] * x * y;
	double norm = sqrt(nx * nx + ny * ny + nz * nz);
	if (norm > 0) {
		nx /= norm;
		ny /= norm;
		nz /= norm;
	}
	partWallNormalX[p] = nx;
	partWallNormalY[p] = ny;
	partWallNormalZ[p] = nz;
	return max(L, 0.0) / (1 << c.level); // in the cells of the level 0
}

// Spreads the impulse of the particles to the liquid cells around. The cell gets weight / sum(weight * vof)
// of it and uses it as force * vof / rho (cellForce), so the whole impulse reaches the fluid also next to
// the walls and to the free surface. The force acts in all steps of the level until the next call: a block of
// the level L makes 2^L steps per particle step, its momentum unit is 8^-L of the level 0 one, so it gets the
// impulse * 8^L / 2^L = 4^L per step.
void LagrangeParticles::applyForce() {
#pragma omp parallel for schedule(dynamic)
	for (int s = 0; s < nBlockMax; s++) {
		if (!blocks[s] || blocks[s]->ghost != GHOST_NONE) {
			continue;
		}
		LagBlock* b = (LagBlock*)blocks[s];
		for (int ind = 0; ind < b->ntot; ind++) {
			b->particleForceX[ind] = 0;
			b->particleForceY[ind] = 0;
			b->particleForceZ[ind] = 0;
		}
	}
	for (int p = 0; p < nParticles; p++) {
		double Fx = Gx[p];
		double Fy = Gy[p];
		double Fz = Gz[p];
		if (Fx == 0 && Fy == 0 && Fz == 0) {
			continue;
		}
		Cell c;
		if (!locate(partX[p], partY[p], partZ[p], c)) {
			continue;
		}
		double sum = 0;
		liquidCells(c, [&](LagBlock* b, int ind, double weight) {
			sum += weight * min(max(b->vof[ind], 0.0), 1.0);
		});
		if (sum <= 0) {
			continue;
		}
		sum = max(sum, BlockFreeSurface::FORCE_KERNEL_MIN) / (double)(1 << (2 * c.level));
		liquidCells(c, [&](LagBlock* b, int ind, double weight) {
			b->particleForceX[ind] += weight / sum * Fx;
			b->particleForceY[ind] += weight / sum * Fy;
			b->particleForceZ[ind] += weight / sum * Fz;
		});
	}
}

// lists of the particles in the cells (collisions, concentration)
void LagrangeParticles::boxParticles() {
#pragma omp parallel for schedule(dynamic)
	for (int s = 0; s < nBlockMax; s++) {
		if (!blocks[s]) {
			continue;
		}
		LagBlock* b = (LagBlock*)blocks[s];
		for (int ind = 0; ind < b->ntot; ind++) {
			b->particleBox[ind * (b->nParticleBox + 1)] = 0;
			b->boxVolume[ind] = 0;
			b->boxArea[ind] = 0;
		}
	}
	int overflow = 0;
	for (int p = 0; p < nParticles; p++) {
		Cell c;
		if (!locate(partX[p], partY[p], partZ[p], c)) {
			continue;
		}
		LagBlock* b = (LagBlock*)blocks[c.bl];
		int ind = (ni * nj) * c.k + ni * c.j + c.i;
		int indp = ind * (b->nParticleBox + 1);
		if (b->particleBox[indp] < b->nParticleBox) {
			b->particleBox[indp]++;
			b->particleBox[indp + b->particleBox[indp]] = p;
			b->boxVolume[ind] += partVolume[p];
			b->boxArea[ind] += partCrossection[p];
		}
		else {
			overflow++;
		}
	}
	static bool warned = false;
	if (overflow > 0 && !warned) {
		warned = true;
		cout << "warning: particle box overflow (" << overflow << " particles), increase maxParticlesInCell" << endl;
	}
}

// repulsive forces between the particles of the neighbouring cells
void LagrangeParticles::particleParticleInteraction() {
#pragma omp parallel for
	for (int p = 0; p < nParticles; p++) {
		Cell c;
		if (!locate(partX[p], partY[p], partZ[p], c)) {
			continue;
		}
		int ind = (ni * nj) * c.k + ni * c.j + c.i;
		// the cells of the level have to cover the range 1.5 d of the force
		int r = max(1, (int)ceil(1.5 * dMax * (1 << c.level)));
		for (int i = -r; i <= r; i++) {
			for (int j = -r; j <= r; j++) {
				for (int k = -r; k <= r; k++) {
					int bl, loc;
					blocks[c.bl]->local2globalShift(c.bl, ind, bl, loc, i, j, k);
					if (bl < 0 || bl >= nBlockMax || !blocks[bl]) {
						continue;
					}
					LagBlock* b = (LagBlock*)blocks[bl];
					const int* box = &b->particleBox[loc * (b->nParticleBox + 1)];
					for (int ps = 0; ps < box[0]; ps++) {
						int q = box[ps + 1];
						if (q == p) {
							continue;
						}
						double nx = partX[p] - partX[q];
						double ny = partY[p] - partY[q];
						double nz = partZ[p] - partZ[q];
						double L = sqrt(nx * nx + ny * ny + nz * nz);
						if (L <= 1e-10) {
							continue;
						}
						nx /= L;
						ny /= L;
						nz /= L;
						double uRel = nx * (veloX[p] - veloX[q]) + ny * (veloY[p] - veloY[q]) + nz * (veloZ[p] - veloZ[q]);
						double avgD = (partD[p] + partD[q]) / 2;
						double f = gammaPP(L, avgD, uRel);
						FinterX[p] += f * nx;
						FinterY[p] += f * ny;
						FinterZ[p] += f * nz;
						nInter[p]++;
					}
				}
			}
		}
	}
}

// 4-point cosine delta function (Peskin)
double LagrangeParticles::delta(double x) {
	x = abs(x);
	double y = 0;
	if (x < 2) {
		y = (1 + cos(PI * x / 2.0)) / 4.0;
	}
	return y;
}

double LagrangeParticles::cd(double d, double uSlip) { // return Cd*uSlip
	double CdUslip = 0.44 * uSlip;
	double Re = d * uSlip / vLB;
	if (Re < 1000) {
		CdUslip = 24.0 / d * vLB * (1 + 0.15 * pow(Re, 0.687));
	}
	return CdUslip;
}

double LagrangeParticles::gammaPP(double L, double avgD, double uRel) {
	double y = 0;
	double Lint = 1.5 * avgD;
	if (L < Lint) {
		y = 10 / Ca * sqrt(Lint - L);
	}
	return y;
}

double LagrangeParticles::gammaWall(double x, double D) {
	x = abs(x);
	double y = 0;
	if (x <= 2 * D) {
		y = 1 * (1 - x / (2 * D));
	}
	return y;
}

void LagrangeParticles::particleBarrierForce(int iter) {
	// temporary barrier holding the particles at the start of the simulation (switched off: tRemove = 0)
	double t = iter * Ct;
	double yBarier = 0.11;
	double tRemove = 0;
	if (t < tRemove) {
#pragma omp parallel for
		for (int p = 0; p < nParticles; p++) {
			double dy = Ch * partY[p] - yBarier;
			if (dy < 0) {
				FinterY[p] += gammaWall(dy, 0.1);
			}
		}
	}
}

// init/particles.txt: the number of particles and the lines "x y z diameter density" in physical units
void LagrangeParticles::readInitParticlesDistribution(string simulPath, double Ch, double Cu) {
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
	resize(n);
	for (int i = 0; i < n; i++) {
		double x, y, z, d, rho;
		if (!getline(myfile, line) || !(istringstream(line) >> x >> y >> z >> d >> rho) || !(d > 0) || !(rho > 0)) {
			throw runtime_error("Error in " + fileName + ": particle " + to_string(i + 1) + " needs x y z diameter density");
		}
		partX[i] = x / Ch;
		partY[i] = y / Ch;
		partZ[i] = z / Ch;
		partD[i] = d / Ch;
		partRho[i] = rho / Crho;
		partCrossection[i] = PI / 4 * partD[i] * partD[i];
		partVolume[i] = PI / 6 * partD[i] * partD[i] * partD[i];
		partMass[i] = partVolume[i] * partRho[i];
		dMax = max(dMax, partD[i]);
	}
	cout << nParticles << " Lagrangian particles were readed!" << endl;
}

// positions in the cells of the level 0 (as the mesh in output_*.vtu)
void LagrangeParticles::saveParticlesPositionsToVTK(string simulPath, int iter) {
	ofstream myfile;
	myfile.open(simulPath + "result/output_particles" + std::to_string((long long)(iter * Ct * 1000000)) + "us.vtk", ios::out);
	myfile << "# vtk DataFile Version 2.6\noutput file\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	myfile << "POINTS " << nParticles << " float\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << partX[i] << " " << partY[i] << " " << partZ[i] << "\n";
	}
	myfile << "\n";

	myfile << "CELLS " << nParticles << " " << 2 * nParticles << "\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << "1 " << i << "\n";
	}
	myfile << "\n";

	myfile << "CELL_TYPES " << nParticles << "\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << "1 " << "\n";
	}
	myfile << "\n";

	myfile << "CELL_DATA " << nParticles << "\n";
	myfile << "SCALARS densityRatio float 1\n";
	myfile << "LOOKUP_TABLE default\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << partRho[i] << "\n";
	}
	myfile << "\n";

	myfile << "SCALARS diameter float 1\n";
	myfile << "LOOKUP_TABLE default\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << partD[i] << "\n";
	}
	myfile << "\n";

	myfile << "SCALARS particleDensity float 1\n";
	myfile << "LOOKUP_TABLE default\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << volumeSolid[i] << "\n";
	}
	myfile << "\n";

	myfile << "VECTORS rychlost float\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << veloX[i] * Cu << " " << veloY[i] * Cu << " " << veloZ[i] * Cu << "\n";
	}
	myfile << "\n";

	myfile << "VECTORS gradP float\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << pX[i] << " " << pY[i] << " " << pZ[i] << "\n";
	}
	myfile << "\n";

	myfile << "VECTORS force float\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << Gx[i] * Cf << " " << Gy[i] * Cf << " " << Gz[i] * Cf << "\n";
	}
	myfile << "\n";
	writeExtraVTK(myfile);

	myfile.close();
}
