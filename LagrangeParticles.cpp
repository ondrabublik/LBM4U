/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "LagrangeParticles.h"
#include "Block.h"
#include "BlockFreeSurfaceLagrangeParticles.h"

using namespace std;

LagrangeParticles::LagrangeParticles(Block** blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, unordered_map<string, string> par) {
	this->blocks = blocks;
	this->nBlock = nBlock;
	this->nBlockMax = nBlockMax;
	this->nBi = nBi;
	this->nBj = nBj;
	this->nBk = nBk;
	this->ni = ni;
	this->nj = nj;
	this->nk = nk;

	Block b;
	Cu = b.getParDouble("uMax", par, true, 0) / 0.1;
	Ch = b.getParDouble("H", par, true, 0) / (this->ni);
	Ct = Ch / Cu;
	Ca = Ch / (Ct * Ct);
	Crho = b.getParDouble("density", par, true, 0);
	Cf = Crho * Ch * Ch * Cu * Cu;
	G = new double[3]{
		b.getParDouble("gx",par, false, 0) / Ca, b.getParDouble("gy",par, false, 0) / Ca, b.getParDouble("gz",par, false, 0) / Ca
	}; //gravity
	G0 = new double[3]{
		b.getParDouble("gx",par, false, 0) / Ca, b.getParDouble("gy",par, false, 0) / Ca, b.getParDouble("gz",par, false, 0) / Ca
	}; //gravity

	Cvis = Ch * Ch / Ct;
	vLB = b.getParDouble("vKin", par, true, 0) / Cvis;

	interpLevel = 2;
}

LagrangeParticles::~LagrangeParticles() {
	delete partX;
	delete partY;
	delete partZ;

	delete veloX;
	delete veloY;
	delete veloZ;
}

void LagrangeParticles::getParticlesInfo() {
	double kinEnergy = 0;
	for (int i = 0; i < nParticles; i++) {
		kinEnergy += sqrt(veloX[i] * veloX[i] + veloY[i] * veloY[i] + veloZ[i] * veloZ[i]);
	}
	cout << "Kinetic energy: " << kinEnergy << endl;
}

void LagrangeParticles::newParticlePosition(double dtP) {
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < nParticles; i++) {
		double wallDist = wallDistance(i, partX[i], partY[i], partZ[i]);
		if (!inGas[i]) {
			// fluid velocity correction on wall
			if (wallDist < 1) {
				double norm = partWallNormalX[i] * veloFX[i] + partWallNormalY[i] * veloFY[i] + partWallNormalZ[i] * veloFZ[i];
				if (norm < 0) {
					veloFX[i] -= norm * partWallNormalX[i];
					veloFY[i] -= norm * partWallNormalY[i];
					veloFZ[i] -= norm * partWallNormalZ[i];
				}
			}

			// compute interaction force					
			double Uslip = sqrt((veloX[i] - veloFX[i]) * (veloX[i] - veloFX[i]) + (veloY[i] - veloFY[i]) * (veloY[i] - veloFY[i]) + (veloZ[i] - veloFZ[i]) * (veloZ[i] - veloFZ[i]));
			double Cd_mult_Uslip = cd(partD[i], Uslip);
			double areaCorr = 1;
			double coef = 0.5 * Cd_mult_Uslip * partCrossection[i] * areaCorr;

			veloX[i] = (veloX[i] * partMass[i] / dtP + coef * veloFX[i] + G[0] * (partMass[i] - partVolume[i]) + FinterX[i]) / (partMass[i] / dtP + coef);
			veloY[i] = (veloY[i] * partMass[i] / dtP + coef * veloFY[i] + G[1] * (partMass[i] - partVolume[i]) + FinterY[i]) / (partMass[i] / dtP + coef);
			veloZ[i] = (veloZ[i] * partMass[i] / dtP + coef * veloFZ[i] + G[2] * (partMass[i] - partVolume[i]) + FinterZ[i]) / (partMass[i] / dtP + coef);

			Gx[i] = coef * (veloX[i] - veloFX[i]);
			Gy[i] = coef * (veloY[i] - veloFY[i]);
			Gz[i] = coef * (veloZ[i] - veloFZ[i]);

			FinterX[i] = 0;
			FinterY[i] = 0;
			FinterZ[i] = 0;
			nInter[i] = 1;
		}
		else {
			veloX[i] += dtP * G[0];
			veloY[i] += dtP * G[1];
			veloZ[i] += dtP * G[2];
		}

		partX[i] += dtP * veloX[i];
		partY[i] += dtP * veloY[i];
		partZ[i] += dtP * veloZ[i];

		// test of the next possition of the particle
		if (wallDist < partD[i] / 2) {
			double vx0 = veloX[i];
			double vy0 = veloY[i];
			double vz0 = veloZ[i];
			double x0 = partX[i] - dtP * vx0;
			double y0 = partY[i] - dtP * vy0;
			double z0 = partZ[i] - dtP * vz0;
			double norm = partWallNormalX[i] * vx0 + partWallNormalY[i] * vy0 + partWallNormalZ[i] * vz0;

			veloX[i] -= 2 * norm * partWallNormalX[i];
			veloY[i] -= 2 * norm * partWallNormalY[i];
			veloZ[i] -= 2 * norm * partWallNormalZ[i];

			int n = 1;
			while (wallDist < partD[i] / 2) {
				partX[i] = x0 + dtP * (vx0 / n - norm * partWallNormalX[i]);
				partY[i] = y0 + dtP * (vy0 / n - norm * partWallNormalY[i]);
				partZ[i] = z0 + dtP * (vz0 / n - norm * partWallNormalZ[i]);
				wallDist = wallDistance(i, partX[i], partY[i], partZ[i]);

				if (n > 4) {
					partX[i] = x0;
					partY[i] = y0;
					partZ[i] = z0;
					veloX[i] = vx0;
					veloY[i] = vy0;
					veloZ[i] = vz0;
					break;
				}
				n *= 2;
			}
		}
	}
}

void LagrangeParticles::interpolateFluidVelocity() {
#pragma omp parallel for
	for (int p = 0; p < nParticles; p++) {
		// save old velocities
		veloFX[p] = 0;
		veloFY[p] = 0;
		veloFZ[p] = 0;
		pX[p] = 0;
		pY[p] = 0;
		pZ[p] = 0;
		volumeSolid[p] = 0;
		areaSolid[p] = 0;

		// compute block indexes of particle position
		double intX = partX[p] / ni;
		if (intX < 0) {
			intX -= 1;
		}
		double intY = partY[p] / nj;
		if (intY < 0) {
			intY -= 1;
		}
		double intZ = partZ[p] / nk;
		if (intZ < 0) {
			intZ -= 1;
		}
		int iB = (int)intX;
		int jB = (int)intY;
		int kB = (int)intZ;
		int bl = nBi * nBj * kB + nBi * jB + iB;
		if (bl > -1 && bl < nBlockMax && blocks[bl]) {
			int nParticleBox = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->nParticleBox;
			int* particleBox = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->particleBox;
			double* boxVolume = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->boxVolume;
			double* boxArea = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->boxArea;
			double* dpX = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->dpX;
			double* dpY = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->dpY;
			double* dpZ = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->dpZ;
			int* flagMap = ((BlockFreeSurface*)blocks[bl])->flagMap;
			bool* boundary = ((BlockFreeSurface*)blocks[bl])->boundary;
			double* U = blocks[bl]->U;
			double* V = blocks[bl]->V;
			double* W = blocks[bl]->W;

			// compute local particle coordinates on block
			double x = (partX[p] - iB * ni);
			double y = (partY[p] - jB * nj);
			double z = (partZ[p] - kB * nk);

			// compute local indexes in block of nearest cell
			int iloc = (int)x;
			int jloc = (int)y;
			int kloc = (int)z;

			int ind = (ni * nj) * kloc + ni * jloc + iloc;
			if (!boundary[ind]) {
				inGas[p] = true;
				if (flagMap[ind] == 1 || flagMap[ind] == 3) {
					inGas[p] = false;

					double sum = 0;
					for (int i = -interpLevel; i < interpLevel + 1; i++) {
						for (int j = -interpLevel; j < interpLevel + 1; j++) {
							for (int k = -interpLevel; k < interpLevel + 1; k++) {
								int ip = iloc + i;
								int jp = jloc + j;
								int kp = kloc + k;
								double weight = delta(ip + 0.5 - x) * delta(jp + 0.5 - y) * delta(kp + 0.5 - z);
								if (ip > -1 && ip < ni && jp > -1 && jp < nj && kp >-1 && kp < nk) { // inside block
									int ind = (ni * nj) * kp + ni * jp + ip;
									if (flagMap[ind] == 1 || flagMap[ind] == 3) {
										veloFX[p] += weight * U[ind];
										veloFY[p] += weight * V[ind];
										veloFZ[p] += weight * W[ind];
										pX[p] += weight * dpX[ind];
										pY[p] += weight * dpY[ind];
										pZ[p] += weight * dpZ[ind];
										volumeSolid[p] += weight * boxVolume[ind];
										areaSolid[p] += weight * boxArea[ind];
										sum += weight;
									}
								}
								else { // outside block
									int ind = (ni * nj) * kloc + ni * jloc + iloc;
									int indBlockShift = 0;
									int indShift = 0;
									blocks[bl]->local2globalShift(bl, ind, indBlockShift, indShift, i, j, k);
									if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
										int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
										if (fm == 1 || fm == 3) {
											veloFX[p] += weight * ((BlockFreeSurface*)blocks[indBlockShift])->U[indShift];
											veloFY[p] += weight * ((BlockFreeSurface*)blocks[indBlockShift])->V[indShift];
											veloFZ[p] += weight * ((BlockFreeSurface*)blocks[indBlockShift])->W[indShift];
											pX[p] += weight * ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->dpX[indShift];
											pY[p] += weight * ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->dpY[indShift];
											pZ[p] += weight * ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->dpZ[indShift];
											volumeSolid[p] += weight * ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->boxVolume[indShift];
											areaSolid[p] += weight * ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->boxArea[indShift];
											sum += weight;
										}
									}
								}
							}
						}
					}

					veloFX[p] /= sum;
					veloFY[p] /= sum;
					veloFZ[p] /= sum;
					pX[p] /= sum;
					pY[p] /= sum;
					pZ[p] /= sum;
					volumeSolid[p] /= sum;
					areaSolid[p] /= sum;
				}
			}
		}
	}
}

double LagrangeParticles::wallDistance(int p, double xin, double yin, double zin) {

	double L = -1;

	// compute block indexes of particle position
	double intX = xin / ni;
	if (intX < 0) {
		intX -= 1;
	}
	double intY = yin / nj;
	if (intY < 0) {
		intY -= 1;
	}
	double intZ = zin / nk;
	if (intZ < 0) {
		intZ -= 1;
	}
	int iB = (int)intX;
	int jB = (int)intY;
	int kB = (int)intZ;
	int bl = nBi * nBj * kB + nBi * jB + iB;
	if (bl > -1 && bl < nBlockMax && blocks[bl]) {
		// compute local particle coordinates on block
		double x = (xin - iB * ni);
		double y = (yin - jB * nj);
		double z = (zin - kB * nk);

		// compute local indexes in block of nearest cell
		int iloc = (int)x;
		int jloc = (int)y;
		int kloc = (int)z;

		int ind = (ni * nj) * kloc + ni * jloc + iloc;

		bool* boundary = ((BlockFreeSurface*)blocks[bl])->boundary;
		if (!boundary[ind]) {
			double* coeff = blocks[bl]->wallDistancePolynome;

			// wall normal
			double* wd = blocks[bl]->wallDistance;
			double* wdg = blocks[bl]->wallDistanceGradient;
			//L = wd[ind] + wdg[3 * ind] * (x - (iloc + 0.5)) + wdg[3 * ind + 1] * (y - (jloc + 0.5)) + wdg[3 * ind + 2] * (z - (kloc + 0.5));

			int ic = 8 * ind;
			x = x - iloc;
			y = y - jloc;
			z = z - kloc;
			L = coeff[ic] + coeff[ic + 1] * x + coeff[ic + 2] * y + coeff[ic + 3] * z + coeff[ic + 4] * x * y + coeff[ic + 5] * y * z + coeff[ic + 6] * z * x + coeff[ic + 7] * x * y * z;

			partWallNormalX[p] = coeff[ic + 1] + coeff[ic + 4] * y + coeff[ic + 6] * z + coeff[ic + 7] * y * z;
			partWallNormalY[p] = coeff[ic + 2] + coeff[ic + 4] * x + coeff[ic + 5] * z + coeff[ic + 7] * x * z;
			partWallNormalZ[p] = coeff[ic + 3] + coeff[ic + 5] * y + coeff[ic + 6] * x + coeff[ic + 7] * x * y;

			double norm = sqrt(partWallNormalX[p] * partWallNormalX[p] + partWallNormalY[p] * partWallNormalY[p] + partWallNormalZ[p] * partWallNormalZ[p]);

			if (norm > 0) {
				partWallNormalX[p] = partWallNormalX[p] / norm;
				partWallNormalY[p] = partWallNormalY[p] / norm;
				partWallNormalZ[p] = partWallNormalZ[p] / norm;
			}
		}
	}
	return L;
}

void LagrangeParticles::applyForce() {
	for (int p = 0; p < nParticles; p++) {
		double intX = partX[p] / ni;
		if (intX < 0) {
			intX -= 1;
		}
		double intY = partY[p] / nj;
		if (intY < 0) {
			intY -= 1;
		}
		double intZ = partZ[p] / nk;
		if (intZ < 0) {
			intZ -= 1;
		}
		int iB = (int)intX;
		int jB = (int)intY;
		int kB = (int)intZ;
		int bl = nBi * nBj * kB + nBi * jB + iB;
		double x = (partX[p] - iB * ni);
		double y = (partY[p] - jB * nj);
		double z = (partZ[p] - kB * nk);
		int iloc = (int)x;
		int jloc = (int)y;
		int kloc = (int)z;

		double Fx = Gx[p];
		double Fy = Gy[p];
		double Fz = Gz[p];

		if (bl > -1 && bl < nBlockMax && blocks[bl]) {
			double* particleFx = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->particleForceX;
			double* particleFy = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->particleForceY;
			double* particleFz = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->particleForceZ;
			int* flagMap = ((BlockFreeSurface*)blocks[bl])->flagMap;
			for (int i = -interpLevel; i < interpLevel + 1; i++) {
				for (int j = -interpLevel; j < interpLevel + 1; j++) {
					for (int k = -interpLevel; k < interpLevel + 1; k++) {
						int ip = iloc + i;
						int jp = jloc + j;
						int kp = kloc + k;
						if (ip > -1 && ip < ni && jp > -1 && jp < nj && kp > -1 && kp < nk) { // inside block
							int ind = (ni * nj) * kp + ni * jp + ip;
							if (flagMap[ind] == 1 || flagMap[ind] == 3) {
								double weight = delta(iloc + 0.5 - x + i) * delta(jloc + 0.5 - y + j) * delta(kloc + 0.5 - z + k);
								particleFx[ind] += weight * Fx;
								particleFy[ind] += weight * Fy;
								particleFz[ind] += weight * Fz;
							}
						}
						else { // outside block
							int ind = (ni * nj) * kloc + ni * jloc + iloc;
							int indBlockShift = 0;
							int indShift = 0;
							blocks[bl]->local2globalShift(bl, ind, indBlockShift, indShift, i, j, k);
							if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
								int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
								if (fm == 1 || fm == 3) {
									double weight = delta(iloc + 0.5 - x + i) * delta(jloc + 0.5 - y + j) * delta(kloc + 0.5 - z + k);
									((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->particleForceX[indShift] += weight * Fx;
									((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->particleForceY[indShift] += weight * Fy;
									((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->particleForceZ[indShift] += weight * Fz;
								}
							}
						}
					}
				}
			}
		}
	}
}

void LagrangeParticles::boxParticles() {
	// put particles into cells box
	for (int p = 0; p < nParticles; p++) {
		double intX = partX[p] / ni;
		if (intX < 0) {
			intX -= 1;
		}
		double intY = partY[p] / nj;
		if (intY < 0) {
			intY -= 1;
		}
		double intZ = partZ[p] / nk;
		if (intZ < 0) {
			intZ -= 1;
		}
		int iB = (int)intX;
		int jB = (int)intY;
		int kB = (int)intZ;
		int bl = nBi * nBj * kB + nBi * jB + iB;
		if (bl > -1 && bl < nBlockMax && blocks[bl]) { // particle is inside block
			int nParticleBox = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->nParticleBox;
			int* particleBox = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->particleBox;
			double* boxVolume = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->boxVolume;
			double* boxArea = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->boxArea;

			// compute local particle coordinates on block
			double x = (partX[p] - iB * ni);
			double y = (partY[p] - jB * nj);
			double z = (partZ[p] - kB * nk);

			// compute local indexes in block of nearest cell
			int iloc = (int)x;
			int jloc = (int)y;
			int kloc = (int)z;

			int ind = (ni * nj) * kloc + ni * jloc + iloc;
			int indp = ind * (nParticleBox + 1);
			if (particleBox[indp] < nParticleBox) {
				particleBox[indp]++;
				particleBox[indp + particleBox[indp]] = p;
				boxVolume[ind] += partVolume[p];
				boxArea[ind] += partCrossection[p];
			}
			else {
				cout << "particle box overflow" << endl;
			}
		}
	}
}

void LagrangeParticles::particleWallInteraction() {
	for (int p = 0; p < nParticles; p++) {
		double L = wallDistance(p, partX[p], partY[p], partZ[p]);
		double wallGap = partD[p];
		if (L > 0 && L < wallGap) {
			double xi = 1 - abs(L - partD[p] / 2) / (wallGap - partD[p] / 2);
			double Fwall = 0.05 * xi;
			FwallX[p] = Fwall * partWallNormalX[p];
			FwallY[p] = Fwall * partWallNormalY[p];
			FwallZ[p] = Fwall * partWallNormalZ[p];
		}
		else {
			FwallX[p] = 0;
			FwallY[p] = 0;
			FwallZ[p] = 0;
		}

	}
}

void LagrangeParticles::particleParticleInteraction() {
	// get particles from surrounding boxes
#pragma omp parallel for
	for (int p = 0; p < nParticles; p++) {
		double intX = partX[p] / ni;
		if (intX < 0) {
			intX -= 1;
		}
		double intY = partY[p] / nj;
		if (intY < 0) {
			intY -= 1;
		}
		double intZ = partZ[p] / nk;
		if (intZ < 0) {
			intZ -= 1;
		}
		int iB = (int)intX;
		int jB = (int)intY;
		int kB = (int)intZ;
		int indBlock = nBi * nBj * kB + nBi * jB + iB;
		if (indBlock > -1 && indBlock < nBlockMax && blocks[indBlock]) { // particle is inside of the block
			// compute local particle coordinates on block
			double x = (partX[p] - iB * ni);
			double y = (partY[p] - jB * nj);
			double z = (partZ[p] - kB * nk);

			// compute local indexes in block of nearest cell
			int iloc = (int)x;
			int jloc = (int)y;
			int kloc = (int)z;
			int ind = (ni * nj) * kloc + ni * jloc + iloc;

			int indBlockShift = 0;
			int indShift = 0;
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					for (int k = -1; k <= 1; k++) {
						blocks[indBlock]->local2globalShift(indBlock, ind, indBlockShift, indShift, i, j, k);
						if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) { // block exists
							int nParticleBox = ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->nParticleBox;
							int indShiftp = indShift * (nParticleBox + 1);
							for (int ps = 0; ps < ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->particleBox[indShiftp]; ps++) {
								int q = ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->particleBox[indShiftp + ps + 1];
								//if (p > q) { // for simulation speed up
								double nx = partX[p] - partX[q];
								double ny = partY[p] - partY[q];
								double nz = partZ[p] - partZ[q];
								double L = sqrt(nx * nx + ny * ny + nz * nz);
								if (L > 1e-10) {
									nx /= L;
									ny /= L;
									nz /= L;
								}
								else {
									nx = 0;
									ny = 0;
									nz = 0;
								}
								double uRel = nx * (veloX[p] - veloX[q]) + ny * (veloY[p] - veloY[q]) + nz * (veloZ[p] - veloZ[q]);
								double avgD = (partD[p] + partD[q]) / 2;
								double f = gammaPP(L, avgD, uRel);
								double fx = f * nx;
								double fy = f * ny;
								double fz = f * nz;

								FinterX[p] += fx;
								FinterY[p] += fy;
								FinterZ[p] += fz;
								nInter[p]++;

								//FinterX[q] -= fx;
								//FinterY[q] -= fy;
								//FinterZ[q] -= fz;
							//}
							}
						}
					}
				}
			}
		}
	}
}

double LagrangeParticles::delta(double x) {
	x = abs(x);
	double y = 0;
	if (x <= 2) {
		y = (1 + cos(3.1415 * x / 2.0)) / 4.0;
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
	// get particles from surrounding boxes
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

void LagrangeParticles::readInitParticlesDistribution(string simulPath, double Ch, double Cu) {
	string line;
	string fileName = simulPath + "init/particles.txt";
	try {
		ifstream* myfile = new ifstream(fileName);
		if (myfile->is_open()) {
			getline(*myfile, line);
			nParticles = stoi(line);
			partD = new double[nParticles];
			partRho = new double[nParticles];
			partMass = new double[nParticles];
			partVolume = new double[nParticles];
			partCrossection = new double[nParticles];
			partX = new double[nParticles];
			partY = new double[nParticles];
			partZ = new double[nParticles];
			veloX = new double[nParticles];
			veloY = new double[nParticles];
			veloZ = new double[nParticles];
			veloFX = new double[nParticles];
			veloFY = new double[nParticles];
			veloFZ = new double[nParticles];
			veloX0 = new double[nParticles];
			veloY0 = new double[nParticles];
			veloZ0 = new double[nParticles];
			Gx = new double[nParticles];
			Gy = new double[nParticles];
			Gz = new double[nParticles];
			pX = new double[nParticles];
			pY = new double[nParticles];
			pZ = new double[nParticles];
			FinterX = new double[nParticles];
			FinterY = new double[nParticles];
			FinterZ = new double[nParticles];
			FwallX = new double[nParticles];
			FwallY = new double[nParticles];
			FwallZ = new double[nParticles];
			volumeSolid = new double[nParticles];
			areaSolid = new double[nParticles];
			nInter = new int[nParticles];
			partWallNormalX = new double[nParticles];
			partWallNormalY = new double[nParticles];
			partWallNormalZ = new double[nParticles];
			partWallDistance = new double[nParticles];
			inGas = new bool[nParticles];

			for (int i = 0; i < nParticles; i++) {
				getline(*myfile, line);
				istringstream iss(line);
				vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
				partX[i] = stod(tokens[0]) / Ch;
				partY[i] = stod(tokens[1]) / Ch;
				partZ[i] = stod(tokens[2]) / Ch;
				partD[i] = stod(tokens[3]) / Ch;
				partRho[i] = stod(tokens[4]) / Crho;
				partCrossection[i] = 3.1415 / 4 * partD[i] * partD[i];
				partVolume[i] = 3.1415 / 6 * partD[i] * partD[i] * partD[i];
				partMass[i] = partVolume[i] * partRho[i];
				pX[i] = 0;
				pY[i] = 0;
				pZ[i] = 0;
				veloX[i] = 0;
				veloY[i] = 0;
				veloZ[i] = 0;
				veloX0[i] = 0;
				veloY0[i] = 0;
				veloZ0[i] = 0;
				FinterX[i] = 0;
				FinterY[i] = 0;
				FinterZ[i] = 0;
				nInter[i] = 1;
				inGas[i] = false;
			}
		}
		myfile->close();
		delete myfile;
	}
	catch (const std::exception& ex) {
		cout << "Error in " << fileName << ": " << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	cout << nParticles << " Lagrangian particles were readed!" << endl;
}

void LagrangeParticles::saveParticlesPositionsToVTK(string simulPath, int iter) {
	ofstream myfile;
	myfile.open(simulPath + "result/output_particles" + std::to_string((int)(iter * Ct * 1000000)) + "us.vtk", ios::out);
	myfile << "# vtk DataFile Version 2.6\noutput file\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	myfile << "POINTS " << nParticles << " float\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << partX[i] << " " << partY[i] << " " << partZ[i] << "\n";
	}
	myfile << "\n";

	int s = 0;
	myfile << "CELLS " << nParticles << " " << 2 * nParticles << "\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << "1 " << s << "\n";
		s++;
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

	myfile.close();
}
