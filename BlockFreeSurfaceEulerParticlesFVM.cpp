
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 /*
  * File:   Block.cpp
  * Author: obublik
  *
  * Created on 2. února 2018, 12:27
  */

#include "BlockFreeSurfaceEulerParticlesFVM.h"

using namespace std;

BlockFreeSurfaceEulerParticlesFVM::BlockFreeSurfaceEulerParticlesFVM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurface(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, coarseLevel, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {

	concentration = new double[ntot]();
	fluxConcentration = new double[ntot]();

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

	fcP = new double[ntot * nf]();

	dParticle = getParDouble("particleDiameter", par, true, 0) / Ch;
	rhoParticle = getParDouble("particleDensity", par, true, 0) / Crho;

	double pi = 3.14159265358979323846;
	crossectionParticle = pi * dParticle * dParticle / 4;
	volumeParticle = pi * dParticle * dParticle * dParticle / 6;
	massParticle = rhoParticle * volumeParticle;
	Cmax = 1 / (dParticle * dParticle * dParticle);
	interpLevel = 2;

	exLW = new int[6]{ 1, 0, 0, -1,  0,  0 };
	eyLW = new int[6]{ 0, 1, 0,  0, -1,  0 };
	ezLW = new int[6]{ 0, 0, 1,  0,  0, -1 };
	oppLW = new int[6]{ 3, 4, 5, 0, 1, 2 };
}

BlockFreeSurfaceEulerParticlesFVM::~BlockFreeSurfaceEulerParticlesFVM() {
	delete[] fcP;
	delete[] concentration;
	delete[] fluxConcentration;
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
	delete[] exLW;
	delete[] eyLW;
	delete[] ezLW;
	delete[] oppLW;
	delete[] deltaValue;
}

void BlockFreeSurfaceEulerParticlesFVM::readInitConcentration(int id, string simulPath) {
	string line;
	string fileName = simulPath + "init/concentration" + std::to_string(id) + ".txt";
	try {
		ifstream* myfile = new ifstream(fileName);
		if (myfile->is_open()) {
			for (int i = 0; i < nAll; i++) {
				int ind = indAll[i];

				getline(*myfile, line);
				istringstream iss(line);
				vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
				concentration[ind] = stod(tokens[0]);
			}
		}
		myfile->close();
		delete myfile;
	}
	catch (const std::exception& ex) {
		cout << "Error in " << fileName << ": " << ex.what() << endl;
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < ntot; i++) {
		UP[i] = U[i];
		VP[i] = V[i];
		WP[i] = W[i];
	}

	// values of delta function
	deltaValue = new double[3];
	for (int i = 0; i < 3; i++) {
		deltaValue[i] = delta((double)i);
	}
}

void BlockFreeSurfaceEulerParticlesFVM::cellForce(int ind, double* acc, double& epsF) {
	// smoothed drag is a force density, the exact difference method needs an acceleration
	double rho = Rho[ind] > 0 ? Rho[ind] : 1.0;
	epsF = 1.0;
	acc[0] = G[0] + GxSmooth[ind] * vof[ind] / rho;
	acc[1] = G[1] + GySmooth[ind] * vof[ind] / rho;
	acc[2] = G[2] + GzSmooth[ind] * vof[ind] / rho;
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
	double cA = concentration[ind];
	double cB = bs->concentration[loc];
	double Un = (bs->UP[loc] + UP[ind]) / 2 * exLW[p] + (bs->VP[loc] + VP[ind]) / 2 * eyLW[p] + (bs->WP[loc] + WP[ind]) / 2 * ezLW[p];
	// low order (upwind) and high order (Lax-Wendroff) flux
	double fluxLow = (Un > 0) ? Un * cA : Un * cB;
	double fluxHigh = Un * (cA + cB) / 2 - Un * Un * (cB - cA) / 2;

	// limiter
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

// conservative: each face contributes (F(A->B) - F(B->A))/2, the neighbour computes the same value with opposite sign
void BlockFreeSurfaceEulerParticlesFVM::computeFlux() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		fluxConcentration[ind] = 0;
		if (boundary[ind] || (flagMap[ind] != 1 && flagMap[ind] != 3)) {
			continue;
		}
		for (int p = 0; p < 6; p++) {
			int bl, loc;
			local2globalShift(indBlock, ind, bl, loc, exLW[p], eyLW[p], ezLW[p]);
			if (bl < 0 || bl >= nBlockMax || !blocks[bl] || blocks[bl]->boundary[loc]) {
				continue;
			}
			BlockFreeSurfaceEulerParticlesFVM* bs = (BlockFreeSurfaceEulerParticlesFVM*)blocks[bl];
			if (bs->flagMap[loc] != 1 && bs->flagMap[loc] != 3) {
				continue;
			}
			fluxConcentration[ind] += (faceFlux(ind, p) - bs->faceFlux(loc, oppLW[p])) / 2;
		}
	}
}

void BlockFreeSurfaceEulerParticlesFVM::applyFlux() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		concentration[ind] -= fluxConcentration[ind];
	}
}

void BlockFreeSurfaceEulerParticlesFVM::upwind() {
	computeFlux();
	applyFlux();
}

void BlockFreeSurfaceEulerParticlesFVM::computeNewParticleVelocityField() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (concentration[ind] > 1e-3) {
			double mass = concentration[ind] * massParticle;
			double volume = concentration[ind] * volumeParticle;
			double area = concentration[ind] * crossectionParticle;

			double FinterX = 0*dPPx[ind];
			double FinterY = 0*dPPy[ind];
			double FinterZ = 0*dPPz[ind];

			double Uslip = sqrt((UP[ind] - UF[ind]) * (UP[ind] - UF[ind]) + (VP[ind] - VF[ind]) * (VP[ind] - VF[ind]) + (WP[ind] - WF[ind]) * (WP[ind] - WF[ind]));
			double Cd_mult_Uslip = cd(Uslip);
			double coef = 0.5 * Cd_mult_Uslip * area;

			//UP[ind] = (mass * UP[ind] + coef * UF[ind] + G[0] * mass - 2 * dpX[ind] * volume + FinterX) / (mass + coef);
			//VP[ind] = (mass * VP[ind] + coef * VF[ind] + G[1] * mass - 2 * dpY[ind] * volume + FinterY) / (mass + coef);
			//WP[ind] = (mass * WP[ind] + coef * WF[ind] + G[2] * mass - 2 * dpZ[ind] * volume + FinterZ) / (mass + coef);
			
			UP[ind] = (mass * UP[ind] + coef * UF[ind] + G[0] * (mass - volume)) / (mass + coef);
			VP[ind] = (mass * VP[ind] + coef * VF[ind] + G[1] * (mass - volume)) / (mass + coef);
			WP[ind] = (mass * WP[ind] + coef * WF[ind] + G[2] * (mass - volume)) / (mass + coef);

			Gx[ind] = coef * (UP[ind] - UF[ind]);
			Gy[ind] = coef * (VP[ind] - VF[ind]);
			Gz[ind] = coef * (WP[ind] - WF[ind]);
		}
		else {
			UP[ind] = 0;
			VP[ind] = 0;
			WP[ind] = 0;
		}
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


void BlockFreeSurfaceEulerParticlesFVM::computeGradP() {
	int indBlockShift = 0;
	int indShift = 0;
	double* auxRho = new double[nf];

	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		dpX[ind] = 0;
		dpY[ind] = 0;
		dpZ[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indShift = ind + (ni * nj) * ez[dir] + ni * ey[dir] + ex[dir];
				bool dirExist = false;
				if (!boundary[indShift] && (flagMap[indShift] == 1 || flagMap[indShift] == 3)) {
					auxRho[dir] = Rho[indShift];
					dirExist = true;
				}
				if (!dirExist) {
					auxRho[dir] = Rho[ind];
				}
			}

			double px = 0;
			double py = 0;
			double pz = 0;
			for (int dir = 1; dir < nf; dir++) {
				px += w[dir] * ex[dir] * auxRho[dir];
				py += w[dir] * ey[dir] * auxRho[dir];
				pz += w[dir] * ez[dir] * auxRho[dir];
			}
			dpX[ind] = 1.5 * px / 3;
			dpY[ind] = 1.5 * py / 3;
			dpZ[ind] = 1.5 * pz / 3;
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		dpX[ind] = 0;
		dpY[ind] = 0;
		dpZ[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indBlockShift = nbBlock(ind, dir);
				int indShift = nbLocal(ind, dir);
				bool dirExist = false;
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
					if (!(blocks[indBlockShift]->boundary[indShift]) && (fm == 1 || fm == 3)) {
						auxRho[dir] = ((BlockFreeSurface*)blocks[indBlockShift])->Rho[indShift];
						dirExist = true;
					}
				}
				if (!dirExist) {
					auxRho[dir] = Rho[ind];
				}
			}

			double px = 0;
			double py = 0;
			double pz = 0;
			for (int dir = 1; dir < nf; dir++) {
				px += w[dir] * ex[dir] * auxRho[dir];
				py += w[dir] * ey[dir] * auxRho[dir];
				pz += w[dir] * ez[dir] * auxRho[dir];
			}
			dpX[ind] = 1.5 * px / 3;
			dpY[ind] = 1.5 * py / 3;
			dpZ[ind] = 1.5 * pz / 3;
		}
	}
}

void BlockFreeSurfaceEulerParticlesFVM::interpolateFluidVelocity() {
	for (int c = 0; c < nAll; c++) {
		int ind = indAll[c];
		// save old velocities
		UF[ind] = 0;
		VF[ind] = 0;
		WF[ind] = 0;

		if (!boundary[ind]) {
			if (flagMap[ind] == 1 || flagMap[ind] == 3) {
				double sum = 0;
				for (int i = -interpLevel; i < interpLevel + 1; i++) {
					for (int j = -interpLevel; j < interpLevel + 1; j++) {
						for (int k = -interpLevel; k < interpLevel + 1; k++) {
							int indBlockShift = 0;
							int indShift = 0;
							local2globalShift(indBlock, ind, indBlockShift, indShift, i, j, k);
							if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
								int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
								if (fm == 1 || fm == 3) {
									double weight = deltaValue[abs(i)] * deltaValue[abs(j)] * deltaValue[abs(k)];
									UF[ind] += weight * ((BlockFreeSurface*)blocks[indBlockShift])->U[indShift];
									VF[ind] += weight * ((BlockFreeSurface*)blocks[indBlockShift])->V[indShift];
									WF[ind] += weight * ((BlockFreeSurface*)blocks[indBlockShift])->W[indShift];
									sum += weight;
								}
							}
						}
					}
				}
				UF[ind] /= sum;
				VF[ind] /= sum;
				WF[ind] /= sum;
			}
		}
	}
}

void BlockFreeSurfaceEulerParticlesFVM::smoothForce() {
	for (int c = 0; c < nAll; c++) {
		int ind = indAll[c];
		if (!boundary[ind]) {
			if (flagMap[ind] == 1 || flagMap[ind] == 3) {
				double sum = 0;
				for (int i = -interpLevel; i < interpLevel + 1; i++) {
					for (int j = -interpLevel; j < interpLevel + 1; j++) {
						for (int k = -interpLevel; k < interpLevel + 1; k++) {
							int indBlockShift = 0;
							int indShift = 0;
							local2globalShift(indBlock, ind, indBlockShift, indShift, i, j, k);
							if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
								int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
								if (fm == 1 || fm == 3) {
									double weight = deltaValue[abs(i)] * deltaValue[abs(j)] * deltaValue[abs(k)];
									((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlockShift])->GxSmooth[indShift] += weight * Gx[ind];
									((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlockShift])->GySmooth[indShift] += weight * Gy[ind];
									((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlockShift])->GzSmooth[indShift] += weight * Gz[ind];
								}
							}
						}
					}
				}
			}
		}
	}
}

void BlockFreeSurfaceEulerParticlesFVM::computeParticlePressure() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		PP[ind] = 0*concentration[ind] / Cmax;
	}
}

void BlockFreeSurfaceEulerParticlesFVM::computeParticlePressureGradient() {
	int indBlockShift = 0;
	int indShift = 0;
	double* auxPP = new double[nf];

	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		dPPx[ind] = 0;
		dPPy[ind] = 0;
		dPPz[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indShift = ind + (ni * nj) * ez[dir] + ni * ey[dir] + ex[dir];
				bool dirExist = false;
				if (!boundary[indShift] && (flagMap[indShift] == 1 || flagMap[indShift] == 3)) {
					auxPP[dir] = PP[indShift];
					dirExist = true;
				}
				if (!dirExist) {
					auxPP[dir] = PP[ind];
				}
			}

			double px = 0;
			double py = 0;
			double pz = 0;
			for (int dir = 1; dir < nf; dir++) {
				px += w[dir] * ex[dir] * auxPP[dir];
				py += w[dir] * ey[dir] * auxPP[dir];
				pz += w[dir] * ez[dir] * auxPP[dir];
			}
			dPPx[ind] = 1.5 * px / 3;
			dPPy[ind] = 1.5 * py / 3;
			dPPz[ind] = 1.5 * pz / 3;
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		dPPx[ind] = 0;
		dPPy[ind] = 0;
		dPPz[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indBlockShift = nbBlock(ind, dir);
				int indShift = nbLocal(ind, dir);
				bool dirExist = false;
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
					if (!(blocks[indBlockShift]->boundary[indShift]) && (fm == 1 || fm == 3)) {
						auxPP[dir] = ((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlockShift])->PP[indShift];
						dirExist = true;
					}
				}
				if (!dirExist) {
					auxPP[dir] = PP[ind];
				}
			}

			double px = 0;
			double py = 0;
			double pz = 0;
			for (int dir = 1; dir < nf; dir++) {
				px += w[dir] * ex[dir] * auxPP[dir];
				py += w[dir] * ey[dir] * auxPP[dir];
				pz += w[dir] * ez[dir] * auxPP[dir];
			}
			dPPx[ind] = 1.5 * px / 3;
			dPPy[ind] = 1.5 * py / 3;
			dPPz[ind] = 1.5 * pz / 3;
		}
	}
}

double BlockFreeSurfaceEulerParticlesFVM::delta(double x) {
	x = abs(x);
	double y = 0;
	if (x <= 2) {
		y = (1 + cos(3.1415 * x / 2.0)) / 4.0;
	}
	return y;
}