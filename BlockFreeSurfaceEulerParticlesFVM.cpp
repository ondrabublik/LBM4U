
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

BlockFreeSurfaceEulerParticlesFVM::BlockFreeSurfaceEulerParticlesFVM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int* indAll, int* indIn, int* indBound, bool mrtModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, unordered_map<string, string> par)
	: BlockFreeSurface(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, coarseLevel, indAll, indIn, indBound, mrtModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {

	concentration = new double[ntot];
	fluxConcentration = new double[ntot];

	dpX = new double[ntot];
	dpY = new double[ntot];
	dpZ = new double[ntot];

	UP = new double[ntot];
	VP = new double[ntot];
	WP = new double[ntot];

	PP = new double[ntot]; // particle pressure
	dPPx = new double[ntot];
	dPPy = new double[ntot];
	dPPz = new double[ntot];

	UF = new double[ntot];
	VF = new double[ntot];
	WF = new double[ntot];

	Gx = new double[ntot];
	Gy = new double[ntot];
	Gz = new double[ntot];

	GxSmooth = new double[ntot];
	GySmooth = new double[ntot];
	GzSmooth = new double[ntot];

	fcP = new double[ntot * nf];

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
	oppLW = new int[6]{ 3, 0, 4, 1, 5, 2 };
}

BlockFreeSurfaceEulerParticlesFVM::~BlockFreeSurfaceEulerParticlesFVM() {
	delete[] fcP;
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
	catch (int ex) {
		cout << "Error in " << fileName << endl;
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

void BlockFreeSurfaceEulerParticlesFVM::collision() {
	double force[3];
	double feq[27];
	double feqNew[27];
	if (mrtModel) {
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (flagMap[ind] == 1 || flagMap[ind] == 3) {
				equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind]);
				force[0] = G[0] + GxSmooth[ind] * vof[ind];
				force[1] = G[1] + GySmooth[ind] * vof[ind];
				force[2] = G[2] + GzSmooth[ind] * vof[ind];
				equilibrium(feqNew, Rho[ind], dRho[ind], U[ind] + force[0], V[ind] + force[1], W[ind] + force[2]);
				for (int p = 0; p < nf; p++) {
					double sum = 0;
					int pind = nf * ind + p;
					for (int q = 0; q < nf; q++) {
						int qind = nf * ind + q;
						sum += T[nf * p + q] * (f[qind] - feq[q]);
					}
					double bodyForce = feqNew[p] - feq[p];
					fc[pind] = f[pind] - sum + bodyForce;
				}
			}
			Gx[ind] = 0;
			Gy[ind] = 0;
			Gz[ind] = 0;
			GxSmooth[ind] = 0;
			GySmooth[ind] = 0;
			GzSmooth[ind] = 0;
		}
	}
	else {
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (flagMap[ind] == 1 || flagMap[ind] == 3) {
				equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind]);
				force[0] = G[0] + GxSmooth[ind] * vof[ind];
				force[1] = G[1] + GySmooth[ind] * vof[ind];
				force[2] = G[2] + GzSmooth[ind] * vof[ind];
				equilibrium(feqNew, Rho[ind], dRho[ind], U[ind] + force[0], V[ind] + force[1], W[ind] + force[2]);
				for (int p = 0; p < nf; p++) {
					int pind = nf * ind + p;
					double bodyForce = feqNew[p] - feq[p];
					fc[pind] = f[pind] - 1.0 / tau * (f[pind] - feq[p]) + bodyForce;
				}
			}
			Gx[ind] = 0;
			Gy[ind] = 0;
			Gz[ind] = 0;
			GxSmooth[ind] = 0;
			GySmooth[ind] = 0;
			GzSmooth[ind] = 0;
		}
	}
}

void BlockFreeSurfaceEulerParticlesFVM::upwind() {
	double fluxLow;
	double fluxHigh;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		fluxConcentration[ind] = 0;
	}

	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 0; p < 6; p++) {
				int indShift = ind + (ni * nj) * ezLW[p] + ni * eyLW[p] + exLW[p];
				if ((flagMap[indShift] == 1 || flagMap[indShift] == 3) && !boundary[indShift]) {
					// compute low order flux
					double Un = (UP[indShift] + UP[ind]) / 2 * exLW[p] + (VP[indShift] + VP[ind]) / 2 * eyLW[p] + (WP[indShift] + WP[ind]) / 2 * ezLW[p];
					if (Un > 0) {
						fluxLow = Un * concentration[ind];
					}
					else {
						fluxLow = Un * concentration[indShift];
					}

					// compute high order flux
					fluxHigh = Un * (concentration[ind] + concentration[indShift]) / 2 - Un*Un*(concentration[indShift] - concentration[ind])/2;

					// compute limiter
					int pop = oppLW[p];
					int indShiftOpp = ind + (ni * nj) * ezLW[pop] + ni * eyLW[pop] + exLW[pop];
					double limiter = 0;
					if ((flagMap[indShiftOpp] == 1 || flagMap[indShiftOpp] == 3) && !boundary[indShiftOpp]) {
						if (concentration[indShift] != concentration[ind]) {
							double r = (concentration[ind] - concentration[indShiftOpp]) / (concentration[indShift] - concentration[ind]);
							limiter = (r + abs(r)) / (1 + abs(r));
						}
					}

					double flux = fluxLow - limiter * (fluxLow - fluxHigh);
					fluxConcentration[ind] += flux/2;
					fluxConcentration[indShift] -= flux/2;
				}
			}
		}
	}

	int indBlockShift = 0;
	int indShift = 0;
	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 0; p < 6; p++) {
				local2globalShift(indBlock, ind, indBlockShift, indShift, exLW[p], eyLW[p], ezLW[p]);
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					BlockFreeSurfaceEulerParticlesFVM* bs = ((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlockShift]);
					int flag = bs->flagMap[indShift];
					bool boundarys = bs->boundary[indShift];
					if ((flag == 1 || flag == 3) && !boundarys) {
						double Un = (bs->UP[indShift] + UP[ind]) / 2 * exLW[p] + (bs->VP[indShift] + VP[ind]) / 2 * eyLW[p] + (bs->WP[indShift] + WP[ind]) / 2 * ezLW[p];
						if (Un > 0) {
							fluxLow = Un * concentration[ind];
						}
						else {
							fluxLow = Un * bs->concentration[indShift];
						}
						fluxHigh = Un * (concentration[ind] + concentration[indShift]) / 2 - Un * Un * (concentration[indShift] - concentration[ind]) / 2;

						// compute limiter
						int indBlockShiftOpp = 0;
						int indShiftOpp = 0;
						int pop = oppLW[p];
						double limiter = 0;
						local2globalShift(indBlock, ind, indBlockShiftOpp, indShiftOpp, exLW[pop], eyLW[pop], ezLW[pop]);
						if (indBlockShiftOpp > -1 && indBlockShiftOpp < nBlockMax && blocks[indBlockShiftOpp]) {
							BlockFreeSurfaceEulerParticlesFVM* bsOpp = ((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlockShiftOpp]);
							int flag = bsOpp->flagMap[indShiftOpp];
							bool boundarys = bsOpp->boundary[indShiftOpp];
							if ((flag == 1 || flag == 3) && !boundarys) {
								if (bs->concentration[indShift] != concentration[ind]) {
									double r = (concentration[ind] - bsOpp->concentration[indShiftOpp]) / (bs->concentration[indShift] - concentration[ind]);
									limiter = (r + abs(r)) / (1 + abs(r));
								}
							}
						}

						double flux = fluxLow - limiter * (fluxLow - fluxHigh);
						fluxConcentration[ind] += flux/2;
						fluxConcentration[indShift] -= flux/2;
					}
				}
			}
		}
	}

	double Ctotal = 0;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		concentration[ind] -= fluxConcentration[ind];
		Ctotal += concentration[ind];
	}
	//cout << Ctotal << endl;
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
				int indBlockShift = bShift[nAll * dir + ind];
				int indShift = shift[nAll * dir + ind];
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
				int indBlockShift = bShift[nAll * dir + ind];
				int indShift = shift[nAll * dir + ind];
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