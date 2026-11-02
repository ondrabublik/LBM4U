
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 /*
  * File:   BlockFreeSurfaceEulerParticlesLBM.cpp
  * Author: obublik
  *
  * Created on 2. �nora 2018, 12:27
  */

#include "BlockFreeSurfaceEulerParticlesLBM.h"

using namespace std;

BlockFreeSurfaceEulerParticlesLBM::BlockFreeSurfaceEulerParticlesLBM(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurface(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par) {

	concentration = new double[ntot]();
	gc = new double[ntot * nf]();
	g = new double[ntot * nf]();

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
	interPartFCoef = getParDouble("particleInterForces", par, true, 0) / Cf;
	taus = 1.0 / 2 + 3 * diffusLB;
	double diffusLBoverPack = getParDouble("overpackDiffusivity", par, true, 0);
	tauOverPack = 3 * diffusLBoverPack;
	theta = 1.0 / 1.4;

	positiveDirs = new int[(nf - 1) / 2];
	int* aux = new int[nf];
	for (int p = 0; p < nf; p++) {
		aux[p] = 0;
	}
	int s = 0;
	for (int p = 1; p < nf; p++) {
		if (aux[p] != 1) {
			positiveDirs[s] = p;
			aux[p] = 1;
			aux[opp[p]] = 1;
			s++;
		}
	}

	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		tauDir[nf * ind + 0] = taus;
	}

	double pi = 3.14159265358979323846;
	crossectionParticle = pi * dParticle * dParticle / 4;
	volumeParticle = pi * dParticle * dParticle * dParticle / 6;
	massParticle = rhoParticle * volumeParticle;
	concentrationLimit = 0.74 / (4.0 / 3.0 * 3.1415 * dParticle * dParticle * dParticle / 8);
	interpLevel = 2;

	if (isFirst) {
		cout << "taus = " << taus << endl;
		cout << "maximal allowed value of concetration = " << concentrationLimit << endl;
	}
}

BlockFreeSurfaceEulerParticlesLBM::~BlockFreeSurfaceEulerParticlesLBM() {
	delete[] gc;
	delete[] g;
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
	delete[] deltaValue;
}

void BlockFreeSurfaceEulerParticlesLBM::readInitConcentration(int id, string simulPath) {
	string line;
	string fileName = simulPath + "init/concentration" + fileTag + ".txt";
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
			}
		}
	}
}

void BlockFreeSurfaceEulerParticlesLBM::propagationScalar() {
	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 0; p < nf; p++) {
				int indShift = ind + (ni * nj) * ez[p] + ni * ey[p] + ex[p];
				if (boundary[indShift] || flagMap[indShift] == 2) { // wall boundary
					g[nf * ind + opp[p]] = gc[nf * ind + p];
				}
				else { // standard propagation
					g[nf * ind + opp[p]] = gc[nf * indShift + opp[p]];
				}
			}
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 0; p < nf; p++) {
				int indLocalNeigh = nbLocal(ind, p);
				int indBlockNeigh = nbBlock(ind, p);
				bool noBlockWall = indBlockNeigh > -1 && indBlockNeigh < nBlockMax && blocks[indBlockNeigh];
				if (noBlockWall) {
					if (blocks[indBlockNeigh]->boundary[indLocalNeigh] || ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockNeigh])->flagMap[indLocalNeigh] == 2) { // wall boundary
						g[nf * ind + opp[p]] = gc[nf * ind + p];
					}
					else // standard propagation
						g[nf * ind + opp[p]] = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockNeigh])->gc[nf * (indLocalNeigh)+opp[p]];
				}
				else {  // wall boundary
					g[nf * ind + opp[p]] = gc[nf * ind + p];
				}
			}
		}
	}
}

void BlockFreeSurfaceEulerParticlesLBM::calculConcentration() {
	concentrationSum = 0.0;
	concentrationMax = 0.0;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		concentration[ind] = 0.0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {

			double auxsum = 1.0 /(3.0 * tauDir[nf * ind + 0]) + (theta - 1.0) / tauDir[nf * ind + 0];
			for (int p = 1; p < nf; p++) {
				auxsum += w[p]/tauDir[nf * ind + p];
			}
			for (int p = 0; p < nf; p++) {
				concentration[ind] += g[nf * ind + p] / tauDir[nf * ind + p] / auxsum;
			}
			concentrationSum += concentration[ind];
			if (concentration[ind] > concentrationMax) concentrationMax = concentration[ind];
		}
	}

	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		user_def[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 1; p < nf; p++) {
				if (tauDir[nf * ind + p] > user_def[ind]) {
					user_def[ind] = tauDir[nf * ind + p];
				}
			}
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

void BlockFreeSurfaceEulerParticlesLBM::computeNewParticleVelocityField(double* dist, double* gradDist) {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			if (concentration[ind] > 1e-4) {
				double mass = concentration[ind] * massParticle;
				double volume = concentration[ind] * volumeParticle;
				double area = concentration[ind] * crossectionParticle;

				double Uslip = sqrt((UP[ind] - UF[ind]) * (UP[ind] - UF[ind]) + (VP[ind] - VF[ind]) * (VP[ind] - VF[ind]) + (WP[ind] - WF[ind]) * (WP[ind] - WF[ind]));
				double Cd_mult_Uslip = cd(Uslip);
				double coef = 0.5 * Cd_mult_Uslip * area;

				double up = UP[ind], vp = VP[ind], wp = WP[ind];
				UP[ind] = (mass * (up - up * gUx[ind] - vp * gUy[ind] - wp * gUz[ind]) + coef * UF[ind] + G[0] * (mass - volume)) / (mass + coef);
				VP[ind] = (mass * (vp - up * gVx[ind] - vp * gVy[ind] - wp * gVz[ind]) + coef * VF[ind] + G[1] * (mass - volume)) / (mass + coef);
				WP[ind] = (mass * (wp - up * gWx[ind] - vp * gWy[ind] - wp * gWz[ind]) + coef * WF[ind] + G[2] * (mass - volume)) / (mass + coef);
					
				if (dist[ind] < 1) {
					UP[ind] = UF[ind];
					VP[ind] = VF[ind];
					WP[ind] = WF[ind];
				}

				Gx[ind] = coef * (UP[ind] - UF[ind]);
				Gy[ind] = coef * (VP[ind] - VF[ind]);
				Gz[ind] = coef * (WP[ind] - WF[ind]);
			}
			else {
				UP[ind] = UF[ind];
				VP[ind] = VF[ind];
				WP[ind] = WF[ind];
			}
		}
		else {
			UP[ind] = 0;
			VP[ind] = 0;
			WP[ind] = 0;
		}
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

void BlockFreeSurfaceEulerParticlesLBM::computeGradP() {
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


void BlockFreeSurfaceEulerParticlesLBM::computeGradU() {
	int indBlockShift = 0;
	int indShift = 0;
	double* auxU = new double[nf];
	double* auxV = new double[nf];
	double* auxW = new double[nf];

	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		gUx[ind] = 0;
		gUy[ind] = 0;
		gUz[ind] = 0;
		gVx[ind] = 0;
		gVy[ind] = 0;
		gVz[ind] = 0;
		gWx[ind] = 0;
		gWy[ind] = 0;
		gWz[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indShift = ind + (ni * nj) * ez[dir] + ni * ey[dir] + ex[dir];
				bool dirExist = false;
				if (!boundary[indShift] && (flagMap[indShift] == 1 || flagMap[indShift] == 3)) {
					auxU[dir] = UP[indShift];
					auxV[dir] = VP[indShift];
					auxW[dir] = WP[indShift];
					dirExist = true;
				}
				if (!dirExist) {
					auxU[dir] = UP[ind];
					auxV[dir] = VP[ind];
					auxW[dir] = WP[ind];
				}
			}

			for (int dir = 1; dir < nf; dir++) {
				gUx[ind] += 0.5 * w[dir] * ex[dir] * auxU[dir];
				gUy[ind] += 0.5 * w[dir] * ey[dir] * auxU[dir];
				gUz[ind] += 0.5 * w[dir] * ez[dir] * auxU[dir];

				gVx[ind] += 0.5 * w[dir] * ex[dir] * auxV[dir];
				gVy[ind] += 0.5 * w[dir] * ey[dir] * auxV[dir];
				gVz[ind] += 0.5 * w[dir] * ez[dir] * auxV[dir];

				gWx[ind] += 0.5 * w[dir] * ex[dir] * auxW[dir];
				gWy[ind] += 0.5 * w[dir] * ey[dir] * auxW[dir];
				gWz[ind] += 0.5 * w[dir] * ez[dir] * auxW[dir];
				
			}
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		gUx[ind] = 0;
		gUy[ind] = 0;
		gUz[ind] = 0;
		gVx[ind] = 0;
		gVy[ind] = 0;
		gVz[ind] = 0;
		gWx[ind] = 0;
		gWy[ind] = 0;
		gWz[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indBlockShift = nbBlock(ind, dir);
				int indShift = nbLocal(ind, dir);
				bool dirExist = false;
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
					if (!(blocks[indBlockShift]->boundary[indShift]) && (fm == 1 || fm == 3)) {
						auxU[dir] = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->UP[indShift];
						auxV[dir] = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->VP[indShift];
						auxW[dir] = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->WP[indShift];
						dirExist = true;
					}
				}
				if (!dirExist) {
					auxU[dir] = UP[ind];
					auxV[dir] = VP[ind];
					auxW[dir] = WP[ind];
				}
			}

			for (int dir = 1; dir < nf; dir++) {
				gUx[ind] += 0.5 * w[dir] * ex[dir] * auxU[dir];
				gUy[ind] += 0.5 * w[dir] * ey[dir] * auxU[dir];
				gUz[ind] += 0.5 * w[dir] * ez[dir] * auxU[dir];

				gVx[ind] += 0.5 * w[dir] * ex[dir] * auxV[dir];
				gVy[ind] += 0.5 * w[dir] * ey[dir] * auxV[dir];
				gVz[ind] += 0.5 * w[dir] * ez[dir] * auxV[dir];

				gWx[ind] += 0.5 * w[dir] * ex[dir] * auxW[dir];
				gWy[ind] += 0.5 * w[dir] * ey[dir] * auxW[dir];
				gWz[ind] += 0.5 * w[dir] * ez[dir] * auxW[dir];
			}
		}
	}
}

void BlockFreeSurfaceEulerParticlesLBM::computeGradC() {
	int indBlockShift = 0;
	int indShift = 0;
	double* auxC = new double[nf];

	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		dcX[ind] = 0;
		dcY[ind] = 0;
		dcZ[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indShift = ind + (ni * nj) * ez[dir] + ni * ey[dir] + ex[dir];
				bool dirExist = false;
				if (!boundary[indShift] && (flagMap[indShift] == 1 || flagMap[indShift] == 3)) {
					auxC[dir] = concentration[indShift];
					dirExist = true;
				}
				if (!dirExist) {
					auxC[dir] = concentration[ind];
				}
			}

			double cx = 0;
			double cy = 0;
			double cz = 0;
			for (int dir = 1; dir < nf; dir++) {
				cx += w[dir] * ex[dir] * auxC[dir];
				cy += w[dir] * ey[dir] * auxC[dir];
				cz += w[dir] * ez[dir] * auxC[dir];
			}
			dcX[ind] = 1.5 * cx / 3;
			dcY[ind] = 1.5 * cy / 3;
			dcZ[ind] = 1.5 * cz / 3;
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		dcX[ind] = 0;
		dcY[ind] = 0;
		dcZ[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int dir = 1; dir < nf; dir++) {
				int indBlockShift = nbBlock(ind, dir);
				int indShift = nbLocal(ind, dir);
				bool dirExist = false;
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
					if (!(blocks[indBlockShift]->boundary[indShift]) && (fm == 1 || fm == 3)) {
						auxC[dir] = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->concentration[indShift];
						dirExist = true;
					}
				}
				if (!dirExist) {
					auxC[dir] = concentration[ind];
				}
			}

			double cx = 0;
			double cy = 0;
			double cz = 0;
			for (int dir = 1; dir < nf; dir++) {
				cx += w[dir] * ex[dir] * auxC[dir];
				cy += w[dir] * ey[dir] * auxC[dir];
				cz += w[dir] * ez[dir] * auxC[dir];
			}
			dcX[ind] = 1.5 * cx / 3;
			dcY[ind] = 1.5 * cy / 3;
			dcZ[ind] = 1.5 * cz / 3;
		}
	}
}

void BlockFreeSurfaceEulerParticlesLBM::computeGradVOFp() {
	int indBlockShift = 0;
	int indShift = 0;
	double* auxVOFp = new double[nf];

	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		gVOFpX[ind] = 0;
		gVOFpY[ind] = 0;
		gVOFpZ[ind] = 0;
		for (int dir = 1; dir < nf; dir++) {
			int indShift = ind + (ni * nj) * ez[dir] + ni * ey[dir] + ex[dir];
			bool dirExist = false;
			if (!boundary[indShift]) {
				auxVOFp[dir] = concentration[indShift];
				dirExist = true;
			}
			if (!dirExist) {
				auxVOFp[dir] = concentration[ind];
			}
		}
		double wtot = 0.0;
		for (int dir = 1; dir < nf; dir++) {
			gVOFpX[ind] += w[dir] * ex[dir] * auxVOFp[dir];
			gVOFpY[ind] += w[dir] * ey[dir] * auxVOFp[dir];
			gVOFpZ[ind] += w[dir] * ez[dir] * auxVOFp[dir];
			wtot += w[dir];
		}
		gVOFpX[ind] /= (concentrationLimit * wtot);
		gVOFpY[ind] /= (concentrationLimit * wtot);
		gVOFpZ[ind] /= (concentrationLimit * wtot);
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		gVOFpX[ind] = 0;
		gVOFpY[ind] = 0;
		gVOFpZ[ind] = 0;
			for (int dir = 1; dir < nf; dir++) {
				int indBlockShift = nbBlock(ind, dir);
				int indShift = nbLocal(ind, dir);
				bool dirExist = false;
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (!(blocks[indBlockShift]->boundary[indShift])) {
						auxVOFp[dir] = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->concentration[indShift];
						dirExist = true;
					}
				}
				if (!dirExist) {
					auxVOFp[dir] = concentration[ind];
				}
			}
			double wtot = 0.0;
			for (int dir = 1; dir < nf; dir++) {
				gVOFpX[ind] += w[dir] * ex[dir] * auxVOFp[dir];
				gVOFpY[ind] += w[dir] * ey[dir] * auxVOFp[dir];
				gVOFpZ[ind] += w[dir] * ez[dir] * auxVOFp[dir];
				wtot += w[dir];
			}
			gVOFpX[ind] /= (concentrationLimit * wtot);
			gVOFpY[ind] /= (concentrationLimit * wtot);
			gVOFpZ[ind] /= (concentrationLimit * wtot);
	}
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

void BlockFreeSurfaceEulerParticlesLBM::interpolateFluidVelocity() {
	for (int c = 0; c < nAll; c++) {
		int ind = indAll[c];
		// save old velocities
		UF[ind] = 0;
		VF[ind] = 0;
		WF[ind] = 0;

		if (!boundary[ind]) {
			if (flagMap[ind] == 1) {
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

void BlockFreeSurfaceEulerParticlesLBM::smoothForce() {
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
									// ???????????????????????????????? nechapu GxSmooth
									//((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->GxSmooth[indShift] += weight * Gx[ind];
									//((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->GySmooth[indShift] += weight * Gy[ind];
									//((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->GzSmooth[indShift] += weight * Gz[ind];

									GxSmooth[ind] += weight * ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->Gx[indShift];
									GySmooth[ind] += weight * ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->Gy[indShift];
									GzSmooth[ind] += weight * ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShift])->Gz[indShift];

									sum += weight;
								}
							}
						}
					}
				}
				GxSmooth[ind] /= sum;
				GySmooth[ind] /= sum;
				GzSmooth[ind] /= sum;
			}
		}
	}
}

double BlockFreeSurfaceEulerParticlesLBM::delta(double x) {
	x = abs(x);
	double y = 0;
	if (x <= 2) {
		y = (1 + cos(3.1415 * x / 2.0)) / 4.0;
	}
	return y;
}

void BlockFreeSurfaceEulerParticlesLBM::setTau() {
	int nPos = (nf - 1) / 2;
	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 0; p < nPos; p++) {
				int pFor = positiveDirs[p];
				int pBack = opp[pFor];
				int indFor = ind + (ni * nj) * ez[pFor] + ni * ey[pFor] + ex[pFor];
				int indBack = ind + (ni * nj) * ez[pBack] + ni * ey[pBack] + ex[pBack];
				double cm = concentration[ind];
				double cf = concentration[ind];
				double cb = concentration[ind];
				if (!boundary[indFor]) {
					cf = concentration[indFor];
				}
				if (!boundary[indBack]) {
					cb = concentration[indBack];
				}
				
				double tauLoc = computeTau(cf, cm, cb);
				tauDir[nf * ind + pFor] = tauLoc;
				tauDir[nf * ind + pBack] = tauLoc;
			}
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 0; p < nPos; p++) {
				int pFor = positiveDirs[p];
				int pBack = opp[pFor];

				int indBlockShiftFor = nbBlock(ind, pFor);
				int indShiftFor = nbLocal(ind, pFor);

				int indBlockShiftBack = nbBlock(ind, pBack);
				int indShiftBack = nbLocal(ind, pBack);

				double cm = concentration[ind];
				double cf = concentration[ind];
				double cb = concentration[ind];
				if (indBlockShiftFor > -1 && indBlockShiftFor < nBlockMax && blocks[indBlockShiftFor] && !(blocks[indBlockShiftFor]->boundary[indShiftFor])) {
					cf = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShiftFor])->concentration[indShiftFor];
				}
				if (indBlockShiftBack > -1 && indBlockShiftBack < nBlockMax && blocks[indBlockShiftBack] && !(blocks[indBlockShiftBack]->boundary[indShiftBack])) {
					cb = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlockShiftBack])->concentration[indShiftBack];
				}
	
				double tauLoc = computeTau(cf, cm, cb);
				tauDir[nf * ind + pFor] = tauLoc;
				tauDir[nf * ind + pBack] = tauLoc;	
			}
		}
	}
}

double BlockFreeSurfaceEulerParticlesLBM::computeTau(double cf, double cm, double cb) {
	double d2 = cf - 2 * cm + cb;
	//double dVasek = (cf + cb) / 2 - cm;
	double overpack = max((cm - concentrationLimit) / concentrationLimit, 0.0);
	return taus + abs(d2)/1000 + overpack*tauOverPack;
}