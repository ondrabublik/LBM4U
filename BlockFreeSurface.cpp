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

#include "BlockFreeSurface.h"

using namespace std;

BlockFreeSurface::BlockFreeSurface(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int* indAll, int* indIn, int* indBound, bool mrtModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, unordered_map<string, string> par) {
	isBlockEmpty = false;
	this->blocks = blocks;
	this->nBlock = nBlock;
	this->nBlockMax = nBlockMax;
	this->blockMap = blockMap;
	this->nf = nf;
	this->mrtModel = mrtModel;
	this->viscosityModel = viscosityModel;
	this->ex = ex;
	this->ey = ey;
	this->ez = ez;
	this->w = w;
	this->opp = opp;
	this->nBi = nBi;
	this->nBj = nBj;
	this->nBk = nBk;
	this->indBlock = indBlock;
	this->coarseLevel = coarseLevel;
	this->indIn = indIn;
	this->indBound = indBound;
	this->indAll = indAll;
	this->ni = coarseLevel * ni;
	this->nj = coarseLevel * nj;
	this->nk = coarseLevel * nk;
	try {
		xPeriodicity = (stoi(par["xPeriodicity"]) ? true : false);
		yPeriodicity = (stoi(par["yPeriodicity"]) ? true : false);
		zPeriodicity = (stoi(par["zPeriodicity"]) ? true : false);
	}
	catch (const std::exception & e) {
		xPeriodicity = false;
		yPeriodicity = false;
		zPeriodicity = false;
		cout << "Periodic boundary not set!" << endl;
	}
	nAll = ni * nj * nk;
	nIn = (ni - 2) * (nj - 2) * (nk - 2);
	nBound = nAll - nIn;
	massErr = 0.0;
	Cu = getParDouble("uMax", par, true, 0) / 0.1;
	Ch = getParDouble("H", par, true, 0) / (this->ni);
	Ct = Ch / Cu;
	Crho = getParDouble("density", par, true, 0);
	Cp = Crho / 3 * Cu * Cu;
	Cst = Crho * Cu * Cu * Ch; // Cp*Ch
	Cvis = Ch * Ch / Ct;
	vLB = getParDouble("vKin", par, true, 0) / Cvis;
	Cf = Ch / (Ct * Ct);
	this->tau = 1.0 / 2 + 3 * vLB;
	surfTension = getParDouble("surfaceTension", par, false, 0) / Cst;
	G = new double[3]{
		getParDouble("gx",par, false, 0) / Cf, getParDouble("gy",par, false, 0) / Cf, getParDouble("gz",par, false, 0) / Cf
	}; //gravity
	G0 = new double[3]{
		getParDouble("gx",par, false, 0) / Cf, getParDouble("gy",par, false, 0) / Cf, getParDouble("gz",par, false, 0) / Cf
	}; //gravity
	velWall = new double[3]{
		getParDouble("uwall",par, false, 0) / Cu, getParDouble("vwall",par, false, 0) / Cu, getParDouble("wwall",par, false, 0) / Cu
	}; //velocity at the wall of cavity

	rho0 = 0;   //0 == normal, 1 == special treatment for higher precision

	try {
		wetAngleEq = stod(par["wetAngleEq"]) * 3.14159 / 180.0; // wetting angle - Equlibrium [rad] 
		wetAngleAdv = stod(par["wetAngleAdv"]) * 3.14159 / 180.0; // wetting angle - Advancing  [rad]
		wetAngleRec = stod(par["wetAngleRec"]) * 3.14159 / 180.0; // wetting angle - Receding   [rad]  
	}
	catch (const std::exception & e) {
		cout << "Contact angle not set!" << endl;
		wetAngleEq = 0;
		wetAngleAdv = 0;
		wetAngleRec = 0;
	}


	if (isFirst) {
		cout << "tau = " << tau << endl;
	}

	ntot = ni * nj * nk;
	boundary = new bool[ntot];
	fc = new double[ntot * nf];
	fcOld = new double[ntot * nf];
	f = new double[ntot * nf];
	Rho = new double[ntot];
	dRho = new double[ntot];
	U = new double[ntot];
	V = new double[ntot];
	W = new double[ntot];
	mass = new double[ntot];
	vof = new double[ntot];
	vofS = new double[ntot];
	vofS1 = new double[ntot];
	normal = new double[ntot * 3];
	curvature = new double[ntot];
	flagMap = new int[ntot];
	filled = new int[ntot];
	emptied = new int[ntot];
	airToFs = new int[ntot];
	fluidToFs = new int[ntot];
	artefact1 = new int[ntot];
	artefact2 = new int[ntot];
	surfaceLevels = new int[ntot];
	wallDistance = new double[ntot];
	wallDistanceGradient = new double[3 * ntot];
	wallDistancePolynome = new double[8 * ntot];

	if (mrtModel) {
		defineMRTmatrix();
	}

	if (viscosityModel > 0) {
		viscosity = new double[ntot];
		shearRate = new double[ntot];
		try {
			fluidModelParam = getParDoubleArray("fluidModelParameters", par);
			if (isFirst) {
				if (viscosityModel > 0) {
					switch (viscosityModel) {
					case 1: // power-law
						cout << "Power-Law model: v_kinematic = b[0]*gamma^(b[1]-1)" << endl;
						cout << "b[0] = " << fluidModelParam[0] << endl;
						cout << "b[1] = " << fluidModelParam[1] << endl;
						break;
					case 2: // bingham
						cout << "Bingham model: v_kinematic = b[0]+b[1]*gamma^(b[2]-1)" << endl;
						cout << "b[0] = " << fluidModelParam[0] << endl;
						cout << "b[1] = " << fluidModelParam[1] << endl;
						cout << "b[2] = " << fluidModelParam[2] << endl;
						break;
					case 3: // cross
						cout << "Cross model: v_kinematic = b[0]+(b[1]-b[0])/(1+b[2]*gamma^b[3])" << endl;
						cout << "b[0] = " << fluidModelParam[0] << endl;
						cout << "b[1] = " << fluidModelParam[1] << endl;
						cout << "b[2] = " << fluidModelParam[2] << endl;
						cout << "b[3] = " << fluidModelParam[3] << endl;
						break;
					case 4: // Poly
						cout << "Polynomial model: v_kinematic = b[0]+b[1]*gamma+b[2]*gamma^2+b[3]*gamma^3)" << endl;
						cout << "b[0] = " << fluidModelParam[0] << endl;
						cout << "b[1] = " << fluidModelParam[1] << endl;
						cout << "b[2] = " << fluidModelParam[2] << endl;
						cout << "b[3] = " << fluidModelParam[3] << endl;
						break;
					case 5: // Modified power-Law
						cout << "Modified power-Law: v_kinematic = b[0]*gamma^(b[1]-1 + b[3]*log(gamma))" << endl;
						cout << "b[0] = " << fluidModelParam[0] << endl;
						cout << "b[1] = " << fluidModelParam[1] << endl;
						cout << "b[2] = " << fluidModelParam[2] << endl;
						break;
					}
					cout << "tau(gamma = 0) = " << 1.0 / 2 + 3 * nonNewtonianViscosityModel(0) / Cvis << endl;
					cout << "tau(gamma = 0.1) = " << 1.0 / 2 + 3 * nonNewtonianViscosityModel(0.1) / Cvis << endl;
					cout << "tau(gamma = 1) = " << 1.0 / 2 + 3 * nonNewtonianViscosityModel(1) / Cvis << endl;
					cout << "tau(gamma = 10) = " << 1.0 / 2 + 3 * nonNewtonianViscosityModel(10) / Cvis << endl;
					cout << "tau(gamma = 1e4) = " << 1.0 / 2 + 3 * nonNewtonianViscosityModel(1e4) / Cvis << endl;
				}
				else {
					cout << "tau = " << 1.0 / 2 + 3 * vLB / Cvis << endl;
				}
			}
		}
		catch (const std::exception & e) {
			cout << "Please specify non-Newtonian model constants!" << endl;
			exit(EXIT_FAILURE);
		}
	}
}

BlockFreeSurface::~BlockFreeSurface() {
	delete[] G;
	delete[] velWall;
	delete[] boundary;
	delete[] fc;
	delete[] f;
	delete[] Rho;
	delete[] dRho;
	delete[] U;
	delete[] V;
	delete[] W;
	delete[] flagMap;
	delete[] mass;
	delete[] vof;
	delete[] vofS;
	delete[] vofS1;
	delete[] normal;
	delete[] curvature;
	delete[] filled;
	delete[] emptied;
	delete[] airToFs;
	delete[] fluidToFs;
	delete[] artefact1;
	delete[] artefact2;
	delete[] surfaceLevels;
	delete[] wallDistance;
	delete[] wallDistanceGradient;
}

void BlockFreeSurface::collision() {
	double force[3];
	double feq[27];
	double feqNew[27];
	double fAux[27];
	if (mrtModel) {
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (viscosityModel > 0) {
				shearRate[ind] = 0;
				viscosity[ind] = 0;
			}
			if (flagMap[ind] == 1 || flagMap[ind] == 3) {
				equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind]);
				force[0] = G[0];
				force[1] = G[1];
				force[2] = G[2];
				equilibrium(feqNew, Rho[ind], dRho[ind], U[ind] + force[0], V[ind] + force[1], W[ind] + force[2]);
				if (viscosityModel > 0) {
					shearRate[ind] = computeShearRate(ind) / Ct;
					viscosity[ind] = nonNewtonianViscosityModel(shearRate[ind]);
					double tauLoc = 0.5 + 3 * viscosity[ind] / Cvis;
					if (tauLoc > 1.5) {
						tauLoc = 1.5;
					}
					if (tauLoc < 0.5 + 1e-4) {
						tauLoc = 0.5 + 1e-4;
					}
					/*for (int p = 0; p < nf; p++) {
						double sum = 0;
						for (int q = 0; q < nf; q++) {
							int qind = nf * ind + q;
							sum += M[nf * p + q] * (f[qind] - feq[q]);
						}
						fAux[p] = sum;
					}
					setMRTdiagonal(tauLoc);
					for (int p = 0; p < nf; p++) {
						fAux[p] *= S[p];
					}
					for (int p = 0; p < nf; p++) {
						int pind = nf * ind + p;
						double sum = 0;
						for (int q = 0; q < nf; q++) {
							sum += iM[nf * p + q] * fAux[q];
						}
						double bodyForce = feqNew[p] - feq[p];
						fc[pind] = f[pind] - sum + bodyForce;
					}*/

					double omega = 1 / tauLoc;
					for (int p = 0; p < nf; p++) {
						double sum = 0;
						int pind = nf * ind + p;
						for (int q = 0; q < nf; q++) {
							int qind = nf * ind + q;
							sum += (T0[nf * p + q] + omega * T1[nf * p + q]) * (f[qind] - feq[q]);
						}
						double bodyForce = feqNew[p] - feq[p];
						fc[pind] = f[pind] - sum + bodyForce;
					}
				}
				else {
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
			}
		}
	}
	else {
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (flagMap[ind] == 1 || flagMap[ind] == 3) {
				equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind]);
				force[0] = G[0];
				force[1] = G[1];
				force[2] = G[2];
				equilibrium(feqNew, Rho[ind], dRho[ind], U[ind] + force[0], V[ind] + force[1], W[ind] + force[2]);
				for (int p = 0; p < nf; p++) {
					int pind = nf * ind + p;
					double bodyForce = feqNew[p] - feq[p];
					fc[pind] = f[pind] - 1.0 / tau * (f[pind] - feq[p]) + bodyForce;
				}
			}
		}
	}
}

void BlockFreeSurface::propagation() {
	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 0; p < nf; p++) {
				int indShift = ind + (ni * nj) * ez[p] + ni * ey[p] + ex[p];
				double weight = normal[3 * ind] * ex[p] + normal[3 * ind + 1] * ey[p] + normal[3 * ind + 2] * ez[p];
				if (flagMap[indShift] == 2 || (flagMap[ind] == 3 && boundary[indShift] && weight > 0)) {
					double u2 = U[ind] * U[ind] + V[ind] * V[ind] + W[ind] * W[ind];
					double x = ex[p] * U[ind] + ey[p] * V[ind] + ez[p] * W[ind];
					double xOpp = ex[opp[p]] * U[ind] + ey[opp[p]] * V[ind] + ez[opp[p]] * W[ind];
					double RhoS = 1.0 - 6.0 * surfTension * curvature[ind];
					double feq = w[p] * (RhoS * (3.0 * x + 9.0 / 2.0 * x * x - 3.0 / 2.0 * u2) + (RhoS - rho0));
					double feqOpp = w[opp[p]] * (RhoS * (3.0 * xOpp + 9.0 / 2.0 * xOpp * xOpp - 3.0 / 2.0 * u2) + (RhoS - rho0));
					f[nf * ind + opp[p]] = feq + feqOpp - fc[nf * ind + p];
				}
				else if (boundary[indShift]) {
					f[nf * ind + opp[p]] = fc[nf * ind + p];
				}
				else { // standart propagation
					f[nf * ind + opp[p]] = fc[nf * indShift + opp[p]];
				}
			}
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 0; p < nf; p++) {
				int indLocalNeigh = shift[nAll * p + ind];
				int indBlockNeigh = bShift[nAll * p + ind];
				bool noBlockWall = indBlockNeigh > -1 && indBlockNeigh < nBlockMax && blocks[indBlockNeigh];
				if (noBlockWall) {
					double weight = normal[3 * ind + 0] * ex[p] + normal[3 * ind + 1] * ey[p] + normal[3 * ind + 2] * ez[p];
					if (((BlockFreeSurface*)blocks[indBlockNeigh])->flagMap[indLocalNeigh] == 2 || (flagMap[ind] == 3 && blocks[indBlockNeigh]->boundary[indLocalNeigh] && weight > 0)) {
						double u2 = U[ind] * U[ind] + V[ind] * V[ind] + W[ind] * W[ind];
						double x = ex[p] * U[ind] + ey[p] * V[ind] + ez[p] * W[ind];
						double xOpp = ex[opp[p]] * U[ind] + ey[opp[p]] * V[ind] + ez[opp[p]] * W[ind];
						double RhoS = 1.0 - 6.0 * surfTension * curvature[ind];
						double feq = w[p] * (RhoS * (3.0 * x + 9.0 / 2.0 * x * x - 3.0 / 2.0 * u2) + (RhoS - rho0));
						double feqOpp = w[opp[p]] * (RhoS * (3.0 * xOpp + 9.0 / 2.0 * xOpp * xOpp - 3.0 / 2.0 * u2) + (RhoS - rho0));
						f[nf * ind + opp[p]] = feq + feqOpp - fc[nf * ind + p];
					}
					else if (blocks[indBlockNeigh]->boundary[indLocalNeigh]) {
						f[nf * ind + opp[p]] = fc[nf * ind + p]; // bounce back
					}
					else // standart propagation
						f[nf * ind + opp[p]] = blocks[indBlockNeigh]->fc[nf * (indLocalNeigh)+opp[p]];
				}
				else {
					double weight = normal[3 * ind + 0] * ex[p] + normal[3 * ind + 1] * ey[p] + normal[3 * ind + 2] * ez[p];
					if (flagMap[ind] == 3 && weight > 0) {
						double u2 = U[ind] * U[ind] + V[ind] * V[ind] + W[ind] * W[ind];
						double x = ex[p] * U[ind] + ey[p] * V[ind] + ez[p] * W[ind];
						double xOpp = ex[opp[p]] * U[ind] + ey[opp[p]] * V[ind] + ez[opp[p]] * W[ind];
						double RhoS = 1.0 - 6.0 * surfTension * curvature[ind];
						double feq = w[p] * (RhoS * (3.0 * x + 9.0 / 2.0 * x * x - 3.0 / 2.0 * u2) + (RhoS - rho0));
						double feqOpp = w[opp[p]] * (RhoS * (3.0 * xOpp + 9.0 / 2.0 * xOpp * xOpp - 3.0 / 2.0 * u2) + (RhoS - rho0));
						f[nf * ind + opp[p]] = feq + feqOpp - fc[nf * ind + p];
					}
					else { // bounce back

						// Pro ulohu proudeni v kavite
						int indBlock = bShift[ind];
						int kB = indBlock / (nBi * nBj);
						int k = ind / (ni * nj);
						if (kB == nBk-1 && (k + ez[p] == nk)) {
							f[nf * ind + opp[p]] = fc[nf * ind + p] - 6.0 * Rho[ind] * w[p] * (velWall[0] * ex[p] + velWall[1] * ey[p] + velWall[2] * ez[p]);
						}
						else {
							f[nf * ind + opp[p]] = fc[nf * ind + p];
						}

						//f[nf * ind + opp[p]] = fc[nf * ind + p];
					}
				}
			}
		}
	}
}

void BlockFreeSurface::makro() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		Rho[ind] = 0;
		dRho[ind] = 0;
		U[ind] = 0;
		V[ind] = 0;
		W[ind] = 0;
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			for (int p = 0; p < nf; p++) {
				double fVal = f[nf * ind + p];
				dRho[ind] += fVal;
				U[ind] += ex[p] * fVal;
				V[ind] += ey[p] * fVal;
				W[ind] += ez[p] * fVal;
			}
			Rho[ind] = dRho[ind] + rho0;
			U[ind] /= Rho[ind];
			V[ind] /= Rho[ind];
			W[ind] /= Rho[ind];
		}
	}
}

void BlockFreeSurface::initZeroCondition() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		flagMap[ind] = 4;
		Rho[ind] = 1.0;
		dRho[ind] = Rho[ind] - rho0;
		U[ind] = 0;
		V[ind] = 0;
		W[ind] = 0;
		mass[ind] = 0;
		vof[ind] = 0;
		curvature[ind] = 0;
	}
}

void BlockFreeSurface::readInitCondition(int id, string simulPath) {
	string line;
	string fileName = simulPath + "init/blockInit" + std::to_string(id) + ".txt";
	try {
		ifstream* myfile = new ifstream(fileName);
		if (myfile->is_open()) {
			for (int i = 0; i < nAll; i++) {
				int ind = indAll[i];
				getline(*myfile, line);
				istringstream iss(line);
				vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
				Rho[ind] = stod(tokens[0]);
				U[ind] = stod(tokens[1]) / Cu;
				V[ind] = stod(tokens[2]) / Cu;
				W[ind] = stod(tokens[3]) / Cu;
				dRho[ind] = Rho[ind] - rho0;
				if (!boundary[ind]) {
					flagMap[ind] = stoi(tokens[4]);
				}
			}
		}
		myfile->close();
		delete myfile;
	}
	catch (int ex) {
		cout << "Error in " << fileName << endl;
	}
}

void BlockFreeSurface::writeOutCondition(int id, string simulPath) {
	string line;
	string fileName = simulPath + "out/blockInit" + std::to_string(id) + ".txt";
	try {
		ofstream myfile;
		myfile.open(fileName);
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			myfile << Rho[ind] << " " << U[ind] * Cu << " " << V[ind] * Cu << " " << W[ind] * Cu << " " << flagMap[ind] << "\n";
		}
		myfile.close();
	}
	catch (int ex) {
		cout << "Error in " << fileName << endl;
	}
}


void BlockFreeSurface::computeInitDistribution() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		for (int p = 0; p < nf; p++) {
			double x = ex[p] * U[ind] + ey[p] * V[ind] + ez[p] * W[ind];
			double u2 = U[ind] * U[ind] + V[ind] * V[ind] + W[ind] * W[ind];
			fc[nf * ind + p] = w[p] * (Rho[ind] * (3.0 * x + 9.0 / 2.0 * x * x - 3.0 / 2.0 * u2) + dRho[ind]);
		}
	}
}

void BlockFreeSurface::initFlagMap() {
	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		if (flagMap[ind] == 1) {
			for (int dir = 1; dir < nf; dir++) {
				int indShift = ind + (ni * nj) * ez[dir] + ni * ey[dir] + ex[dir];
				if (flagMap[indShift] == 2) {
					flagMap[ind] = 3;
				}
			}
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		if (flagMap[ind] == 1) {
			for (int dir = 1; dir < nf; dir++) {
				int indBlockShift = bShift[nAll * dir + ind];
				int indShift = shift[nAll * dir + ind];
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 2) {
						flagMap[ind] = 3;
					}
				}
			}
		}
	}
}

void BlockFreeSurface::iniVofAndMass() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3) {
			vof[ind] = 0.5;
			mass[ind] = Rho[ind] * vof[ind];
		}
		else if (flagMap[ind] == 1) {
			vof[ind] = 1.0;
			mass[ind] = Rho[ind];
		}
	}
}

void BlockFreeSurface::massExchange() {
	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		if (flagMap[ind] == 3) {
			for (int p = 0; p < nf; p++) {
				int pind = nf * ind + p; // index pro hodnotu distribuce p v bode ijk
				int indShift = ind + (ni * nj) * ez[p] + ni * ey[p] + ex[p];
				if (flagMap[indShift] == 3) {
					mass[ind] -= (fc[pind] - fc[nf * indShift + opp[p]]) * (vof[ind] + vof[indShift]) / 2;
				}
				else if (flagMap[indShift] == 1) {
					mass[ind] -= (fc[pind] - fc[nf * indShift + opp[p]]);
				}
			}
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		if (flagMap[ind] == 3) {
			for (int p = 0; p < nf; p++) {
				int pind = nf * ind + p; // index pro hodnotu distribuce p v bode ijk
				int indBlockShift = bShift[nAll * p + ind];
				int indShift = shift[nAll * p + ind];
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 3) {
						mass[ind] -= (fc[pind] - blocks[indBlockShift]->fc[nf * (indShift)+opp[p]]) * (vof[ind] + ((BlockFreeSurface*)blocks[indBlockShift])->vof[indShift]) / 2;
					}
					else if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 1) {
						mass[ind] -= (fc[pind] - blocks[indBlockShift]->fc[nf * (indShift)+opp[p]]);
					}
				}
			}
		}
	}
}

void BlockFreeSurface::findSurfaceLevels(int level) {
	if (level == 1) {
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (flagMap[ind] == 3) {
				surfaceLevels[ind] = 1;
			}
			else {
				surfaceLevels[ind] = 0;
			}
		}
	}
	else {
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (surfaceLevels[ind] == level - 1) {
				for (int dir = 1; dir < nf; dir++) {
					int indBlockShift = bShift[nAll * dir + ind];
					int indShift = shift[nAll * dir + ind];
					if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift] && !blocks[indBlockShift]->boundary[indShift]) {
						if (((BlockFreeSurface*)blocks[indBlockShift])->surfaceLevels[indShift] == 0) {
							((BlockFreeSurface*)blocks[indBlockShift])->surfaceLevels[indShift] = level;
						}
					}
				}
			}
		}
	}
}

void BlockFreeSurface::findFilledEmptied(int iter) {
	double tolerance = 0.05;
	int fsStatus = iter % 2;
	nFilled = 0;
	nEmptied = 0;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3) { // free surface cell
			double rho = rho0;
			for (int p = 0; p < nf; p++) {
				double fVal = fc[nf * ind + p];
				rho += fVal;
			}
			vof[ind] = mass[ind] / rho;
			if (fsStatus == 0) { //filling 
				if (vof[ind] > 1.0 + tolerance) {
					filled[nFilled] = ind; // cell is marked as filled
					++nFilled;
				}
			}
			else if (fsStatus == 1) {//emptying
				if (vof[ind] < 0.0 - tolerance) {
					emptied[nEmptied] = ind; // cell is marked as emptied
					++nEmptied;
				}
			}
		}
	}
}

void BlockFreeSurface::vicinityFilledEmptied() { // kriticka sekce dochazi k pristupu vice vlaken
	for (int cell = 0; cell < nFilled; cell++) {
		int ind = filled[cell];
		for (int dir = 1; dir < nf; dir++) {
			int indBlockShift = bShift[nAll * dir + ind];
			int indShift = shift[nAll * dir + ind];
			if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift] && !blocks[indBlockShift]->boundary[indShift]) {
				if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 2) {
#pragma omp critical(CriticalSection)
					{
						((BlockFreeSurface*)blocks[indBlockShift])->airToFs[((BlockFreeSurface*)blocks[indBlockShift])->nAirToFs++] = indShift; // cell is marked as airToFs
					}
				}
			}
		}
	}

	//nFluidToFs = 0;
	for (int cell = 0; cell < nEmptied; cell++) {
		int ind = emptied[cell];
		for (int dir = 1; dir < nf; dir++) {
			int indBlockShift = bShift[nAll * dir + ind];
			int indShift = shift[nAll * dir + ind];
			if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift] && !blocks[indBlockShift]->boundary[indShift]) {
				if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 1) {
#pragma omp critical(CriticalSection)
					{
						((BlockFreeSurface*)blocks[indBlockShift])->fluidToFs[((BlockFreeSurface*)blocks[indBlockShift])->nFluidToFs++] = indShift; // cell is marked as fluidToFs
					}
				}
			}
		}
	}
}

void BlockFreeSurface::calculateNormalVectors() {
	double normTol = 1e-5;
	int indBlockShift = 0;
	int indShift = 0;
	double* auxVof = new double[nf];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (surfaceLevels[ind] > 0 && surfaceLevels[ind] < 3) {
			for (int dir = 1; dir < nf; dir++) {
				indBlockShift = bShift[nAll * dir + ind];
				indShift = shift[nAll * dir + ind];
				bool dirExist = false;
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (!(blocks[indBlockShift]->boundary[indShift])) {
						auxVof[dir] = ((BlockFreeSurface*)blocks[indBlockShift])->vofS[indShift];
						dirExist = true;
					}
				}
				if (!dirExist) {
					indBlockShift = bShift[nAll * opp[dir] + ind];
					indShift = shift[nAll * opp[dir] + ind];
					bool oposDirExist = false;
					if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
						if (!(blocks[indBlockShift]->boundary[indShift])) {
							auxVof[dir] = 2 * vofS[ind] - ((BlockFreeSurface*)blocks[indBlockShift])->vofS[indShift];
							oposDirExist = true;
						}
					}
					if (!oposDirExist) {
						auxVof[dir] = vofS[ind];
					}
				}
			}

			double normX = 0;
			double normY = 0;
			double normZ = 0;
			for (int dir = 1; dir < nf; dir++) {
				normX += w[dir] * ex[dir] * auxVof[dir];
				normY += w[dir] * ey[dir] * auxVof[dir];
				normZ += w[dir] * ez[dir] * auxVof[dir];
			}
			double nn = sqrt(normX * normX + normY * normY + normZ * normZ);
			if (nn > normTol) {
				normal[3 * ind + 0] = normX / -nn;
				normal[3 * ind + 1] = normY / -nn;
				normal[3 * ind + 2] = normZ / -nn;
			}
			else {
				normal[3 * ind + 0] = 0.0;
				normal[3 * ind + 1] = 0.0;
				normal[3 * ind + 2] = 0.0;
			}
		}
		else {
			normal[3 * ind + 0] = 0.0;
			normal[3 * ind + 1] = 0.0;
			normal[3 * ind + 2] = 0.0;
		}
	}
	delete[] auxVof;
}

void BlockFreeSurface::setWetAngle() {
	if (wetAngleEq != 0.0 && surfTension != 0.0) {
		double wetAngleDyn;
		double eps = 1e-4;
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (!boundary[ind] && surfaceLevels[ind] > 0) {
				int tmpWallNormal[3] = { 0, 0, 0 };
				int numWallCells = 0;
				int numFillCellsNearWall = 0;
				for (int dir = 1; dir < nf; dir++) {
					int indBlockShift = bShift[nAll * dir + ind];
					int indShift = shift[nAll * dir + ind];
					if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
						if (blocks[indBlockShift]->boundary[indShift]) {
							tmpWallNormal[0] += ex[dir];
							tmpWallNormal[1] += ey[dir];
							tmpWallNormal[2] += ez[dir];
							numWallCells++;
						}
						else if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 1) {
							for (int dir2 = 1; dir2 < nf; dir2++) {
								int indBlockShift2 = ((BlockFreeSurface*)blocks[indBlockShift])->bShift[nAll * dir2 + indShift];
								int indShift2 = ((BlockFreeSurface*)blocks[indBlockShift])->shift[nAll * dir2 + indShift];
								if (!(indBlockShift2 > -1 && indBlockShift2 < nBlockMax && blocks[indBlockShift2])) {
									numFillCellsNearWall++;
									break;
								}
								else if (blocks[indBlockShift2]->boundary[indShift2]) {
									numFillCellsNearWall++;
									break;
								}
							}
						}
					}
					else {
						tmpWallNormal[0] += ex[dir];
						tmpWallNormal[1] += ey[dir];
						tmpWallNormal[2] += ez[dir];
						numWallCells++;
					}
				}

				double wallNormal[3];
				bool isContactCell = false;
				//if (numWallCells != 0 && numFillCellsNearWall != 0) { //&& numEmptyCells != 0) {
				if (numWallCells != 0) {
					int norm2tmpWallNormal = tmpWallNormal[0] * tmpWallNormal[0] + tmpWallNormal[1] * tmpWallNormal[1] + tmpWallNormal[2] * tmpWallNormal[2];
					if (norm2tmpWallNormal != 0) {
						double tmpNormWallNormal = sqrt((double)norm2tmpWallNormal);
						wallNormal[0] = -(double)tmpWallNormal[0] / tmpNormWallNormal;
						wallNormal[1] = -(double)tmpWallNormal[1] / tmpNormWallNormal;
						wallNormal[2] = -(double)tmpWallNormal[2] / tmpNormWallNormal;
						isContactCell = true;
					}
				}

				if (isContactCell) {
					double wallTang[3];
					double cosPhi = normal[3 * ind] * wallNormal[0] + normal[3 * ind + 1] * wallNormal[1] + normal[3 * ind + 2] * wallNormal[2];
					wallTang[0] = normal[3 * ind] - cosPhi * wallNormal[0];
					wallTang[1] = normal[3 * ind + 1] - cosPhi * wallNormal[1];
					wallTang[2] = normal[3 * ind + 2] - cosPhi * wallNormal[2];
					double normWallTang = wallTang[0] * wallTang[0] + wallTang[1] * wallTang[1] + wallTang[2] * wallTang[2];
					if (normWallTang > 0.00001) {
						double sgrtWT = sqrt(normWallTang);
						wallTang[0] = wallTang[0] / sgrtWT;
						wallTang[1] = wallTang[1] / sgrtWT;
						wallTang[2] = wallTang[2] / sgrtWT;

						double normalVelMag = 0.0;
						if (flagMap[ind] == 2) {
							int count = 0;
							double weighTot = 0.0;
							for (int dir = 1; dir < nf; dir++) {
								int indBlockShift = bShift[nAll * dir + ind];
								int indShift = shift[nAll * dir + ind];
								if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
									if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 1 ||
										((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 3) {
										count++;
										double weight = max(((BlockFreeSurface*)blocks[indBlockShift])->vof[indShift], 0.0);
										normalVelMag += (blocks[indBlockShift]->U[indShift] * wallTang[0] +
											blocks[indBlockShift]->V[indShift] * wallTang[1] +
											blocks[indBlockShift]->W[indShift] * wallTang[2]) * weight;
										weighTot += weight;
									}
								}
							}
							if (weighTot != 0.0) {
								normalVelMag /= weighTot;
							}
							else {
								normalVelMag = 0.0;
							}
						}
						else {
							normalVelMag = U[ind] * wallTang[0] + V[ind] * wallTang[1] + W[ind] * wallTang[2];
						}

						double Ca = vLB * normalVelMag / surfTension;
						double ka = 9e-9;
						double kr = 9e-8;
						if (normalVelMag >= 0.0) {
							wetAngleDyn = min(wetAngleEq + cbrt(Ca / ka) * (3.14159265358979323844 / 180.0), wetAngleAdv);
						}
						else {
							wetAngleDyn = max(wetAngleEq + cbrt(Ca / kr) * (3.14159265358979323844 / 180.0), wetAngleRec);
						}
						normal[3 * ind] = cos(wetAngleDyn) * wallNormal[0] + sin(wetAngleDyn) * wallTang[0];
						normal[3 * ind + 1] = cos(wetAngleDyn) * wallNormal[1] + sin(wetAngleDyn) * wallTang[1];
						normal[3 * ind + 2] = cos(wetAngleDyn) * wallNormal[2] + sin(wetAngleDyn) * wallTang[2];
					}
				}
			}
		}
	}
}

void BlockFreeSurface::smoothVOF(int smoothStep, int nSmoothStep) {
	switch (nSmoothStep) {
	case 0:
		vofS = vof;
		break;
	case 1:
		switch (smoothStep) {
		case 0:
			break;
		case 1:
			smooth3D(vof, vofS, 1);
			break;
		}
		break;
	case 2:
		switch (smoothStep) {
		case 0:
			break;
		case 1:
			smooth3D(vof, vofS1, 1);
			break;
		case 2:
			smooth3D(vofS1, vofS, 2);
			break;
		}
		break;
	}
}

void BlockFreeSurface::smooth3D(double* arrayin, double* arrayout, int kod) {
	int indBlockShift = 0;
	int indShift = 0;
	double* auxVof = new double[nf];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (surfaceLevels[ind] > 0) {
			for (int dir = 0; dir < nf; dir++) {
				indBlockShift = bShift[nAll * dir + ind];
				indShift = shift[nAll * dir + ind];
				bool dirExist = false;
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (!(blocks[indBlockShift]->boundary[indShift])) {
						if (kod == 1) {
							auxVof[dir] = ((BlockFreeSurface*)blocks[indBlockShift])->vof[indShift];
						}
						else {
							auxVof[dir] = ((BlockFreeSurface*)blocks[indBlockShift])->vofS1[indShift];
						}
						dirExist = true;
					}
				}
				if (!dirExist) {
					indBlockShift = bShift[nAll * opp[dir] + ind];
					indShift = shift[nAll * opp[dir] + ind];
					bool oposDirExist = false;
					if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
						if (!(blocks[indBlockShift]->boundary[indShift])) {
							if (kod == 1) {
								auxVof[dir] = 2 * vof[ind] - ((BlockFreeSurface*)blocks[indBlockShift])->vof[indShift];
							}
							else {
								auxVof[dir] = 2 * vofS1[ind] - ((BlockFreeSurface*)blocks[indBlockShift])->vofS1[indShift];
							}
							oposDirExist = true;
						}
					}
					if (!oposDirExist) {
						if (kod == 1) {
							auxVof[dir] = vof[ind];
						}
						else {
							auxVof[dir] = vofS1[ind];
						}
					}
				}
			}
			double avg = 0.0;
			double count = 0;
			for (int dir = 0; dir < nf; dir++) {
				avg += w[dir] * auxVof[dir];
				count += w[dir];
			}
			arrayout[ind] = avg / count;
		}
	}
}

void BlockFreeSurface::calculateCurvature() {
	if (surfTension != 0.0) {
		int indBlockShift = 0;
		int indShift = 0;
		double* auxNx = new double[nf];
		double* auxNy = new double[nf];
		double* auxNz = new double[nf];
		double curv;
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (flagMap[ind] == 3) {
				for (int dir = 1; dir < nf; dir++) {
					indBlockShift = bShift[nAll * dir + ind];
					indShift = shift[nAll * dir + ind];
					bool dirExist = false;
					if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
						if (!(blocks[indBlockShift]->boundary[indShift])) {
							auxNx[dir] = ((BlockFreeSurface*)blocks[indBlockShift])->normal[3 * indShift];
							auxNy[dir] = ((BlockFreeSurface*)blocks[indBlockShift])->normal[3 * indShift + 1];
							auxNz[dir] = ((BlockFreeSurface*)blocks[indBlockShift])->normal[3 * indShift + 2];
							dirExist = true;
						}
					}
					if (!dirExist) {
						indBlockShift = bShift[nAll * opp[dir] + ind];
						indShift = shift[nAll * opp[dir] + ind];
						bool oposDirExist = false;
						if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
							if (!(blocks[indBlockShift]->boundary[indShift])) {
								auxNx[dir] = 2 * normal[3 * ind] - ((BlockFreeSurface*)blocks[indBlockShift])->normal[3 * indShift];
								auxNy[dir] = 2 * normal[3 * ind + 1] - ((BlockFreeSurface*)blocks[indBlockShift])->normal[3 * indShift + 1];
								auxNz[dir] = 2 * normal[3 * ind + 2] - ((BlockFreeSurface*)blocks[indBlockShift])->normal[3 * indShift + 2];
								double norm = sqrt(auxNx[dir] * auxNx[dir] + auxNy[dir] * auxNy[dir] + auxNz[dir] * auxNz[dir]);
								if (norm > 1e-5) {
									auxNx[dir] /= norm;
									auxNy[dir] /= norm;
									auxNz[dir] /= norm;
								}
								oposDirExist = true;
							}
						}
						if (!oposDirExist) {
							auxNx[dir] = normal[3 * ind];
							auxNy[dir] = normal[3 * ind + 1];
							auxNz[dir] = normal[3 * ind + 2];
						}
					}
				}

				curv = 0.0;
				for (int dir = 1; dir < nf; dir++) {
					curv += w[dir] * (ex[dir] * auxNx[dir] + ey[dir] * auxNy[dir] + ez[dir] * auxNz[dir]);
				}
				curv *= -3.0;

				if (fabs(curv) > 4.0) {
					if (curv < 0.0) {
						curv = -4.0;
					}
					else {
						curv = 4.0;
					}
				}
				curvature[ind] = curv / 2.0; /// ?????????????
			}
		}
		delete[] auxNx;
		delete[] auxNy;
		delete[] auxNz;
	}
}

void BlockFreeSurface::newFsGeometry() {
	int indBlockShift = 0;
	int indShift = 0;
	double feq[27];
	double* avgRho = new double[nAirToFs];
	double* avgU = new double[nAirToFs];
	double* avgV = new double[nAirToFs];
	double* avgW = new double[nAirToFs];
	for (int cell = 0; cell < nAirToFs; cell++) {
		double count = 0.0;
		avgRho[cell] = 0.0;
		avgU[cell] = 0.0;
		avgV[cell] = 0.0;
		avgW[cell] = 0.0;
		for (int dir = 1; dir < nf; dir++) {
			int ind = airToFs[cell];
			indBlockShift = bShift[nAll * dir + ind];
			indShift = shift[nAll * dir + ind];
			if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
				if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 1 || ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 3) {
					double weight = w[dir] * ((BlockFreeSurface*)blocks[indBlockShift])->vof[indShift];
					count += weight;
					avgRho[cell] += ((BlockFreeSurface*)blocks[indBlockShift])->Rho[indShift] * weight;
					avgU[cell] += ((BlockFreeSurface*)blocks[indBlockShift])->U[indShift] * weight;
					avgV[cell] += ((BlockFreeSurface*)blocks[indBlockShift])->V[indShift] * weight;
					avgW[cell] += ((BlockFreeSurface*)blocks[indBlockShift])->W[indShift] * weight;
				}
			}
		}
		if (count == 0) {
			count = 1;
			cout << "warning: changing air to fs error " << endl;
			//cin >> count;
			avgRho[cell] = 1.0;
			avgU[cell] = 0;
			avgV[cell] = 0;
			avgW[cell] = 0;
		}
		else {
			avgRho[cell] /= count;
			avgU[cell] /= count;
			avgV[cell] /= count;
			avgW[cell] /= count;
		}
	}

	for (int cell = 0; cell < nAirToFs; cell++) {
		int ind = airToFs[cell];
		flagMap[ind] = 3;
		vof[ind] = 0.0;
		mass[ind] = 0.0;
		Rho[ind] = avgRho[cell];
		dRho[ind] = Rho[ind] - rho0;
		U[ind] = avgU[cell];
		V[ind] = avgV[cell];
		W[ind] = avgW[cell];
		equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind]);
		for (int p = 0; p < nf; p++) { // initializace equlibriem
			int pind = nf * ind + p;
			fc[pind] = feq[p];
		}
	}

	for (int cell = 0; cell < nFilled; cell++) {
		int ind = filled[cell];
		flagMap[ind] = 1;
		vof[ind] = 1.0;
		curvature[ind] = 0.0;
	}

	for (int cell = 0; cell < nEmptied; cell++) {
		int ind = emptied[cell];
		flagMap[ind] = 2;
		vof[ind] = 0.0;
		Rho[ind] = 0.0;
		//        U [ind] = 0.0;
		//        V [ind] = 0.0;
		//        W [ind] = 0.0;
		curvature[ind] = 0.0;
	}

	for (int cell = 0; cell < nFluidToFs; cell++) {
		int ind = fluidToFs[cell];
		flagMap[ind] = 3;
		mass[ind] = Rho[ind];
	}

	delete[] avgRho;
	delete[] avgU;
	delete[] avgV;
	delete[] avgW;
}

void BlockFreeSurface::findArtefactsEmptied() {
	nArtefacts1 = 0;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3) {
			bool isArtefact = true;
			int nAir = 0;
			for (int dir = 1; dir < nf; dir++) {
				int indBlockShift = bShift[nAll * dir + ind];
				int indShift = shift[nAll * dir + ind];
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 1) {
						isArtefact = false;
					}
					else if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 2) {
						nAir++;
					}
				}
			}
			if (isArtefact == true) {
				for (int j = 0; j < nArtefacts2; j++) { // if not filled artefact
					if (ind == artefact2[j]) {
						isArtefact = false;
						break;
					}
				}
			}
			if (isArtefact == true && nAir > 0) {
				artefact1[nArtefacts1] = ind; // cell is marked as artefact
				++nArtefacts1;
			}
			else {
			}
		}
	}
}

void BlockFreeSurface::findArtefactsFilled() {
	nArtefacts2 = 0;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3) {
			bool isArtefact = true;
			for (int dir = 1; dir < nf; dir++) {
				double weight = normal[3 * ind + 0] * ex[dir] + normal[3 * ind + 1] * ey[dir] + normal[3 * ind + 2] * ez[dir];
				int indBlockShift = bShift[nAll * dir + ind];
				int indShift = shift[nAll * dir + ind];
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 2) {
						isArtefact = false;
					}
					else if (((BlockFreeSurface*)blocks[indBlockShift])->boundary[ind]) {
						if (weight > 0.6) {
							isArtefact = false;
						}
					}
				}
				else {
					if (weight > 0.6) {
						isArtefact = false;
					}
				}
			}
			if (isArtefact == true) {
				artefact2[nArtefacts2] = ind; // cell is marked as artefact
				++nArtefacts2;
			}
		}
	}
}

void BlockFreeSurface::destroyArtefacts() {
	for (int cell = 0; cell < nArtefacts2; cell++) {
		int ind = artefact2[cell];
		flagMap[ind] = 1;
		vof[ind] = 1.0;
		curvature[ind] = 0.0;
	}
	for (int cell = 0; cell < nArtefacts1; cell++) {
		int ind = artefact1[cell];
		flagMap[ind] = 2;
		vof[ind] = 0.0;
		Rho[ind] = 0.0;
		curvature[ind] = 0.0;
	}
}

void BlockFreeSurface::massRedistribution() {
	for (int cell = 0; cell < nEmptied; cell++) {
		int ind = emptied[cell];
		redistribute(ind, mass[ind]);
	}
	for (int cell = 0; cell < nFilled; cell++) {
		int ind = filled[cell];
		redistribute(ind, mass[ind] - Rho[ind]);
	}
	for (int cell = 0; cell < nArtefacts1; cell++) {
		int ind = artefact1[cell];
		redistribute(ind, mass[ind]);
	}
	for (int cell = 0; cell < nArtefacts2; cell++) {
		int ind = artefact2[cell];
		redistribute(ind, mass[ind] - Rho[ind]);
	}
}

void BlockFreeSurface::redistribute(int ind, double massRedist) {
	int indNeigh[26];
	int indBlockNeigh[26];
	double listWi[26];
	int countNeigh = 0;
	double wtot = 0.0;
	int indBlockShift = 0;
	int indShift = 0;

	for (int dx = -1; dx < 2; dx++) {
		for (int dy = -1; dy < 2; dy++) {
			for (int dz = -1; dz < 2; dz++) {
				local2globalShift(indBlock, ind, indBlockShift, indShift, dx, dy, dz);
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					double wi = 1; //max(0.0, (U[ind] * dx + V[ind] * dy + W[ind] * dz));
					if (((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 3 && wi > 0.0) {
						indNeigh[countNeigh] = indShift;
						indBlockNeigh[countNeigh] = indBlockShift;
						listWi[countNeigh] = wi;
						countNeigh++;
						wtot += wi;
					}
				}
			}
		}
	}
	if (countNeigh == 0) {
		massErr += massRedist;
	}
	else {
		for (int ng = 0; ng < countNeigh; ng++) {
			((BlockFreeSurface*)blocks[indBlockNeigh[ng]])->mass[indNeigh[ng]] += (listWi[ng] / wtot) * massRedist;
		}
	}
}

void BlockFreeSurface::calculMassTotal() {
	massBlock = 0.0;
	activeCellsNumBlock = 0;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3) {
			massBlock += mass[ind];
			activeCellsNumBlock++;
		}
		else if (flagMap[ind] == 1) {
			massBlock += Rho[ind];
			activeCellsNumBlock++;
		}
	}
	massBlock *= Crho * Ch * Ch * Ch;
}

void BlockFreeSurface::checkIfBlockEmpty() {
	isBlockEmpty = true;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3 || flagMap[ind] == 1) {
			isBlockEmpty = false;
			return;
		}
	}
}

void BlockFreeSurface::touchNeighBlocks() {
	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		if (surfaceLevels[ind] > 0) {
			for (int dir = 1; dir < nf; dir++) {
				int indBlockShift = bShift[nAll * dir + ind];
				int indShift = shift[nAll * dir + ind];
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (indBlockShift != indBlock) {
						((BlockFreeSurface*)blocks[indBlockShift])->isBlockTouched = true;
					}
				}
			}
		}
	}
}

double BlockFreeSurface::computeShearRate(int ind) {
	double dudx = 0;
	double dudy = 0;
	double dudz = 0;
	double dvdx = 0;
	double dvdy = 0;
	double dvdz = 0;
	double dwdx = 0;
	double dwdy = 0;
	double dwdz = 0;
	int indBlockShift;
	int indShift;
	double gamma = 0;
	if (flagMap[ind] == 1 || flagMap[ind] == 3) {
		for (int dir = 1; dir < nf; dir++) {
			indBlockShift = bShift[nAll * dir + ind];
			indShift = shift[nAll * dir + ind];
			if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
				if (!(blocks[indBlockShift]->boundary[indShift]) && ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] != 2) {
					double wu = 3 * w[dir] * ((BlockFreeSurface*)blocks[indBlockShift])->U[indShift];
					dudx += wu * ex[dir];
					dudy += wu * ey[dir];
					dudz += wu * ez[dir];

					double wv = 3 * w[dir] * ((BlockFreeSurface*)blocks[indBlockShift])->V[indShift];
					dvdx += wv * ex[dir];
					dvdy += wv * ey[dir];
					dvdz += wv * ez[dir];

					double ww = 3 * w[dir] * ((BlockFreeSurface*)blocks[indBlockShift])->W[indShift];
					dwdx += ww * ex[dir];
					dwdy += ww * ey[dir];
					dwdz += ww * ez[dir];
				}
				else if (!(blocks[indBlockShift]->boundary[indShift]) && ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift] == 2) {
					double wu = 3 * w[dir] * U[ind];
					dudx += wu * ex[dir];
					dudy += wu * ey[dir];
					dudz += wu * ez[dir];

					double wv = 3 * w[dir] * V[ind];
					dvdx += wv * ex[dir];
					dvdy += wv * ey[dir];
					dvdz += wv * ez[dir];

					double ww = 3 * w[dir] * W[ind];
					dwdx += ww * ex[dir];
					dwdy += ww * ey[dir];
					dwdz += ww * ez[dir];
				}
				else {
					double wu = -3 * w[dir] * U[ind];
					dudx += wu * ex[dir];
					dudy += wu * ey[dir];
					dudz += wu * ez[dir];

					double wv = -3 * w[dir] * V[ind];
					dvdx += wv * ex[dir];
					dvdy += wv * ey[dir];
					dvdz += wv * ez[dir];

					double ww = -3 * w[dir] * W[ind];
					dwdx += ww * ex[dir];
					dwdy += ww * ey[dir];
					dwdz += ww * ez[dir];
				}
			}
			else {
				double wu = -3 * w[dir] * U[ind];
				dudx += wu * ex[dir];
				dudy += wu * ey[dir];
				dudz += wu * ez[dir];

				double wv = -3 * w[dir] * V[ind];
				dvdx += wv * ex[dir];
				dvdy += wv * ey[dir];
				dvdz += wv * ez[dir];

				double ww = -3 * w[dir] * W[ind];
				dwdx += ww * ex[dir];
				dwdy += ww * ey[dir];
				dwdz += ww * ez[dir];
			}
		}
		gamma = sqrt(2 * (dudx * dudx + dvdy * dvdy + dwdz * dwdz) + (dudy + dvdx) * (dudy + dvdx) + (dudz + dwdx) * (dudz + dwdx) + (dvdz + dwdy) * (dvdz + dwdy));
	}

	return gamma;
}

void BlockFreeSurface::writeBlock(int id) {
	ofstream myfile;
	myfile.open("results/results" + std::to_string(id) + ".txt", ios::out);
	myfile << ni << " " << nj << " " << nk << "\n";
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			for (int k = 0; k < nk; k++) {
				int ind = (ni * nj) * k + ni * j + i;
				myfile << Rho[ind] << " " << U[ind] << " " << V[ind] << " " << W[ind] << " " << vof[ind] << " " << boundary[ind] << " " << flagMap[ind] << "\n";
			}
		}
	}
	myfile.close();
}

void BlockFreeSurface::defineMRTmatrix() {
	if (nf == 15) { //D3Q15
		M = new double[nf * nf]{
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


		iM = new double[nf * nf]{
			0.0666666666666667, -0.1111111111111111, 0.0444444444444444, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000,
			0.0666666666666667, -0.0555555555555556, -0.0111111111111111, 0.1000000000000000, -0.1000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.1666666666666667, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000,
			0.0666666666666667, -0.0555555555555556, -0.0111111111111111, -0.1000000000000000, 0.1000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.1666666666666667, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000,
			0.0666666666666667, -0.0555555555555555, -0.0111111111111111, 0.0000000000000000, 0.0000000000000000, 0.1000000000000000, -0.1000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0833333333333333, 0.2500000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000,
			0.0666666666666667, -0.0555555555555556, -0.0111111111111111, -0.0000000000000000, 0.0000000000000000, -0.1000000000000000, 0.1000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0833333333333333, 0.2500000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000,
			0.0666666666666667, -0.0555555555555555, -0.0111111111111111, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.1000000000000000, -0.1000000000000000, -0.0833333333333333, -0.2500000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000,
			0.0666666666666667, -0.0555555555555556, -0.0111111111111111, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.1000000000000000, 0.1000000000000000, -0.0833333333333333, -0.2500000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000,
			0.0666666666666666, 0.0555555555555555, 0.0027777777777778, 0.1000000000000000, 0.0250000000000000, 0.1000000000000001, 0.0250000000000000, 0.1000000000000000, 0.0250000000000000, 0.0000000000000000, 0.0000000000000000, 0.1250000000000000, 0.1250000000000000, 0.1250000000000000, 0.1250000000000000,
			0.0666666666666667, 0.0555555555555556, 0.0027777777777778, -0.1000000000000000, -0.0250000000000000, 0.1000000000000000, 0.0250000000000000, 0.1000000000000000, 0.0250000000000000, -0.0000000000000000, 0.0000000000000000, -0.1250000000000000, 0.1250000000000000, -0.1250000000000000, -0.1250000000000000,
			0.0666666666666667, 0.0555555555555556, 0.0027777777777778, 0.1000000000000000, 0.0250000000000000, -0.1000000000000001, -0.0250000000000000, 0.1000000000000000, 0.0250000000000000, 0.0000000000000000, 0.0000000000000000, -0.1250000000000000, -0.1250000000000000, 0.1250000000000000, -0.1250000000000000,
			0.0666666666666666, 0.0555555555555555, 0.0027777777777778, -0.1000000000000000, -0.0250000000000000, -0.1000000000000000, -0.0250000000000000, 0.1000000000000000, 0.0250000000000000, 0.0000000000000000, 0.0000000000000000, 0.1249999999999999, -0.1250000000000000, -0.1250000000000000, 0.1250000000000000,
			0.0666666666666667, 0.0555555555555556, 0.0027777777777778, 0.1000000000000000, 0.0250000000000000, 0.1000000000000000, 0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.0000000000000000, 0.0000000000000000, 0.1250000000000000, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000,
			0.0666666666666667, 0.0555555555555556, 0.0027777777777778, -0.1000000000000000, -0.0250000000000000, 0.1000000000000000, 0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.0000000000000000, 0.0000000000000000, -0.1250000000000000, -0.1250000000000000, 0.1250000000000000, 0.1250000000000000,
			0.0666666666666667, 0.0555555555555556, 0.0027777777777778, 0.1000000000000000, 0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.0000000000000000, 0.0000000000000000, -0.1250000000000000, 0.1250000000000000, -0.1250000000000000, 0.1250000000000000,
			0.0666666666666667, 0.0555555555555556, 0.0027777777777778, -0.1000000000000000, -0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.0000000000000000, 0.0000000000000000, 0.1250000000000000, 0.1250000000000000, 0.1250000000000000, -0.1250000000000000
		};


		double s1 = 1.19;
		double s2 = 1.4;
		double s4 = 1.2;
		double s14 = 1.98;
		double s = 1 / tau;

		S = new double[nf] {
			0, s1, s2, 0, s4, 0, s4, 0, s4, s, s, s, s, s, s14
		};


		T = new double[nf * nf];
		double* aux = new double[nf * nf];
		for (int i = 0; i < nf; i++) {
			for (int j = 0; j < nf; j++) {
				aux[nf * i + j] = S[i] * M[nf * i + j];
			}
		}

		for (int i = 0; i < nf; i++) {
			for (int j = 0; j < nf; j++) {
				T[nf * i + j] = 0;
				for (int k = 0; k < nf; k++) {
					T[nf * i + j] += iM[nf * i + k] * aux[nf * k + j];
				}
			}
		}
		delete[] aux;

		T0 = new double[nf * nf]{
1.2600000000000000, -0.1166666666666667, -0.1166666666666667, -0.1166666666666667, -0.1166666666666667, -0.1166666666666667, -0.1166666666666667, -0.0700000000000000, -0.0700000000000000, -0.0700000000000000, -0.0700000000000000, -0.0700000000000000, -0.0700000000000000, -0.0700000000000000, -0.0700000000000000,
-0.1166666666666666, 0.6083333333333332, -0.3516666666666666, 0.1283333333333333, 0.1283333333333333, 0.1283333333333333, 0.1283333333333333, -0.2016666666666666, 0.0383333333333333, -0.2016666666666666, 0.0383333333333333, -0.2016666666666666, 0.0383333333333333, -0.2016666666666666, 0.0383333333333333,
-0.1166666666666666, -0.3516666666666667, 0.6083333333333333, 0.1283333333333333, 0.1283333333333333, 0.1283333333333333, 0.1283333333333333, 0.0383333333333333, -0.2016666666666667, 0.0383333333333333, -0.2016666666666667, 0.0383333333333333, -0.2016666666666667, 0.0383333333333333, -0.2016666666666667,
-0.1166666666666667, 0.1283333333333333, 0.1283333333333333, 0.6083333333333333, -0.3516666666666667, 0.1283333333333333, 0.1283333333333333, -0.2016666666666667, -0.2016666666666667, 0.0383333333333334, 0.0383333333333334, -0.2016666666666667, -0.2016666666666666, 0.0383333333333334, 0.0383333333333333,
-0.1166666666666666, 0.1283333333333333, 0.1283333333333333, -0.3516666666666667, 0.6083333333333333, 0.1283333333333333, 0.1283333333333333, 0.0383333333333333, 0.0383333333333333, -0.2016666666666667, -0.2016666666666667, 0.0383333333333333, 0.0383333333333333, -0.2016666666666667, -0.2016666666666667,
-0.1166666666666667, 0.1283333333333333, 0.1283333333333333, 0.1283333333333333, 0.1283333333333333, 0.6083333333333333, -0.3516666666666667, -0.2016666666666667, -0.2016666666666667, -0.2016666666666667, -0.2016666666666667, 0.0383333333333334, 0.0383333333333333, 0.0383333333333334, 0.0383333333333334,
-0.1166666666666666, 0.1283333333333333, 0.1283333333333333, 0.1283333333333333, 0.1283333333333333, -0.3516666666666667, 0.6083333333333333, 0.0383333333333333, 0.0383333333333333, 0.0383333333333333, 0.0383333333333333, -0.2016666666666667, -0.2016666666666667, -0.2016666666666667, -0.2016666666666667,
-0.0699999999999999, -0.2016666666666666, 0.0383333333333335, -0.2016666666666666, 0.0383333333333335, -0.2016666666666666, 0.0383333333333335, 0.4075000000000000, -0.1475000000000001, -0.1475000000000001, 0.2874999999999999, -0.1475000000000001, 0.2874999999999999, 0.2874999999999999, -0.2675000000000001,
-0.0700000000000000, 0.0383333333333334, -0.2016666666666667, -0.2016666666666667, 0.0383333333333334, -0.2016666666666667, 0.0383333333333334, -0.1475000000000000, 0.4075000000000000, 0.2875000000000000, -0.1475000000000000, 0.2875000000000000, -0.1475000000000000, -0.2675000000000000, 0.2875000000000000,
-0.0700000000000000, -0.2016666666666667, 0.0383333333333332, 0.0383333333333334, -0.2016666666666668, -0.2016666666666667, 0.0383333333333333, -0.1475000000000000, 0.2875000000000000, 0.4075000000000000, -0.1474999999999999, 0.2875000000000000, -0.2675000000000000, -0.1474999999999999, 0.2875000000000000,
-0.0700000000000000, 0.0383333333333333, -0.2016666666666666, 0.0383333333333333, -0.2016666666666666, -0.2016666666666666, 0.0383333333333334, 0.2875000000000000, -0.1475000000000000, -0.1475000000000000, 0.4075000000000000, -0.2675000000000000, 0.2875000000000000, 0.2875000000000000, -0.1475000000000001,
-0.0700000000000001, -0.2016666666666667, 0.0383333333333333, -0.2016666666666667, 0.0383333333333332, 0.0383333333333333, -0.2016666666666668, -0.1475000000000000, 0.2875000000000000, 0.2875000000000000, -0.2675000000000000, 0.4075000000000000, -0.1474999999999999, -0.1474999999999999, 0.2875000000000001,
-0.0700000000000000, 0.0383333333333333, -0.2016666666666667, -0.2016666666666667, 0.0383333333333333, 0.0383333333333333, -0.2016666666666667, 0.2875000000000000, -0.1475000000000000, -0.2675000000000000, 0.2875000000000000, -0.1475000000000000, 0.4075000000000000, 0.2875000000000000, -0.1475000000000000,
-0.0700000000000000, -0.2016666666666667, 0.0383333333333333, 0.0383333333333333, -0.2016666666666667, 0.0383333333333333, -0.2016666666666667, 0.2875000000000000, -0.2675000000000000, -0.1475000000000000, 0.2875000000000000, -0.1475000000000000, 0.2875000000000000, 0.4075000000000000, -0.1475000000000000,
-0.0700000000000000, 0.0383333333333333, -0.2016666666666667, 0.0383333333333333, -0.2016666666666667, 0.0383333333333334, -0.2016666666666667, -0.2675000000000001, 0.2875000000000000, 0.2875000000000000, -0.1475000000000001, 0.2875000000000001, -0.1475000000000000, -0.1475000000000000, 0.4075000000000001
		};

		T1 = new double[nf * nf]{
0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
0.0000000000000000, 0.3333333333333334, 0.3333333333333334, -0.1666666666666667, -0.1666666666666667, -0.1666666666666667, -0.1666666666666667, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
0.0000000000000000, 0.3333333333333334, 0.3333333333333334, -0.1666666666666667, -0.1666666666666667, -0.1666666666666667, -0.1666666666666667, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
0.0000000000000000, -0.1666666666666667, -0.1666666666666667, 0.3333333333333334, 0.3333333333333333, -0.1666666666666667, -0.1666666666666667, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
0.0000000000000000, -0.1666666666666667, -0.1666666666666667, 0.3333333333333333, 0.3333333333333334, -0.1666666666666667, -0.1666666666666667, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
0.0000000000000000, -0.1666666666666667, -0.1666666666666667, -0.1666666666666667, -0.1666666666666667, 0.3333333333333334, 0.3333333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000,
0.0000000000000000, -0.1666666666666667, -0.1666666666666667, -0.1666666666666667, -0.1666666666666667, 0.3333333333333333, 0.3333333333333334, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.3750000000000000, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000, 0.3750000000000000,
0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.1250000000000000, 0.3750000000000000, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000, 0.3750000000000000, -0.1250000000000000,
0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.1250000000000000, -0.1250000000000000, 0.3749999999999999, -0.1249999999999999, -0.1249999999999999, 0.3750000000000000, -0.1250000000000000, -0.1250000000000000,
0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.1250000000000001, -0.1249999999999999, -0.1249999999999999, 0.3750000000000000, 0.3749999999999999, -0.1250000000000000, -0.1249999999999999, -0.1250000000000001,
0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000, 0.3750000000000000, 0.3749999999999999, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000,
0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.1250000000000000, -0.1250000000000000, 0.3750000000000000, -0.1250000000000000, -0.1250000000000000, 0.3749999999999999, -0.1250000000000000, -0.1250000000000000,
0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.1250000000000000, 0.3749999999999999, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000, -0.1250000000000000, 0.3749999999999999, -0.1250000000000000,
0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.3750000000000001, -0.1250000000000000, -0.1250000000000000, -0.1250000000000001, -0.1250000000000001, -0.1250000000000001, -0.1250000000000000, 0.3750000000000000
		};


	}

	if (nf == 19) { //D3Q19
		M = new double[nf * nf]{
			1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
			-30.000, -11.000, -11.000, -11.000, -11.000, -11.000, -11.000, 8.000, 8.000, 8.000, 8.000, 8.000, 8.000, 8.000, 8.000, 8.000, 8.000, 8.000, 8.000,
			12.000, -4.000, -4.000, -4.000, -4.000, -4.000, -4.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
			0.000, 1.000, -1.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, 1.000, -1.000, 1.000, -1.000, 1.000, -1.000, 0.000, 0.000, 0.000, 0.000,
			0.000, -4.000, 4.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, 1.000, -1.000, 1.000, -1.000, 1.000, -1.000, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.000, 0.000, 1.000, -1.000, 0.000, 0.000, 1.000, 1.000, -1.000, -1.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, 1.000, -1.000,
			0.000, 0.000, 0.000, -4.000, 4.000, 0.000, 0.000, 1.000, 1.000, -1.000, -1.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, 1.000, -1.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, 0.000, 0.000, 0.000, 0.000, 1.000, 1.000, -1.000, -1.000, 1.000, 1.000, -1.000, -1.000,
			0.000, 0.000, 0.000, 0.000, 0.000, -4.000, 4.000, 0.000, 0.000, 0.000, 0.000, 1.000, 1.000, -1.000, -1.000, 1.000, 1.000, -1.000, -1.000,
			0.000, 2.000, 2.000, -1.000, -1.000, -1.000, -1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, -2.000, -2.000, -2.000, -2.000,
			0.000, -4.000, -4.000, 2.000, 2.000, 2.000, 2.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, -2.000, -2.000, -2.000, -2.000,
			0.000, 0.000, 0.000, 1.000, 1.000, -1.000, -1.000, 1.000, 1.000, 1.000, 1.000, -1.000, -1.000, -1.000, -1.000, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.000, 0.000, -2.000, -2.000, 2.000, 2.000, 1.000, 1.000, 1.000, 1.000, -1.000, -1.000, -1.000, -1.000, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, -1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, -1.000, 1.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, -1.000, 1.000, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, 1.000, -1.000, -1.000, 1.000, -1.000, 1.000, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -1.000, -1.000, 1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 1.000, -1.000, 1.000, -1.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 1.000, -1.000, -1.000, -1.000, -1.000, 1.000, 1.000
		};

		iM = new double[nf * nf]{
			0.0526315789473684, -0.0125313283208020, 0.0476190476190476, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000,
			0.0526315789473684, -0.0045948203842941, -0.0158730158730159, 0.1000000000000000, -0.1000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0555555555555556, -0.0555555555555556, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
			0.0526315789473684, -0.0045948203842941, -0.0158730158730159, -0.1000000000000000, 0.1000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0555555555555556, -0.0555555555555556, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000,
			0.0526315789473684, -0.0045948203842941, -0.0158730158730159, 0.0000000000000000, 0.0000000000000000, 0.1000000000000000, -0.1000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0277777777777778, 0.0277777777777778, 0.0833333333333333, -0.0833333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000,
			0.0526315789473684, -0.0045948203842941, -0.0158730158730159, 0.0000000000000000, 0.0000000000000000, -0.1000000000000000, 0.1000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0277777777777778, 0.0277777777777778, 0.0833333333333333, -0.0833333333333333, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000,
			0.0526315789473684, -0.0045948203842941, -0.0158730158730159, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.1000000000000000, -0.1000000000000000, -0.0277777777777778, 0.0277777777777778, -0.0833333333333333, 0.0833333333333333, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000,
			0.0526315789473684, -0.0045948203842941, -0.0158730158730159, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.1000000000000000, 0.1000000000000000, -0.0277777777777778, 0.0277777777777778, -0.0833333333333333, 0.0833333333333333, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000,
			0.0526315789473683, 0.0033416875522139, 0.0039682539682540, 0.1000000000000000, 0.0250000000000000, 0.1000000000000000, 0.0250000000000000, 0.0000000000000000, 0.0000000000000000, 0.0277777777777778, 0.0138888888888889, 0.0833333333333333, 0.0416666666666667, 0.2500000000000000, -0.0000000000000001, -0.0000000000000001, 0.1250000000000000, -0.1250000000000001, 0.0000000000000000,
			0.0526315789473684, 0.0033416875522139, 0.0039682539682540, -0.1000000000000000, -0.0250000000000000, 0.1000000000000000, 0.0250000000000000, 0.0000000000000000, 0.0000000000000000, 0.0277777777777778, 0.0138888888888889, 0.0833333333333334, 0.0416666666666667, -0.2500000000000000, -0.0000000000000000, 0.0000000000000000, -0.1250000000000000, -0.1250000000000000, -0.0000000000000000,
			0.0526315789473685, 0.0033416875522139, 0.0039682539682540, 0.1000000000000000, 0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.0000000000000000, -0.0000000000000000, 0.0277777777777778, 0.0138888888888889, 0.0833333333333333, 0.0416666666666667, -0.2500000000000000, 0.0000000000000000, 0.0000000000000001, 0.1250000000000001, 0.1250000000000000, 0.0000000000000000,
			0.0526315789473685, 0.0033416875522139, 0.0039682539682540, -0.1000000000000000, -0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.0000000000000000, -0.0000000000000000, 0.0277777777777778, 0.0138888888888889, 0.0833333333333333, 0.0416666666666667, 0.2500000000000000, 0.0000000000000000, 0.0000000000000000, -0.1250000000000000, 0.1250000000000000, 0.0000000000000000,
			0.0526315789473684, 0.0033416875522139, 0.0039682539682540, 0.1000000000000001, 0.0250000000000000, -0.0000000000000000, -0.0000000000000000, 0.1000000000000000, 0.0250000000000000, 0.0277777777777777, 0.0138888888888889, -0.0833333333333333, -0.0416666666666667, 0.0000000000000000, 0.0000000000000000, 0.2499999999999999, -0.1250000000000001, 0.0000000000000000, 0.1250000000000001,
			0.0526315789473684, 0.0033416875522139, 0.0039682539682540, -0.1000000000000000, -0.0250000000000000, 0.0000000000000000, 0.0000000000000000, 0.1000000000000000, 0.0250000000000000, 0.0277777777777778, 0.0138888888888889, -0.0833333333333333, -0.0416666666666667, 0.0000000000000000, 0.0000000000000000, -0.2500000000000000, 0.1250000000000000, -0.0000000000000000, 0.1250000000000000,
			0.0526315789473684, 0.0033416875522139, 0.0039682539682540, 0.1000000000000000, 0.0250000000000000, 0.0000000000000000, 0.0000000000000000, -0.1000000000000000, -0.0250000000000000, 0.0277777777777778, 0.0138888888888889, -0.0833333333333333, -0.0416666666666667, 0.0000000000000000, -0.0000000000000000, -0.2499999999999999, -0.1250000000000000, 0.0000000000000000, -0.1250000000000000,
			0.0526315789473684, 0.0033416875522139, 0.0039682539682540, -0.1000000000000000, -0.0250000000000000, 0.0000000000000000, 0.0000000000000000, -0.1000000000000000, -0.0250000000000000, 0.0277777777777778, 0.0138888888888889, -0.0833333333333333, -0.0416666666666667, 0.0000000000000000, -0.0000000000000000, 0.2499999999999999, 0.1250000000000000, 0.0000000000000000, -0.1250000000000000,
			0.0526315789473684, 0.0033416875522139, 0.0039682539682540, -0.0000000000000000, -0.0000000000000000, 0.1000000000000000, 0.0250000000000000, 0.1000000000000000, 0.0250000000000000, -0.0555555555555556, -0.0277777777777778, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.2500000000000001, 0.0000000000000000, -0.0000000000000000, 0.1250000000000000, -0.1250000000000000,
			0.0526315789473684, 0.0033416875522139, 0.0039682539682540, 0.0000000000000000, 0.0000000000000000, -0.1000000000000000, -0.0250000000000000, 0.1000000000000000, 0.0250000000000000, -0.0555555555555556, -0.0277777777777778, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.2500000000000001, 0.0000000000000000, 0.0000000000000000, -0.1250000000000000, -0.1250000000000000,
			0.0526315789473684, 0.0033416875522139, 0.0039682539682540, -0.0000000000000000, -0.0000000000000000, 0.1000000000000000, 0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.0555555555555555, -0.0277777777777778, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.2500000000000000, -0.0000000000000000, -0.0000000000000000, 0.1250000000000000, 0.1250000000000000,
			0.0526315789473684, 0.0033416875522139, 0.0039682539682540, 0.0000000000000000, 0.0000000000000000, -0.1000000000000000, -0.0250000000000000, -0.1000000000000000, -0.0250000000000000, -0.0555555555555556, -0.0277777777777778, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.2500000000000000, -0.0000000000000000, 0.0000000000000000, -0.1250000000000000, 0.1250000000000000
		};

		double s1 = 1.19;
		double s2 = 1.4;
		double s4 = 1.2;
		double s10 = s2;
		double s16 = 1.98;
		double s = 1 / tau;

		S = new double[nf] {
			0, s1, s2, 0, s4, 0, s4, 0, s4, s, s10, s, s10, s, s, s, s16, s16, s16
		};

		T = new double[nf * nf];
		double* aux = new double[nf * nf];
		for (int i = 0; i < nf; i++) {
			for (int j = 0; j < nf; j++) {
				aux[nf * i + j] = S[i] * M[nf * i + j];
			}
		}

		for (int i = 0; i < nf; i++) {
			for (int j = 0; j < nf; j++) {
				T[nf * i + j] = 0;
				for (int k = 0; k < nf; k++) {
					T[nf * i + j] += iM[nf * i + k] * aux[nf * k + j];
				}
			}
		}
		delete[] aux;

		T0 = new double[nf * nf]{
1.2473684210526312, -0.1026315789473684, -0.1026315789473684, -0.1026315789473684, -0.1026315789473684, -0.1026315789473684, -0.1026315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684, -0.0526315789473684,
-0.1026315789473685, 0.9401461988304092, -0.0198538011695907, -0.0065204678362573, -0.0065204678362573, -0.0065204678362573, -0.0065204678362573, -0.2637426900584795, -0.0237426900584795, -0.2637426900584795, -0.0237426900584795, -0.2637426900584795, -0.0237426900584795, -0.2637426900584795, -0.0237426900584795, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538,
-0.1026315789473685, -0.0198538011695907, 0.9401461988304092, -0.0065204678362573, -0.0065204678362573, -0.0065204678362573, -0.0065204678362573, -0.0237426900584795, -0.2637426900584795, -0.0237426900584795, -0.2637426900584795, -0.0237426900584795, -0.2637426900584795, -0.0237426900584795, -0.2637426900584795, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538,
-0.1026315789473685, -0.0065204678362573, -0.0065204678362573, 0.9401461988304092, -0.0198538011695907, -0.0065204678362573, -0.0065204678362574, -0.2637426900584795, -0.2637426900584795, -0.0237426900584795, -0.0237426900584795, 0.0895906432748538, 0.0895906432748539, 0.0895906432748538, 0.0895906432748539, -0.2637426900584795, -0.0237426900584795, -0.2637426900584795, -0.0237426900584795,
-0.1026315789473684, -0.0065204678362573, -0.0065204678362573, -0.0198538011695906, 0.9401461988304094, -0.0065204678362573, -0.0065204678362574, -0.0237426900584795, -0.0237426900584795, -0.2637426900584795, -0.2637426900584795, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538, -0.0237426900584796, -0.2637426900584795, -0.0237426900584795, -0.2637426900584795,
-0.1026315789473684, -0.0065204678362573, -0.0065204678362573, -0.0065204678362573, -0.0065204678362573, 0.9401461988304094, -0.0198538011695906, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538, -0.2637426900584795, -0.2637426900584795, -0.0237426900584795, -0.0237426900584795, -0.2637426900584795, -0.2637426900584795, -0.0237426900584795, -0.0237426900584795,
-0.1026315789473684, -0.0065204678362573, -0.0065204678362573, -0.0065204678362573, -0.0065204678362573, -0.0198538011695906, 0.9401461988304094, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538, 0.0895906432748538, -0.0237426900584795, -0.0237426900584795, -0.2637426900584795, -0.2637426900584795, -0.0237426900584795, -0.0237426900584795, -0.2637426900584795, -0.2637426900584795,
-0.0526315789473683, -0.2637426900584794, -0.0237426900584794, -0.2637426900584795, -0.0237426900584795, 0.0895906432748538, 0.0895906432748539, 0.6701461988304094, 0.1151461988304094, 0.1151461988304092, -0.4398538011695908, -0.2190204678362574, 0.2159795321637426, -0.2190204678362574, 0.2159795321637426, -0.2190204678362574, 0.2159795321637428, -0.2190204678362575, 0.2159795321637428,
-0.0526315789473684, -0.0237426900584795, -0.2637426900584795, -0.2637426900584795, -0.0237426900584796, 0.0895906432748538, 0.0895906432748539, 0.1151461988304094, 0.6701461988304094, -0.4398538011695907, 0.1151461988304093, 0.2159795321637427, -0.2190204678362575, 0.2159795321637427, -0.2190204678362574, -0.2190204678362574, 0.2159795321637428, -0.2190204678362574, 0.2159795321637427,
-0.0526315789473685, -0.2637426900584795, -0.0237426900584796, -0.0237426900584795, -0.2637426900584796, 0.0895906432748538, 0.0895906432748537, 0.1151461988304095, -0.4398538011695907, 0.6701461988304095, 0.1151461988304093, -0.2190204678362573, 0.2159795321637429, -0.2190204678362575, 0.2159795321637428, 0.2159795321637427, -0.2190204678362573, 0.2159795321637427, -0.2190204678362572,
-0.0526315789473685, -0.0237426900584796, -0.2637426900584796, -0.0237426900584795, -0.2637426900584795, 0.0895906432748538, 0.0895906432748537, -0.4398538011695906, 0.1151461988304094, 0.1151461988304093, 0.6701461988304093, 0.2159795321637428, -0.2190204678362572, 0.2159795321637427, -0.2190204678362573, 0.2159795321637426, -0.2190204678362573, 0.2159795321637427, -0.2190204678362572,
-0.0526315789473684, -0.2637426900584795, -0.0237426900584794, 0.0895906432748537, 0.0895906432748537, -0.2637426900584796, -0.0237426900584796, -0.2190204678362574, 0.2159795321637428, -0.2190204678362574, 0.2159795321637428, 0.6701461988304096, 0.1151461988304093, 0.1151461988304093, -0.4398538011695909, -0.2190204678362573, -0.2190204678362574, 0.2159795321637429, 0.2159795321637428,
-0.0526315789473684, -0.0237426900584795, -0.2637426900584796, 0.0895906432748538, 0.0895906432748538, -0.2637426900584795, -0.0237426900584795, 0.2159795321637427, -0.2190204678362573, 0.2159795321637427, -0.2190204678362573, 0.1151461988304092, 0.6701461988304094, -0.4398538011695906, 0.1151461988304095, -0.2190204678362573, -0.2190204678362572, 0.2159795321637426, 0.2159795321637426,
-0.0526315789473684, -0.2637426900584796, -0.0237426900584795, 0.0895906432748538, 0.0895906432748538, -0.0237426900584795, -0.2637426900584796, -0.2190204678362573, 0.2159795321637427, -0.2190204678362573, 0.2159795321637427, 0.1151461988304093, -0.4398538011695907, 0.6701461988304094, 0.1151461988304094, 0.2159795321637428, 0.2159795321637427, -0.2190204678362573, -0.2190204678362574,
-0.0526315789473684, -0.0237426900584795, -0.2637426900584795, 0.0895906432748537, 0.0895906432748538, -0.0237426900584795, -0.2637426900584795, 0.2159795321637427, -0.2190204678362573, 0.2159795321637427, -0.2190204678362573, -0.4398538011695906, 0.1151461988304094, 0.1151461988304093, 0.6701461988304093, 0.2159795321637426, 0.2159795321637426, -0.2190204678362572, -0.2190204678362572,
-0.0526315789473684, 0.0895906432748538, 0.0895906432748538, -0.2637426900584795, -0.0237426900584795, -0.2637426900584795, -0.0237426900584796, -0.2190204678362573, -0.2190204678362573, 0.2159795321637427, 0.2159795321637427, -0.2190204678362572, -0.2190204678362573, 0.2159795321637427, 0.2159795321637426, 0.6701461988304093, 0.1151461988304093, 0.1151461988304094, -0.4398538011695906,
-0.0526315789473684, 0.0895906432748538, 0.0895906432748538, -0.0237426900584795, -0.2637426900584796, -0.2637426900584796, -0.0237426900584795, 0.2159795321637427, 0.2159795321637427, -0.2190204678362573, -0.2190204678362573, -0.2190204678362574, -0.2190204678362573, 0.2159795321637427, 0.2159795321637428, 0.1151461988304094, 0.6701461988304094, -0.4398538011695907, 0.1151461988304093,
-0.0526315789473684, 0.0895906432748537, 0.0895906432748537, -0.2637426900584795, -0.0237426900584794, -0.0237426900584795, -0.2637426900584795, -0.2190204678362574, -0.2190204678362574, 0.2159795321637427, 0.2159795321637427, 0.2159795321637427, 0.2159795321637427, -0.2190204678362573, -0.2190204678362573, 0.1151461988304094, -0.4398538011695907, 0.6701461988304094, 0.1151461988304093,
-0.0526315789473685, 0.0895906432748538, 0.0895906432748538, -0.0237426900584795, -0.2637426900584796, -0.0237426900584796, -0.2637426900584796, 0.2159795321637428, 0.2159795321637427, -0.2190204678362573, -0.2190204678362573, 0.2159795321637427, 0.2159795321637427, -0.2190204678362573, -0.2190204678362573, -0.4398538011695907, 0.1151461988304094, 0.1151461988304093, 0.6701461988304094
		};

		T1 = new double[nf * nf]{
		0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
		0.0000000000000000, 0.1111111111111112, 0.1111111111111111, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, 0.0555555555555555, 0.0555555555555555, 0.0555555555555556, 0.0555555555555556, 0.0555555555555555, 0.0555555555555556, 0.0555555555555555, 0.0555555555555556, -0.1111111111111111, -0.1111111111111112, -0.1111111111111111, -0.1111111111111112,
		0.0000000000000000, 0.1111111111111111, 0.1111111111111112, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, 0.0555555555555555, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111,
		0.0000000000000000, -0.0555555555555556, -0.0555555555555556, 0.1111111111111112, 0.1111111111111111, -0.0555555555555555, -0.0555555555555555, 0.0555555555555556, 0.0555555555555555, 0.0555555555555555, 0.0555555555555556, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111, 0.0555555555555555, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556,
		0.0000000000000000, -0.0555555555555555, -0.0555555555555555, 0.1111111111111111, 0.1111111111111112, -0.0555555555555556, -0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555555, 0.0555555555555555, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111, 0.0555555555555555, 0.0555555555555555, 0.0555555555555555, 0.0555555555555555,
		0.0000000000000000, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555, -0.0555555555555555, 0.1111111111111112, 0.1111111111111111, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111, 0.0555555555555555, 0.0555555555555555, 0.0555555555555555, 0.0555555555555555, 0.0555555555555555, 0.0555555555555556, 0.0555555555555556, 0.0555555555555555,
		0.0000000000000000, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555, -0.0555555555555555, 0.1111111111111111, 0.1111111111111112, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111, -0.1111111111111111, 0.0555555555555555, 0.0555555555555555, 0.0555555555555555, 0.0555555555555555, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556,
		0.0000000000000000, 0.0555555555555555, 0.0555555555555555, 0.0555555555555556, 0.0555555555555556, -0.1111111111111111, -0.1111111111111111, 0.3611111111111112, -0.1388888888888889, -0.1388888888888889, 0.3611111111111111, -0.0555555555555557, -0.0555555555555555, -0.0555555555555555, -0.0555555555555557, -0.0555555555555556, -0.0555555555555555, -0.0555555555555554, -0.0555555555555556,
		0.0000000000000000, 0.0555555555555555, 0.0555555555555555, 0.0555555555555556, 0.0555555555555556, -0.1111111111111111, -0.1111111111111111, -0.1388888888888889, 0.3611111111111112, 0.3611111111111111, -0.1388888888888888, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555, -0.0555555555555555, -0.0555555555555555, -0.0555555555555555, -0.0555555555555556,
		0.0000000000000000, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, -0.1111111111111111, -0.1111111111111111, -0.1388888888888889, 0.3611111111111111, 0.3611111111111111, -0.1388888888888889, -0.0555555555555555, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555, -0.0555555555555555, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555,
		0.0000000000000000, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555555, -0.1111111111111111, -0.1111111111111111, 0.3611111111111111, -0.1388888888888889, -0.1388888888888889, 0.3611111111111113, -0.0555555555555555, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555,
		0.0000000000000000, 0.0555555555555555, 0.0555555555555555, -0.1111111111111111, -0.1111111111111111, 0.0555555555555556, 0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, 0.3611111111111109, -0.1388888888888889, -0.1388888888888889, 0.3611111111111111, -0.0555555555555555, -0.0555555555555555, -0.0555555555555555, -0.0555555555555554,
		0.0000000000000000, 0.0555555555555556, 0.0555555555555556, -0.1111111111111111, -0.1111111111111111, 0.0555555555555556, 0.0555555555555556, -0.0555555555555556, -0.0555555555555555, -0.0555555555555556, -0.0555555555555555, -0.1388888888888889, 0.3611111111111112, 0.3611111111111112, -0.1388888888888888, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556,
		0.0000000000000000, 0.0555555555555556, 0.0555555555555556, -0.1111111111111111, -0.1111111111111111, 0.0555555555555556, 0.0555555555555556, -0.0555555555555555, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.1388888888888888, 0.3611111111111111, 0.3611111111111112, -0.1388888888888888, -0.0555555555555556, -0.0555555555555555, -0.0555555555555555, -0.0555555555555556,
		0.0000000000000000, 0.0555555555555556, 0.0555555555555555, -0.1111111111111111, -0.1111111111111111, 0.0555555555555556, 0.0555555555555555, -0.0555555555555556, -0.0555555555555555, -0.0555555555555556, -0.0555555555555555, 0.3611111111111111, -0.1388888888888888, -0.1388888888888889, 0.3611111111111108, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555, -0.0555555555555556,
		0.0000000000000000, -0.1111111111111111, -0.1111111111111111, 0.0555555555555555, 0.0555555555555555, 0.0555555555555556, 0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, 0.3611111111111113, -0.1388888888888890, -0.1388888888888890, 0.3611111111111112,
		0.0000000000000000, -0.1111111111111111, -0.1111111111111111, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555555, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.1388888888888889, 0.3611111111111112, 0.3611111111111112, -0.1388888888888890,
		0.0000000000000000, -0.1111111111111111, -0.1111111111111111, 0.0555555555555556, 0.0555555555555555, 0.0555555555555556, 0.0555555555555555, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555, -0.0555555555555556, -0.0555555555555555, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.1388888888888890, 0.3611111111111111, 0.3611111111111112, -0.1388888888888889,
		0.0000000000000000, -0.1111111111111111, -0.1111111111111111, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, 0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555555, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, 0.3611111111111111, -0.1388888888888888, -0.1388888888888889, 0.3611111111111112
		};


	}
}

void BlockFreeSurface::setMRTdiagonal(double tauLoc) {
	if (nf == 15) { //D3Q15
		double s = 1 / tauLoc;
		S[9] = s;
		S[10] = s;
		S[11] = s;
		S[12] = s;
		S[13] = s;
	}

	if (nf == 19) { //D3Q19
		double s = 1 / tauLoc;
		S[9] = s;
		S[11] = s;
		S[13] = s;
		S[14] = s;
		S[15] = s;
	}
}

double BlockFreeSurface::nonNewtonianViscosityModel(double gamma) {
	double* b = fluidModelParam;
	double gammaC, lg;
	switch (viscosityModel) {
	case 1: // power-law
		return b[0] * pow(gamma, b[1] - 1.0);
	case 2: // bingham
		return b[0] + b[1] * pow(gamma, b[2] - 1.0);
	case 3: // cross
		return b[0] + (b[1] - b[0]) / (1 + b[2] * pow(gamma, b[3]));
	case 4: // poly
		return b[0] + b[1] * gamma + b[2] * gamma * gamma + b[3] * gamma * gamma * gamma;
	case 5: // modific power-law
		gammaC = exp((1 - b[1]) / (2 * b[2]));
		if (gamma < gammaC) {
			gamma = gammaC;
		}
		return b[0] * pow(gamma, b[1] - 1.0 + b[2] * log(gamma));
	case 6: // modific power-law 2
		gammaC = exp((-2 * b[2] - sqrt(4 * b[2] * b[2] - 12 * b[3] * (b[1] - 1))) / (6 * b[3]));
		if (gamma < gammaC) {
			gamma = gammaC;
		}
		lg = log(gamma);
		return b[0] * pow(gamma, b[1] - 1.0 + b[2] * lg + b[3] * lg * lg);
	default:
		throw std::invalid_argument(" unknown viscosity model! ");
	}
}
