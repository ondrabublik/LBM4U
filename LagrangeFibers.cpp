/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "LagrangeFibers.h"
#include "Block.h"
#include "BlockFreeSurfaceLagrangeFibers.h"

using namespace std;

LagrangeFibers::LagrangeFibers(Block** blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, unordered_map<string, string> par) {
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

	Cvis = Ch * Ch / Ct;
	vLB = b.getParDouble("vKin", par, true, 0) / Cvis;

	interpLevel = 2;
	nChain = 10;
}

LagrangeFibers::~LagrangeFibers() {
	delete fibreCoords;
	delete fibreVelo;
}

void LagrangeFibers::newFiberPosition(double dtP) {
#pragma omp parallel for
	for (int i = 0; i < nFibers; i++) {
		for (int c = 0; c < nChain; c++) {
			int px = 3 * nChain * i + 3 * c; // fibre chain index
			int py = 3 * nChain * i + 3 * c + 1; // fibre chain index
			int pz = 3 * nChain * i + 3 * c + 2; // fibre chain index
			double wallDist = wallDistance(px, fibreCoords[px], fibreCoords[py], fibreCoords[pz]);
			if (!inGas[nChain * i + c]) {
				// fluid velocity correction on wall
				if (wallDist < 1) {
					double norm = partWallNormal[px] * fibreVelo[px] + partWallNormal[py] * fibreVelo[py] + partWallNormal[pz] * fibreVelo[pz];
					if (norm < 0) {
						fibreFluidVelo[px] -= norm * partWallNormal[px];
						fibreFluidVelo[py] -= norm * partWallNormal[py];
						fibreFluidVelo[pz] -= norm * partWallNormal[pz];
					}
				}

				// forces in fibers
				double cfib = 0.001;
				double fibx = 0;
				double fiby = 0;
				double fibz = 0;
				if (c > 0) {
					double rx = fibreCoords[px-3] - fibreCoords[px];
					double ry = fibreCoords[py-3] - fibreCoords[py];
					double rz = fibreCoords[pz-3] - fibreCoords[pz];
					double r = sqrt(rx * rx + ry * ry + rz * rz);
					if (r < 1e-3) {
						r = 1e-3;
					}
					fibx += cfib * rx / r * (r - rChain);
					fiby += cfib * ry / r * (r - rChain);
					fibz += cfib * rz / r * (r - rChain);
				}
				if (c < nChain - 1) {
					double rx = fibreCoords[px + 3] - fibreCoords[px];
					double ry = fibreCoords[py + 3] - fibreCoords[py];
					double rz = fibreCoords[pz + 3] - fibreCoords[pz];
					double r = sqrt(rx * rx + ry * ry + rz * rz);
					if (r < 1e-3) {
						r = 1e-3;
					}
					fibx += cfib * rx / r * (r - rChain);
					fiby += cfib * ry / r * (r - rChain);
					fibz += cfib * rz / r * (r - rChain);
				}

				// compute interaction force
				double Uslip = sqrt((fibreVelo[px] - fibreFluidVelo[px]) * (fibreVelo[px] - fibreFluidVelo[px]) + (fibreVelo[py] - fibreFluidVelo[py]) * (fibreVelo[py] - fibreFluidVelo[py]) + (fibreVelo[pz] - fibreFluidVelo[pz]) * (fibreVelo[pz] - fibreFluidVelo[pz]));
				double Cd_mult_Uslip = cd(partD[i], Uslip);
				double coef = 0.5 * Cd_mult_Uslip * partCrossection[i];

				fibreVelo[px] = (fibreVelo[px] * partMass[i] / dtP + coef * fibreFluidVelo[px] + G[0] * (partMass[i] - partVolume[i]) + 0 * Fwall[px] + fibx) / (partMass[i] / dtP + coef);
				fibreVelo[py] = (fibreVelo[py] * partMass[i] / dtP + coef * fibreFluidVelo[py] + G[1] * (partMass[i] - partVolume[i]) + 0 * Fwall[py] + fiby) / (partMass[i] / dtP + coef);
				fibreVelo[pz] = (fibreVelo[pz] * partMass[i] / dtP + coef * fibreFluidVelo[pz] + G[2] * (partMass[i] - partVolume[i]) + 0 * Fwall[pz] + fibz) / (partMass[i] / dtP + coef);

				fibreFluidForce[px] = coef * (fibreVelo[px] - fibreFluidVelo[px]);
				fibreFluidForce[py] = coef * (fibreVelo[py] - fibreFluidVelo[py]);
				fibreFluidForce[pz] = coef * (fibreVelo[pz] - fibreFluidVelo[pz]);
			}
			else {
				fibreVelo[px] += dtP * G[0];
				fibreVelo[py] += dtP * G[1];
				fibreVelo[pz] += dtP * G[2];
			}

			fibreCoords[px] += dtP * fibreVelo[px];
			fibreCoords[py] += dtP * fibreVelo[py];
			fibreCoords[pz] += dtP * fibreVelo[pz];

			// test of the next possition of the fibre chain
			if (wallDist < partD[i] / 2) {
				double vx0 = fibreVelo[px];
				double vy0 = fibreVelo[py];
				double vz0 = fibreVelo[pz];
				double x0 = fibreCoords[px] - dtP * vx0;
				double y0 = fibreCoords[py] - dtP * vy0;
				double z0 = fibreCoords[pz] - dtP * vz0;
				double norm = partWallNormal[px] * vx0 + partWallNormal[py] * vy0 + partWallNormal[pz] * vz0;

				fibreVelo[px] -= 2 * norm * partWallNormal[px];
				fibreVelo[py] -= 2 * norm * partWallNormal[py];
				fibreVelo[pz] -= 2 * norm * partWallNormal[pz];

				int n = 1;
				while (wallDist < partD[i] / 2) {
					fibreCoords[px] = x0 + dtP * (vx0 / n - norm * partWallNormal[px]);
					fibreCoords[py] = y0 + dtP * (vy0 / n - norm * partWallNormal[py]);
					fibreCoords[pz] = z0 + dtP * (vz0 / n - norm * partWallNormal[pz]);
					wallDist = wallDistance(i, fibreCoords[px], fibreCoords[py], fibreCoords[pz]);

					if (n > 4) {
						fibreCoords[px] = x0;
						fibreCoords[py] = y0;
						fibreCoords[pz] = z0;
						fibreVelo[px] = vx0;
						fibreVelo[py] = vy0;
						fibreVelo[pz] = vz0;
						break;
					}
					n *= 2;
				}
			}
		}
	}
}

void LagrangeFibers::interpolateFluidVelocity() {
#pragma omp parallel for
	for (int p = 0; p < nFibers; p++) {
		for (int c = 0; c < nChain; c++) {
			// save old velocities
			int px = 3 * nChain * p + 3 * c; // fibre chain index
			int py = 3 * nChain * p + 3 * c + 1; // fibre chain index
			int pz = 3 * nChain * p + 3 * c + 2; // fibre chain index
			fibreFluidVelo[px] = 0;
			fibreFluidVelo[py] = 0;
			fibreFluidVelo[pz] = 0;

			// compute block indexes of particle position
			double intX = fibreCoords[px] / ni;
			if (intX < 0) {
				intX -= 1;
			}
			double intY = fibreCoords[py] / nj;
			if (intY < 0) {
				intY -= 1;
			}
			double intZ = fibreCoords[pz] / nk;
			if (intZ < 0) {
				intZ -= 1;
			}
			int iB = (int)intX;
			int jB = (int)intY;
			int kB = (int)intZ;
			int bl = nBi * nBj * kB + nBi * jB + iB;
			if (bl > -1 && bl < nBlockMax && blocks[bl]) {
				int* flagMap = ((BlockFreeSurface*)blocks[bl])->flagMap;
				bool* boundary = ((BlockFreeSurface*)blocks[bl])->boundary;
				double* U = blocks[bl]->U;
				double* V = blocks[bl]->V;
				double* W = blocks[bl]->W;

				// compute local particle coordinates on block
				double x = (fibreCoords[px] - iB * ni);
				double y = (fibreCoords[py] - jB * nj);
				double z = (fibreCoords[pz] - kB * nk);

				// compute local indexes in block of nearest cell
				int iloc = (int)x;
				int jloc = (int)y;
				int kloc = (int)z;

				int ind = (ni * nj) * kloc + ni * jloc + iloc;
				if (!boundary[ind]) {
					inGas[nChain * p + c] = true;
					if (flagMap[ind] == 1 || flagMap[ind] == 3) {
						inGas[nChain * p + c] = false;

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
											fibreFluidVelo[px] += weight * U[ind];
											fibreFluidVelo[py] += weight * V[ind];
											fibreFluidVelo[pz] += weight * W[ind];
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
												fibreFluidVelo[px] += weight * ((BlockFreeSurface*)blocks[indBlockShift])->U[indShift];
												fibreFluidVelo[py] += weight * ((BlockFreeSurface*)blocks[indBlockShift])->V[indShift];
												fibreFluidVelo[pz] += weight * ((BlockFreeSurface*)blocks[indBlockShift])->W[indShift];
												sum += weight;
											}
										}
									}
								}
							}
						}
						fibreFluidVelo[px] /= sum;
						fibreFluidVelo[py] /= sum;
						fibreFluidVelo[pz] /= sum;
					}
				}
			}
		}
	}
}

double LagrangeFibers::wallDistance(int p, double xin, double yin, double zin) {

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

			partWallNormal[p + 0] = coeff[ic + 1] + coeff[ic + 4] * y + coeff[ic + 6] * z + coeff[ic + 7] * y * z;
			partWallNormal[p + 1] = coeff[ic + 2] + coeff[ic + 4] * x + coeff[ic + 5] * z + coeff[ic + 7] * x * z;
			partWallNormal[p + 2] = coeff[ic + 3] + coeff[ic + 5] * y + coeff[ic + 6] * x + coeff[ic + 7] * x * y;

			double norm = sqrt(partWallNormal[p] * partWallNormal[p] + partWallNormal[p + 1] * partWallNormal[p + 1] + partWallNormal[p + 2] * partWallNormal[p + 2]);

			if (norm > 0) {
				partWallNormal[p] = partWallNormal[p] / norm;
				partWallNormal[p + 1] = partWallNormal[p + 1] / norm;
				partWallNormal[p + 2] = partWallNormal[p + 2] / norm;
			}
		}
	}
	return L;
}

void LagrangeFibers::applyForce() {
	for (int p = 0; p < nFibers; p++) {
		for (int c = 0; c < nChain; c++) {
			int px = 3 * nChain * p + 3 * c; // fibre chain index
			int py = 3 * nChain * p + 3 * c + 1; // fibre chain index
			int pz = 3 * nChain * p + 3 * c + 2; // fibre chain index

			double intX = fibreCoords[px] / ni;
			if (intX < 0) {
				intX -= 1;
			}
			double intY = fibreCoords[py] / nj;
			if (intY < 0) {
				intY -= 1;
			}
			double intZ = fibreCoords[pz] / nk;
			if (intZ < 0) {
				intZ -= 1;
			}
			int iB = (int)intX;
			int jB = (int)intY;
			int kB = (int)intZ;
			int bl = nBi * nBj * kB + nBi * jB + iB;
			double x = (fibreCoords[px] - iB * ni);
			double y = (fibreCoords[py] - jB * nj);
			double z = (fibreCoords[pz] - kB * nk);
			int iloc = (int)x;
			int jloc = (int)y;
			int kloc = (int)z;

			double Fx = fibreFluidForce[px];
			double Fy = fibreFluidForce[py];
			double Fz = fibreFluidForce[pz];

			if (bl > -1 && bl < nBlockMax && blocks[bl]) {
				double* particleFx = ((BlockFreeSurfaceLagrangeFibers*)blocks[bl])->particleForceX;
				double* particleFy = ((BlockFreeSurfaceLagrangeFibers*)blocks[bl])->particleForceY;
				double* particleFz = ((BlockFreeSurfaceLagrangeFibers*)blocks[bl])->particleForceZ;
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
										((BlockFreeSurfaceLagrangeFibers*)blocks[indBlockShift])->particleForceX[indShift] += weight * Fx;
										((BlockFreeSurfaceLagrangeFibers*)blocks[indBlockShift])->particleForceY[indShift] += weight * Fy;
										((BlockFreeSurfaceLagrangeFibers*)blocks[indBlockShift])->particleForceZ[indShift] += weight * Fz;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


void LagrangeFibers::fiberWallInteraction() {
	for (int i = 0; i < nFibers; i++) {
		for (int c = 0; c < nChain; c++) {
			int px = 3 * nChain * i + 3 * c; // fibre chain index
			int py = 3 * nChain * i + 3 * c + 1; // fibre chain index
			int pz = 3 * nChain * i + 3 * c + 2; // fibre chain index
			double L = wallDistance(px, fibreCoords[px], fibreCoords[py], fibreCoords[pz]);
			double wallGap = partD[i];
			if (L > 0 && L < wallGap) {
				double xi = 1 - abs(L - partD[i] / 2) / (wallGap - partD[i] / 2);
				double Fw = 0.05 * xi;
				Fwall[px] = Fw * partWallNormal[px];
				Fwall[py] = Fw * partWallNormal[py];
				Fwall[pz] = Fw * partWallNormal[pz];
			}
			else {
				Fwall[px] = 0;
				Fwall[py] = 0;
				Fwall[pz] = 0;
			}
		}
	}
}

double LagrangeFibers::delta(double x) {
	x = abs(x);
	double y = 0;
	if (x <= 2) {
		y = (1 + cos(3.1415 * x / 2.0)) / 4.0;
	}
	return y;
}

double LagrangeFibers::cd(double d, double uSlip) { // return Cd*uSlip
	double CdUslip = 0.44 * uSlip;
	double Re = d * uSlip / vLB;
	if (Re < 1000) {
		CdUslip = 24.0 / d * vLB * (1 + 0.15 * pow(Re, 0.687));
	}
	return CdUslip;
}


void LagrangeFibers::readInitFibersDistribution(string simulPath, double Ch, double Cu) {
	string line;
	string fileName = simulPath + "init/particles.txt";
	try {
		ifstream* myfile = new ifstream(fileName);
		if (myfile->is_open()) {
			getline(*myfile, line);
			nFibers = stoi(line);
			partD = new double[nFibers];
			partRho = new double[nFibers];
			partMass = new double[nFibers];
			partVolume = new double[nFibers];
			partCrossection = new double[nFibers];
			int nfc = nFibers * nChain * 3;
			fibreCoords = new double[nfc];
			fibreVelo = new double[nfc];
			fibreFluidVelo = new double[nfc];
			fibreFluidForce = new double[nfc];
			Fwall = new double[nfc];
			volumeSolid = new double[nFibers];
			areaSolid = new double[nFibers];
			partWallNormal = new double[nfc];
			inGas = new bool[nFibers * nChain];
			for (int i = 0; i < nFibers; i++) {
				getline(*myfile, line);
				istringstream iss(line);
				vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
				double xc = stod(tokens[0]) / Ch;
				double yc = stod(tokens[1]) / Ch;
				double zc = stod(tokens[2]) / Ch;
				partD[i] = stod(tokens[3]) / Ch;
				partRho[i] = stod(tokens[4]) / Crho;
				partCrossection[i] = 3.1415 / 4 * partD[i] * partD[i];
				partVolume[i] = 3.1415 / 6 * partD[i] * partD[i] * partD[i];
				partMass[i] = partVolume[i] * partRho[i];

				double dirx = rand();
				double diry = rand();
				double dirz = rand();
				double ndir = sqrt(dirx * dirx + diry * diry + dirz * dirz);
				for (int c = 0; c < nChain; c++) {
					int px = 3 * nChain * i + 3 * c; // fibre chain index
					int py = 3 * nChain * i + 3 * c + 1; // fibre chain index
					int pz = 3 * nChain * i + 3 * c + 2; // fibre chain index
					fibreCoords[px] = xc + 0.5 * c * dirx / ndir;
					fibreCoords[py] = yc + 0.5 * c * diry / ndir;
					fibreCoords[pz] = zc + 0.5 * c * dirz / ndir;
					fibreVelo[px] = 0;
					fibreVelo[py] = 0;
					fibreVelo[pz] = 0;
					inGas[nChain * i + c] = false;
				}
			}
		}
		if (nChain > 1) {
			double dx = fibreCoords[3] - fibreCoords[0];
			double dy = fibreCoords[4] - fibreCoords[1];
			double dz = fibreCoords[5] - fibreCoords[2];
			rChain = sqrt(dx*dx + dy*dy + dz*dz);
		}
		myfile->close();
		delete myfile;
	}
	catch (int ex) {
		cout << "Error in " << fileName << endl;
	}
	cout << nFibers << " Lagrangian fibers were readed!" << endl;
}

void LagrangeFibers::saveFibersPositionsToVTK(string simulPath, int iter) {
	ofstream myfile;
	myfile.open(simulPath + "result/output_particles" + std::to_string((int)(iter * Ct * 1000000)) + "us.vtk", ios::out);
	myfile << "# vtk DataFile Version 2.6\noutput file\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	myfile << "POINTS " << nFibers * nChain << " float\n";
	for (int i = 0; i < nFibers; i++) {
		for (int c = 0; c < nChain; c++) {
			int px = 3 * nChain * i + 3 * c; // fibre chain index
			int py = 3 * nChain * i + 3 * c + 1; // fibre chain index
			int pz = 3 * nChain * i + 3 * c + 2; // fibre chain index
			myfile << fibreCoords[px] << " " << fibreCoords[py] << " " << fibreCoords[pz] << "\n";
		}
	}
	myfile << "\n";

	int s = 0;
	myfile << "CELLS " << nFibers << " " << nFibers * (nChain + 1) << "\n";
	for (int i = 0; i < nFibers; i++) {
		myfile << nChain;
		for (int c = 0; c < nChain; c++) {
			myfile << " " << s;
			s++;
		}
		myfile << "\n";
	}
	myfile << "\n";

	myfile << "CELL_TYPES " << nFibers << "\n";
	for (int i = 0; i < nFibers; i++) {
		myfile << "4" << "\n";
	}
	myfile << "\n";

	myfile << "POINT_DATA " << nFibers * nChain << "\n";
	myfile << "VECTORS rychlost float\n";
	for (int i = 0; i < nFibers; i++) {
		for (int c = 0; c < nChain; c++) {
			int px = 3 * nChain * i + 3 * c; // fibre chain index
			int py = 3 * nChain * i + 3 * c + 1; // fibre chain index
			int pz = 3 * nChain * i + 3 * c + 2; // fibre chain index
			myfile << fibreVelo[px] * Cu << " " << fibreVelo[py] * Cu << " " << fibreVelo[pz] * Cu << "\n";
		}
	}
	myfile << "\n";

	myfile.close();
}
