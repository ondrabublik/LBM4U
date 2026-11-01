/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "LagrangeParticlesVar2.h"
#include "Block.h"
#include "BlockFreeSurfaceLagrangeParticles.h"

using namespace std;

LagrangeParticlesVar2::LagrangeParticlesVar2(Block** blocks, int nBlock, int nBlockMax, int nBi, int nBj, int nBk, int ni, int nj, int nk, int nf, unordered_map<string, string> par) {
	this->blocks = blocks;
	this->nBlock = nBlock;
	this->nBlockMax = nBlockMax;
	this->nBi = nBi;
	this->nBj = nBj;
	this->nBk = nBk;
	this->ni = ni;
	this->nj = nj;
	this->nk = nk;
	this->nf = nf;

	Cu = stod(par["uMax"]) / 0.1;
	Ch = stod(par["H"]) / (this->ni);
	Ct = Ch / Cu;
	Cf = Ch / (Ct * Ct);
	CRho = stod(par["density"]);
	G = new double[3]{
		stod(par["gx"]) / Cf, stod(par["gy"]) / Cf, stod(par["gz"]) / Cf
	}; //gravity
	Cvis = Ch * Ch / Ct;
	vLB = stod(par["vKin"]) / Cvis;
}

LagrangeParticlesVar2::~LagrangeParticlesVar2() {
	delete[] Gx;
	delete[] Gy;
	delete[] Gz;

	delete[] partX;
	delete[] partY;
	delete[] partZ;

	delete[] veloX;
	delete[] veloY;
	delete[] veloZ;

	delete[] veloX0;
	delete[] veloY0;
	delete[] veloZ0;
}

void LagrangeParticlesVar2::newParticlePosition() {
#pragma omp parallel for
	for (int i = 0; i < nParticles; i++) {
		veloX[i] += Gx[i] / partMass[i] + (1 - 1 / partRho[i]) * G[0] + FinterX[i] / partMass[i];
		veloY[i] += Gy[i] / partMass[i] + (1 - 1 / partRho[i]) * G[1] + FinterY[i] / partMass[i];
		veloZ[i] += Gz[i] / partMass[i] + (1 - 1 / partRho[i]) * G[2] + FinterZ[i] / partMass[i];

		partX[i] += veloX[i];
		partY[i] += veloY[i];
		partZ[i] += veloZ[i];

		Gx[i] = 0;
		Gy[i] = 0;
		Gz[i] = 0;

		FinterX[i] = 0;
		FinterY[i] = 0;
		FinterZ[i] = 0;
	}
}

void LagrangeParticlesVar2::getParticleForce() {

#pragma omp parallel for
	for (int p = 0; p < nParticles; p++) {
		// compute block indexes of particle position
		int iB = (int)(partX[p] / ni);
		int jB = (int)(partY[p] / nj);
		int kB = (int)(partZ[p] / nk);
		int bl = nBi * nBj * kB + nBi * jB + iB;
		if (bl > -1 && bl < nBlockMax && blocks[bl]) { // particle is inside of the block
			int nf = blocks[bl]->nf;
			int* opp = blocks[bl]->opp;
			double* w = blocks[bl]->w;
			int* ex = blocks[bl]->ex;
			int* ey = blocks[bl]->ey;
			int* ez = blocks[bl]->ez;
			int* flagMap = ((BlockFreeSurface*)blocks[bl])->flagMap;
			double* f = blocks[bl]->f;
			bool* boundary = blocks[bl]->boundary;

			// compute local particle coordinates on block
			double x = (partX[p] - iB * ni);
			double y = (partY[p] - jB * nj);
			double z = (partZ[p] - kB * nk);

			// compute local indexes in block of nearest cell
			int iloc = (int)x;
			int jloc = (int)y;
			int kloc = (int)z;

			int indCenter = (ni * nj) * kloc + ni * jloc + iloc;

			// force from fluid ======================================================================================
			if ((flagMap[indCenter] == 1 || flagMap[indCenter] == 3) && !boundary[indCenter]) {
				int approxSize = 1;
				int iStart = iloc - approxSize;
				int iEnd = iloc + approxSize;
				int jStart = jloc - approxSize;
				int jEnd = jloc + approxSize;
				int kStart = kloc - approxSize;
				int kEnd = kloc + approxSize;

				// define Lagrange basis
				double* L = new double[(2 * approxSize + 1) * (2 * approxSize + 1) * (2 * approxSize + 1)];
				L[0] = 1;
				int s = 0;
				for (int i = iStart; i < iEnd + 1; i++) {
					for (int j = jStart; j < jEnd + 1; j++) {
						for (int k = kStart; k < kEnd + 1; k++) {
							double Lx = 1;
							for (int l = iStart; l < iEnd + 1; l++) {
								if (l != i)
									Lx *= (x - (l + 0.5)) / (i - l);
							}
							double Ly = 1;
							for (int m = jStart; m < jEnd + 1; m++) {
								if (m != j)
									Ly *= (y - (m + 0.5)) / (j - m);
							}
							double Lz = 1;
							for (int n = kStart; n < kEnd + 1; n++) {
								if (n != k)
									Lz *= (z - (n + 0.5)) / (k - n);
							}
							L[s] = Lx * Ly * Lz;
							s++;
						}
					}
				}

				// interpolation
				double fa[27];
				for (int p = 0; p < nf; p++) {
					fa[p] = 0;
				}
				s = 0;
				for (int i = iStart; i < iEnd + 1; i++) {
					for (int j = jStart; j < jEnd + 1; j++) {
						for (int k = kStart; k < kEnd + 1; k++) {
							bool gasForce = false;
							int ind = (ni * nj) * k + ni * j + i;
							if (i > -1 && i < ni && j > -1 && j < nj && k >-1 && k < nk) { // inside block
								if (!boundary[ind]) { // cell is not boundary
									if (flagMap[ind] == 1 || flagMap[ind] == 3) {
										for (int r = 0; r < nf; r++) {
											int pind = nf * ind + r;
											fa[r] += L[s] * f[pind];
										}
									}
									else { // point is outside of fluid
										for (int r = 0; r < nf; r++) {
											int pind = nf * indCenter + r;
											fa[r] += L[s] * f[pind];
										}
										gasForce = true;
									}
								}
								else { // point is outside of geometry
									for (int r = 0; r < nf; r++) {
										int pind = nf * indCenter + r;
										fa[r] += L[s] * f[pind];
									}
								}
							}
							else { // outside block
								int indBlockShift = 0;
								int indShift = 0;
								blocks[bl]->local2globalShift(bl, ind, indBlockShift, indShift, i - iStart - approxSize, j - jStart - approxSize, k - kStart - approxSize);
								if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
									int fm = ((BlockFreeSurface*)blocks[indBlockShift])->flagMap[indShift];
									double* fs = ((BlockFreeSurface*)blocks[indBlockShift])->f;
									bool* boundarys = ((BlockFreeSurface*)blocks[indBlockShift])->boundary;
									if (!boundarys[indShift]) { // cell is not boundary
										if (fm == 1 || fm == 3) {
											for (int r = 0; r < nf; r++) {
												int pind = nf * indShift + r;
												fa[r] += L[s] * fs[pind];
											}
										}
										else {
											for (int r = 0; r < nf; r++) {
												int pind = nf * indCenter + r;
												fa[r] += L[s] * f[pind];
											}
											gasForce = true;
										}
									}
									else {
										for (int r = 0; r < nf; r++) {
											int pind = nf * indCenter + r;
											fa[r] += L[s] * f[pind];
										}
									}
								}
								else {
									for (int r = 0; r < nf; r++) {
										int pind = nf * indCenter + r;
										fa[r] += L[s] * f[pind];
									}
								}
							}

							if (gasForce) {

							}
							s++;
						}
					}
				}
				delete[] L;

				double rho = 0;
				double ulb = 0;
				double vlb = 0;
				double wlb = 0;
				for (int r = 0; r < nf; r++) {
				    rho += fa[r];
					ulb += fa[r] * ex[r];
					vlb += fa[r] * ey[r];
					wlb += fa[r] * ez[r];
				}
				double Urel = sqrt((veloX[p] - ulb) * (veloX[p] - ulb) + (veloY[p] - vlb) * (veloY[p] - vlb) + (veloZ[p] - wlb) * (veloZ[p] - vlb));
				//double Up = sqrt(veloX[p]* veloX[p] + veloY[p] * veloY[p] + veloZ[p] * veloZ[p]);
				double Re = partD[p] * Urel / vLB;
				double Cd = cd(Re);

				for (int r = 0; r < nf; r++) {
					int op = opp[r];
					double fnew = fa[op] - 6 * w[op] * rho * (ex[op] * veloX[p] + ey[op] * veloY[p] + ez[op] * veloZ[p]);
					Gx[p] += ex[op] * (fnew - fa[r]) * partCrossection[p] * Cd*Urel/4;
					Gy[p] += ey[op] * (fnew - fa[r]) * partCrossection[p] * Cd*Urel/4;
					Gz[p] += ez[op] * (fnew - fa[r]) * partCrossection[p] * Cd*Urel/4;
				}

				/*
				for (int r = 0; r < nf; r++) {
					int op = opp[r];
					Gx[p] += ((ex[r] - veloX[p]) * fa[r] - (ex[op] - veloX[p]) * fa[op]) * partVolume[p];
					Gy[p] += ((ey[r] - veloY[p]) * fa[r] - (ey[op] - veloY[p]) * fa[op]) * partVolume[p];
					Gz[p] += ((ez[r] - veloZ[p]) * fa[r] - (ez[op] - veloZ[p]) * fa[op]) * partVolume[p];
				}*/

				/*
				for (int r = 0; r < nf; r++) {
					double fnew = 2*fa[r] - 6 * w[r] * rho * (ex[r] * veloX[p] + ey[r] * veloY[p] + ez[r] * veloZ[p]);
					Gx[p] += ex[r] * fnew * partVolume[p];
					Gy[p] += ey[r] * fnew * partVolume[p];
					Gz[p] += ez[r] * fnew * partVolume[p];
				}*/

			}
			else { // particle is out of fluid
				//cout << "outside fluid" << endl;
				Gx[p] = 0;
				Gy[p] = 0;
				Gz[p] = 0;
			}

			// force from wall ====================================================================================================================
			double wallDist = blocks[bl]->wallDistance[indCenter] + blocks[bl]->wallDistanceGradient[3 * indCenter] * (x - (iloc + 0.5)) + blocks[bl]->wallDistanceGradient[3 * indCenter + 1] * (y - (jloc + 0.5)) + blocks[bl]->wallDistanceGradient[3 * indCenter + 2] * (z - (kloc + 0.5));
			if (wallDist < 2) {
				double fWall = gammaWall(wallDist);
				FinterX[p] += fWall * blocks[bl]->wallDistanceGradient[3 * indCenter] * partMass[p];
				FinterY[p] += fWall * blocks[bl]->wallDistanceGradient[3 * indCenter + 1] * partMass[p];
				FinterZ[p] += fWall * blocks[bl]->wallDistanceGradient[3 * indCenter + 2] * partMass[p];
			}
		}
		else { // particle is outside of the block
			cout << "outside block" << endl;
			Gx[p] = 0;
			Gy[p] = 0;
			Gz[p] = 0;
		}
	}
}

void LagrangeParticlesVar2::applyForce() {
	for (int p = 0; p < nParticles; p++) {
		int iB = (int)(partX[p] / ni);
		int jB = (int)(partY[p] / nj);
		int kB = (int)(partZ[p] / nk);
		int bl = nBi * nBj * kB + nBi * jB + iB;
		double x = (partX[p] - iB * ni);
		double y = (partY[p] - jB * nj);
		double z = (partZ[p] - kB * nk);
		int iloc = (int)x;
		int jloc = (int)y;
		int kloc = (int)z;

		if (bl > -1 && bl < nBlockMax && blocks[bl]) {
			// force from fluid
			double* particleFx = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->particleForceX;
			double* particleFy = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->particleForceY;
			double* particleFz = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->particleForceZ;
			int* flagMap = ((BlockFreeSurface*)blocks[bl])->flagMap;
			bool* boundary = blocks[bl]->boundary;

			int indCenter = (ni * nj) * kloc + ni * jloc + iloc;
			if ((flagMap[indCenter] == 1 || flagMap[indCenter] == 3) && !boundary[indCenter]) {
				int interpLevelForce = 2;
				for (int i = -interpLevelForce; i < interpLevelForce + 1; i++) {
					for (int j = -interpLevelForce; j < interpLevelForce + 1; j++) {
						for (int k = -interpLevelForce; k < interpLevelForce + 1; k++) {
							int ip = iloc + i;
							int jp = jloc + j;
							int kp = kloc + k;
							if (ip > -1 && ip < ni && jp > -1 && jp < nj && kp > -1 && kp < nk) { // inside block
								int ind = (ni * nj) * kp + ni * jp + ip;
								if (flagMap[ind] == 1 || flagMap[ind] == 3) {
									double weight = delta(iloc + 0.5 - x + i) * delta(jloc + 0.5 - y + j) * delta(kloc + 0.5 - z + k);
									particleFx[ind] -= weight * Gx[p];
									particleFy[ind] -= weight * Gy[p];
									particleFz[ind] -= weight * Gz[p];
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
										((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->particleForceX[indShift] -= weight * Gx[p];
										((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->particleForceY[indShift] -= weight * Gy[p];
										((BlockFreeSurfaceLagrangeParticles*)blocks[indBlockShift])->particleForceZ[indShift] -= weight * Gz[p];

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

double LagrangeParticlesVar2::delta(double x) {
	x = abs(x);
	double y = 0;
	if (x <= 2) {
		y = (1 + cos(3.1415 * x / 2.0)) / 4.0;
	}
	return y;
}

double LagrangeParticlesVar2::gammaPP(double gap, double avgD) {
	double y = 0;
	double L = 2 * avgD;
	if (abs(gap) <= L) {
		y = 1e-3 * sin(3.1415 * gap / L);
	}
	return y;
}

double LagrangeParticlesVar2::gammaWall(double x) {
	x = abs(x);
	double y = 0;
	if (x <= 1.5) {
		y = 1 * (1 - x / 1.5) * (1 - x / 1.5);
	}
	return y;
}

double LagrangeParticlesVar2::cd(double Re) {
	double Cd = 0.44;
	Re = max(Re, 0.3);
	if (Re < 1000) {
		if (Re > 0.5) {
			Cd = 24.0 / Re * (1 + 0.15 * pow(Re, 0.687));
		}
		else {
			Cd = 24.0 / Re;
		}
	}
	return Cd;
}

void LagrangeParticlesVar2::particleParticleInteraction() {
	// put particles into cells box
	for (int p = 0; p < nParticles; p++) {
		int iB = (int)(partX[p] / ni);
		int jB = (int)(partY[p] / nj);
		int kB = (int)(partZ[p] / nk);
		int bl = nBi * nBj * kB + nBi * jB + iB;
		if (bl > -1 && bl < nBlockMax && blocks[bl]) { // particle is inside of the block
			int nParticleBox = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->nParticleBox;
			int* particleBox = ((BlockFreeSurfaceLagrangeParticles*)blocks[bl])->particleBox;

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
			}
			else {
				cout << "out of particle box" << endl;
			}
		}
	}

	// get particles from surrounding boxes
#pragma omp parallel for
	for (int p = 0; p < nParticles; p++) {
		int iB = (int)partX[p] / ni;
		int jB = (int)partY[p] / nj;
		int kB = (int)partZ[p] / nk;
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
								double sx = partX[p] - partX[q];
								double sy = partY[p] - partY[q];
								double sz = partZ[p] - partZ[q];
								double L = sqrt(sx * sx + sy * sy + sz * sz);
								if (L > 1e-10) {
									sx /= L;
									sy /= L;
									sz /= L;
								}
								else {
									sx = 0;
									sy = 0;
									sz = 0;
								}
								double avgD = (partD[p] + partD[q]) / 2;
								double gap = L - avgD;
								double f = gammaPP(gap, avgD);
								double fx = f * sx;
								double fy = f * sy;
								double fz = f * sz;

								FinterX[p] += fx;
								FinterY[p] += fy;
								FinterZ[p] += fz;
							}
						}
					}
				}
			}
		}
	}
}

void LagrangeParticlesVar2::readInitParticlesDistribution(string simulPath, double Ch, double Cu) {
	string line;
	string fileName = simulPath + "init/particles.txt";
	int s = 0;
	try {
		ifstream* myfile = new ifstream(fileName);
		if (myfile->is_open()) {
			getline(*myfile, line);
			nParticles = stod(line);
			partMass = new double[nParticles];
			partCrossection = new double[nParticles];
			partD = new double[nParticles];
			partRho = new double[nParticles];
			partX = new double[nParticles];
			partY = new double[nParticles];
			partZ = new double[nParticles];
			veloX = new double[nParticles];
			veloY = new double[nParticles];
			veloZ = new double[nParticles];
			veloX0 = new double[nParticles];
			veloY0 = new double[nParticles];
			veloZ0 = new double[nParticles];
			Gx = new double[nParticles];
			Gy = new double[nParticles];
			Gz = new double[nParticles];
			FinterX = new double[nParticles];
			FinterY = new double[nParticles];
			FinterZ = new double[nParticles];

			for (int i = 0; i < nParticles; i++) {
				getline(*myfile, line);
				istringstream iss(line);
				vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
				partX[s] = stod(tokens[0]) / Ch;
				partY[s] = stod(tokens[1]) / Ch;
				partZ[s] = stod(tokens[2]) / Ch;
				partD[s] = stod(tokens[3]) / Ch;
				partRho[s] = stod(tokens[4]) / CRho;
				partCrossection[s] = 3.1415 / 4 * partD[s] * partD[s];
				partMass[s] = 3.1415 / 6 * partD[s] * partD[s] * partD[s] * partRho[s];
				Gx[s] = 0;
				Gy[s] = 0;
				Gz[s] = 0;
				FinterX[s] = 0;
				FinterY[s] = 0;
				FinterZ[s] = 0;
				veloX[s] = 0;
				veloY[s] = 0;
				veloZ[s] = 0;
				veloX0[s] = 0;
				veloY0[s] = 0;
				veloZ0[s] = 0;
				s++;
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
	for (int i = 0; i < nParticles; i++) {
		veloX[i] = 0;
		veloY[i] = 0;
		veloZ[i] = 0;
	}
}

void LagrangeParticlesVar2::saveParticlesPositionsToVTK(string simulPath, int iter) {
	ofstream myfile;
	myfile.open(simulPath + "result/output_particles" + std::to_string((int)(iter * Ct * 1000000)) + "us.vtk", ios::out);
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
		myfile << 1 << "\n";
	}
	myfile << "\n";

	myfile << "CELL_DATA " << nParticles << "\n";
	myfile << "SCALARS densityRatio float 1\n";
	myfile << "LOOKUP_TABLE default\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << 0 << "\n";
	}
	myfile << "\n";

	myfile << "VECTORS rychlost float\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << Cu * veloX[i] << " " << Cu * veloY[i] << " " << Cu * veloZ[i] << "\n";
	}
	myfile << "\n";

	myfile << "VECTORS force float\n";
	for (int i = 0; i < nParticles; i++) {
		myfile << Gx[i] * Cf << " " << Gy[i] * Cf << " " << Gz[i] * Cf << "\n";
	}
	myfile << "\n";
	myfile.close();
}

