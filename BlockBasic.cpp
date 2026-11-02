/*
 * File:   BlockBasic.cpp
 * Author: obublik
 *
 * Created on 2. února 2018, 12:27
 */

#include "BlockBasic.h"
#include <stdexcept>

using namespace std;

BlockBasic::BlockBasic(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, bool isTauVariable, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par) {
	initCommon(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
}

BlockBasic::~BlockBasic() {
}

// density and velocity of a cell, the velocity includes the half force correction (Guo et al. 2002)
void BlockBasic::moments(int ind, double& rho, double* u) {
	double d = 0, jx = 0, jy = 0, jz = 0;
	const pop_t* fp = &f[nf * ind];
	for (int p = 0; p < nf; p++) {
		double fVal = fp[p];
		d += fVal;
		jx += ex[p] * fVal;
		jy += ey[p] * fVal;
		jz += ez[p] * fVal;
	}
	rho = d + rho0;
	u[0] = jx / rho + 0.5 * G[0];
	u[1] = jy / rho + 0.5 * G[1];
	u[2] = jz / rho + 0.5 * G[2];
}

void BlockBasic::collision() {
	double feq[27], fneq[27], delta[27], source[27];
	double u[3];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind]) {
			continue;
		}
		double rho;
		moments(ind, rho, u);
		Rho[ind] = rho;
		dRho[ind] = rho - rho0;
		U[ind] = u[0];
		V[ind] = u[1];
		W[ind] = u[2];
		equilibrium(feq, rho, rho - rho0, u[0], u[1], u[2]);
		double F[3] = { rho * G[0], rho * G[1], rho * G[2] };
		pop_t* fp = &f[nf * ind];
		for (int p = 0; p < nf; p++) {
			double eu = ex[p] * u[0] + ey[p] * u[1] + ez[p] * u[2];
			double eF = ex[p] * F[0] + ey[p] * F[1] + ez[p] * F[2];
			double uF = u[0] * F[0] + u[1] * F[1] + u[2] * F[2];
			source[p] = w[p] * (3.0 * (eF - uF) + 9.0 * eu * eF);
			fneq[p] = fp[p] - feq[p] + 0.5 * source[p];
		}
		// f* = f - Omega (f - feq) + (I - Omega/2) S
		double tauLoc = tau;
		if (smagorinsky > 0) {
			double fn[27];
			for (int p = 0; p < nf; p++) {
				fn[p] = fp[p] - feq[p];
			}
			tauLoc = lesTau(tau, fn, rho);
		}
		coll.relax(fneq, tauLoc, delta);
		pop_t* fcp = &fc[nf * ind];
		for (int p = 0; p < nf; p++) {
			fcp[p] = (pop_t)(fp[p] - delta[p] + source[p]);
		}
	}
}

void BlockBasic::propagation() {
	for (int i = 0; i < nIn; i++) {
		int ind = indIn[i];
		if (!boundary[ind]) {
			for (int p = 0; p < nf; p++) {
				int indShift = ind + offs[p];
				int ct = cellType[indShift];
				if (ct == CT_FLUID) {
					f[nf * ind + opp[p]] = fc[nf * indShift + opp[p]];
				}
				else {
					f[nf * ind + opp[p]] = (pop_t)boundaryIncoming(ind, p, ct);
				}
			}
		}
	}

	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		if (!boundary[ind]) {
			for (int p = 0; p < nf; p++) {
				int bl, loc;
				int ct = neighbourType(ind, p, bl, loc);
				if (ct == CT_FLUID) {
					f[nf * ind + opp[p]] = blocks[bl]->fc[nf * loc + opp[p]];
				}
				else {
					f[nf * ind + opp[p]] = (pop_t)boundaryIncoming(ind, p, ct);
				}
			}
		}
	}
}

void BlockBasic::propagateCells(const vector<int>& cells) {
	for (int ind : cells) {
		if (boundary[ind]) {
			continue;
		}
		for (int p = 0; p < nf; p++) {
			int bl, loc;
			int ct = neighbourType(ind, p, bl, loc);
			if (ct == CT_FLUID) {
				f[nf * ind + opp[p]] = blocks[bl]->fc[nf * loc + opp[p]];
			}
			else {
				f[nf * ind + opp[p]] = (pop_t)boundaryIncoming(ind, p, ct);
			}
		}
	}
}

void BlockBasic::makro() {
	double u[3];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (!boundary[ind]) {
			double rho;
			moments(ind, rho, u);
			Rho[ind] = rho;
			dRho[ind] = rho - rho0;
			U[ind] = u[0];
			V[ind] = u[1];
			W[ind] = u[2];
		}
		else {
			Rho[ind] = 1;
			dRho[ind] = Rho[ind] - rho0;
			U[ind] = 0;
			V[ind] = 0;
			W[ind] = 0;
		}
	}
}

void BlockBasic::initZeroCondition() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		Rho[ind] = 1;
		dRho[ind] = Rho[ind] - rho0;
		U[ind] = 0;
		V[ind] = 0;
		W[ind] = 0;
	}
	computeInitDistribution();
}

void BlockBasic::readInitCondition(int id, string simulPath) {
	string fileName = simulPath + "init/blockInit" + fileTag + ".txt";
	ifstream myfile(fileName);
	if (!myfile.is_open()) {
		return; // zero condition stays
	}
	string line;
	try {
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (!getline(myfile, line)) {
				throw runtime_error("unexpected end of file");
			}
			istringstream iss(line);
			vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
			if (tokens.size() < 4) {
				throw runtime_error("expected 4 values on a line");
			}
			Rho[ind] = stod(tokens[0]);
			U[ind] = stod(tokens[1]) / Cu;
			V[ind] = stod(tokens[2]) / Cu;
			W[ind] = stod(tokens[3]) / Cu;
			dRho[ind] = Rho[ind] - rho0;
		}
	}
	catch (const std::exception& e) {
		throw runtime_error("Error in " + fileName + ": " + e.what());
	}
}

void BlockBasic::computeInitDistribution() {
	double feq[27];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind]);
		for (int p = 0; p < nf; p++) {
			fc[nf * ind + p] = (pop_t)feq[p];
			f[nf * ind + p] = (pop_t)feq[p];
		}
	}
}

void BlockBasic::writeBlock(int id) {
	ofstream myfile;
	myfile.open("results/results" + std::to_string(id) + ".txt", ios::out);
	myfile << ni << " " << nj << " " << nk << "\n";
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			for (int k = 0; k < nk; k++) {
				int ind = (ni * nj) * k + ni * j + i;
				myfile << Rho[ind] << " " << U[ind] << " " << V[ind] << " " << W[ind] << " " << boundary[ind] << "\n";
			}
		}
	}
	myfile.close();
}
