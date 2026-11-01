/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include "Block.h"

using namespace std;

void Block::local2globalShift(int indBlock, int ind, int& indBlockNeigh, int& indLocalNeigh, int ex, int ey, int ez) {
	int k = ind / (ni * nj) + ez;
	int pom = ind % (ni * nj);
	int j = pom / ni + ey;
	int i = pom % ni + ex;

	int kB = indBlock / (nBi * nBj);
	pom = indBlock % (nBi * nBj);
	int jB = pom / nBi;
	int iB = pom % nBi;

	if (i <= -1) {
		i = ni + i;
		iB--;
		if (xPeriodicity && iB == -1) {
			iB = nBi - 1;
		}
	}
	if (i >= ni) {
		i = i - ni;
		iB++;
		if (xPeriodicity && iB == nBi) {
			iB = 0;
		}
	}

	if (j <= -1) {
		j = nj + j;
		jB--;
		if (yPeriodicity && jB == -1) {
			jB = nBj - 1;
		}
	}

	if (j >= nj) {
		j = j - nj;
		jB++;
		if (yPeriodicity && jB == nBj) {
			jB = 0;
		}
	}

	if (k <= -1) {
		k = nk + k;
		kB--;
		if (zPeriodicity && kB == -1) {
			kB = nBk - 1;
		}
	}

	if (k >= nk) {
		k = k - nk;
		kB++;
		if (zPeriodicity && kB == nBk) {
			kB = 0;
		}
	}

	indBlockNeigh = nBi * nBj * kB + nBi * jB + iB;
	if (iB < 0 || iB == nBi || jB < 0 || jB == nBj || kB < 0 || kB == nBk) {
		indBlockNeigh = -1;
	}
	indLocalNeigh = ni * nj * k + ni * j + i;
}

void Block::computeShiftIndexes() {
	shift = new int[nAll * nf];
	bShift = new int[nAll * nf];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		for (int dir = 0; dir < nf; dir++) {
			local2globalShift(indBlock, ind, bShift[nAll * dir + ind], shift[nAll * dir + ind], ex[dir], ey[dir], ez[dir]);
		}
	}
}

void Block::initZeroCondition()
{
}

void Block::readInitCondition(int id, string simulPath)
{
}

void Block::computeInitDistribution()
{
}

void Block::initBoundary(int id, string geomPath) {
	// try to load boundary
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		boundary[ind] = false;
	}
	string line;
	string fileName = geomPath + "mesh/boundary/blockBoundary" + std::to_string(id) + ".txt";
	try {
		ifstream* myfile = new ifstream(fileName);
		if (myfile->is_open()) {
			for (int i = 0; i < nAll; i++) {
				int ind = indAll[i];
				getline(*myfile, line);
				boundary[ind] = (line.compare("1") == 1);
			}
		}
		myfile->close();
		delete myfile;
	}
	catch (int ex) {
		cout << "Error in " << fileName << endl;
	}
}

void Block::equilibrium(double* feq, double Rho, double dRho, double U, double V, double W) {
	for (int p = 0; p < nf; p++) { // vypocet equilibria
		double x = ex[p] * U + ey[p] * V + ez[p] * W;
		double u2 = U * U + V * V + W * W;
		feq[p] = w[p] * (Rho * (3.0 * x + 9.0 / 2.0 * x * x - 3.0 / 2.0 * u2) + dRho);
	}
}

void Block::equilibrium(double* feq, double Rho, double dRho, double U, double V, double W, double epsF) {
	for (int p = 0; p < nf; p++) { // vypocet equilibria
		double x = ex[p] * U + ey[p] * V + ez[p] * W;
		double u2 = U * U + V * V + W * W;
		feq[p] = w[p] * (Rho * (3.0 * x + 9.0 / 2.0 * x * x / epsF - 3.0 / 2.0 * u2 / epsF) + dRho);
	}
}

void Block::writeBlock(int id)
{
}

void Block::propagation()
{
}

void Block::collision()
{
}

void Block::makro()
{
}

int Block::getParInt(string quantity, unordered_map<string, string> par, bool mandatory, int defaultValue) {
	int val = 0;
	try {
		val = stoi(par[quantity]);
	}
	catch (const std::exception & e) {
		if (mandatory) {
			cout << "Error: variable " << quantity << " is not specified!" << endl;
			exit(EXIT_FAILURE);
		}
		else {
			cout << "variable " << quantity << " is set to default value " << defaultValue << endl;
			val = defaultValue;
		}
	}
	return val;
}

double Block::getParDouble(string quantity, unordered_map<string, string> par, bool mandatory, double defaultValue) {
	double val = 0;
	try {
		val = stod(par[quantity]);
	}
	catch (const std::exception & e) {
		if (mandatory) {
			cout << "Error: variable " << quantity << " is not specified!" << endl;
			exit(EXIT_FAILURE);
		}
		else {
			cout << "variable " << quantity << " is set to default value " << defaultValue << endl;
			val = defaultValue;
		}
	}
	return val;
}

double* Block::getParDoubleArray(string quantity, unordered_map<string, string> par) {
	double* out;
	try {
		string var = par[quantity];
		int s = 1;
		for (int i = 0; i < var.length(); i++) {
			if (var.at(i) == ',') {
				s++;
			}
		}

		out = new double[s];
		int startPos = 0;
		int endPos = 0;
		s = 0;
		while (endPos < var.length()) {
			if (var.at(endPos) == ',' || endPos == var.length()-1) {
				out[s] = stod(var.substr(startPos, endPos));
				s++;
				startPos = endPos + 1;
			}
			endPos++;
		}
	}
	catch (const std::exception & e) {
		cout << "Error: variable " << quantity << " is not specified!" << endl;
		exit(EXIT_FAILURE);
	}
	return out;
}