/*
 * File:   Block.cpp
 * Author: obublik
 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <climits>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "Block.h"

using namespace std;

Block::~Block() {
	delete[] boundary;
	delete[] cellType;
	delete[] fc;
	delete[] f;
	delete[] Rho;
	delete[] dRho;
	delete[] U;
	delete[] V;
	delete[] W;
	delete[] wallDistance;
	delete[] wallDistanceGradient;
	delete[] wallDistancePolynome;
	delete[] G;
	delete[] G0;
	delete[] nbIndex;
	delete[] nbBlockTab;
	delete[] nbLocalTab;
}

bool Block::isInterior(int i, int j, int k, int ni, int nj, int nk) {
	// a dimension of size one (2D model) does not restrict the interior
	bool in = (ni == 1 || (i > 0 && i < ni - 1));
	in = in && (nj == 1 || (j > 0 && j < nj - 1));
	in = in && (nk == 1 || (k > 0 && k < nk - 1));
	return in;
}

void Block::countCells(int ni, int nj, int nk, int& nIn, int& nBound) {
	nIn = 0;
	for (int k = 0; k < nk; k++) {
		for (int j = 0; j < nj; j++) {
			for (int i = 0; i < ni; i++) {
				if (isInterior(i, j, k, ni, nj, nk)) {
					nIn++;
				}
			}
		}
	}
	nBound = ni * nj * nk - nIn;
}

void Block::initCommon(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int* indAll, int* indIn, int* indBound, int collisionModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par) {
	if (coarseLevel != 1) {
		throw invalid_argument("grid refinement (coarseLevel != 1) is not implemented");
	}
	this->blocks = blocks;
	this->nBlock = nBlock;
	this->nBlockMax = nBlockMax;
	this->blockMap = blockMap;
	this->nf = nf;
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
	this->ni = ni;
	this->nj = nj;
	this->nk = nk;
	nAll = ni * nj * nk;
	ntot = nAll;
	countCells(ni, nj, nk, nIn, nBound);
	isBlockEmpty = false;

	xPeriodicity = getParInt("xPeriodicity", par, false, 0) != 0;
	yPeriodicity = getParInt("yPeriodicity", par, false, 0) != 0;
	zPeriodicity = getParInt("zPeriodicity", par, false, 0) != 0;
	if (nk == 1 && nBk == 1) {
		zPeriodicity = true; // 2D model: the single layer is its own neighbour in z
	}

	Cu = getParDouble("uMax", par, true, 0) / 0.1;
	Ch = getParDouble("H", par, true, 0) / ni;
	Ct = Ch / Cu;
	Crho = getParDouble("density", par, true, 0);
	Cp = Crho / 3 * Cu * Cu;
	Cst = Crho * Cu * Cu * Ch; // Cp*Ch
	Cvis = Ch * Ch / Ct;
	vLB = getParDouble("vKin", par, true, 0) / Cvis;
	Cf = Ch / (Ct * Ct);
	tau = 1.0 / 2 + 3 * vLB;

	G = new double[3]{
		getParDouble("gx", par, false, 0) / Cf, getParDouble("gy", par, false, 0) / Cf, getParDouble("gz", par, false, 0) / Cf
	};
	G0 = new double[3]{ G[0], G[1], G[2] };
	velWall[0] = getParDouble("uwall", par, false, 0) / Cu;
	velWall[1] = getParDouble("vwall", par, false, 0) / Cu;
	velWall[2] = getParDouble("wwall", par, false, 0) / Cu;
	velInlet[0] = getParDouble("uInlet", par, false, 0) / Cu;
	velInlet[1] = getParDouble("vInlet", par, false, 0) / Cu;
	velInlet[2] = getParDouble("wInlet", par, false, 0) / Cu;
	rhoOutlet = 1.0 + getParDouble("pOutlet", par, false, 0) / Cp;
	movingLidTop = getParInt("movingLidTop", par, false, 1) != 0;

	rho0 = getParInt("ddfShift", par, false, 1) ? 1.0 : 0.0;

	if (hasPar("collisionModel", par)) {
		this->collisionModel = getParInt("collisionModel", par, false, collisionModel);
	}
	else {
		this->collisionModel = collisionModel;
	}
	mrtModel = this->collisionModel == COLLISION_MRT;
	tauMin = getParDouble("tauMin", par, false, 0.5001);
	tauMax = getParDouble("tauMax", par, false, 1.5);
	smagorinsky = getParDouble("smagorinsky", par, false, 0);
	coll.init(this->collisionModel, nf, ex, ey, ez, opp, getParDouble("trtMagic", par, false, 3.0 / 16.0), getParInt("mrtMagicSq", par, false, 1) != 0, tau, getParDouble("tauMinusMax", par, false, 1e30));

	for (int dir = 0; dir < nf; dir++) {
		offs[dir] = (ni * nj) * ez[dir] + ni * ey[dir] + ex[dir];
	}

	boundary = new bool[ntot]();
	cellType = new unsigned char[ntot]();
	fc = new pop_t[ntot * nf]();
	f = new pop_t[ntot * nf]();
	Rho = new double[ntot]();
	dRho = new double[ntot]();
	U = new double[ntot]();
	V = new double[ntot]();
	W = new double[ntot]();
	wallDistance = new double[ntot]();
	wallDistanceGradient = new double[3 * ntot]();
	wallDistancePolynome = new double[8 * ntot]();

	if (isFirst) {
		const char* names[] = { "BGK", "MRT", "TRT" };
		cout << "collision model = " << names[this->collisionModel] << endl;
		cout << "tau = " << tau << endl;
		if (smagorinsky > 0) {
			cout << "Smagorinsky LES, C = " << smagorinsky << endl;
		}
		if (tau < 0.51 && smagorinsky == 0) {
			cout << "warning: tau is close to 0.5, the simulation may be unstable (consider smagorinsky = 0.1)" << endl;
		}
	}
}

void Block::local2globalShift(int indBlock, int ind, int& indBlockNeigh, int& indLocalNeigh, int ex, int ey, int ez) {
	int k = ind / (ni * nj) + ez;
	int pom = ind % (ni * nj);
	int j = pom / ni + ey;
	int i = pom % ni + ex;

	int kB = indBlock / (nBi * nBj);
	pom = indBlock % (nBi * nBj);
	int jB = pom / nBi;
	int iB = pom % nBi;

	while (i < 0) {
		i += ni;
		iB--;
		if (xPeriodicity && iB == -1) {
			iB = nBi - 1;
		}
	}
	while (i >= ni) {
		i -= ni;
		iB++;
		if (xPeriodicity && iB == nBi) {
			iB = 0;
		}
	}
	while (j < 0) {
		j += nj;
		jB--;
		if (yPeriodicity && jB == -1) {
			jB = nBj - 1;
		}
	}
	while (j >= nj) {
		j -= nj;
		jB++;
		if (yPeriodicity && jB == nBj) {
			jB = 0;
		}
	}
	while (k < 0) {
		k += nk;
		kB--;
		if (zPeriodicity && kB == -1) {
			kB = nBk - 1;
		}
	}
	while (k >= nk) {
		k -= nk;
		kB++;
		if (zPeriodicity && kB == nBk) {
			kB = 0;
		}
	}

	indBlockNeigh = nBi * nBj * kB + nBi * jB + iB;
	if (iB < 0 || iB >= nBi || jB < 0 || jB >= nBj || kB < 0 || kB >= nBk) {
		indBlockNeigh = -1;
	}
	indLocalNeigh = ni * nj * k + ni * j + i;
}

void Block::computeShiftIndexes() {
	delete[] nbIndex;
	delete[] nbBlockTab;
	delete[] nbLocalTab;
	nbIndex = new int[ntot];
	nbBlockTab = new int[nBound * nf];
	nbLocalTab = new int[nBound * nf];
	for (int i = 0; i < ntot; i++) {
		nbIndex[i] = -1;
	}
	for (int i = 0; i < nBound; i++) {
		int ind = indBound[i];
		nbIndex[ind] = i;
		for (int dir = 0; dir < nf; dir++) {
			int bl, loc;
			local2globalShift(indBlock, ind, bl, loc, ex[dir], ey[dir], ez[dir]);
			if (bl < 0 || bl >= nBlockMax || !blocks[bl]) {
				bl = -1;
			}
			nbBlockTab[nf * i + dir] = bl;
			nbLocalTab[nf * i + dir] = loc;
		}
	}
}

int Block::neighbourType(int ind, int dir, int& bl, int& loc) const {
	bl = nbBlock(ind, dir);
	loc = nbLocal(ind, dir);
	if (bl < 0 || bl >= nBlockMax || !blocks[bl]) {
		bl = -1;
		if (movingLidTop && ez[dir] > 0) {
			int k = ind / (ni * nj);
			int kB = indBlock / (nBi * nBj);
			if (k == nk - 1 && kB == nBk - 1) {
				return CT_MOVING;
			}
		}
		return CT_WALL;
	}
	return blocks[bl]->cellType[loc];
}

double Block::boundaryIncoming(int ind, int dir, int ctype) const {
	double fout = fc[nf * ind + dir];
	switch (ctype) {
	case CT_INLET:
		return fout - 6.0 * w[dir] * Rho[ind] * (ex[dir] * velInlet[0] + ey[dir] * velInlet[1] + ez[dir] * velInlet[2]);
	case CT_MOVING:
		return fout - 6.0 * w[dir] * Rho[ind] * (ex[dir] * velWall[0] + ey[dir] * velWall[1] + ez[dir] * velWall[2]);
	case CT_OUTLET: {
		// anti bounce-back with prescribed density
		double x = ex[dir] * U[ind] + ey[dir] * V[ind] + ez[dir] * W[ind];
		double u2 = U[ind] * U[ind] + V[ind] * V[ind] + W[ind] * W[ind];
		return -fout + 2.0 * w[dir] * (rhoOutlet * (1.0 + 4.5 * x * x - 1.5 * u2) - rho0);
	}
	default:
		return fout; // bounce-back
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
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		boundary[ind] = false;
		cellType[ind] = CT_FLUID;
	}
	string fileName = geomPath + "mesh/boundary/blockBoundary" + std::to_string(id) + ".txt";
	ifstream myfile(fileName);
	if (!myfile.is_open()) {
		return; // no boundary file: the whole block is fluid
	}
	string line;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (!getline(myfile, line)) {
			throw runtime_error("unexpected end of file " + fileName);
		}
		int type = 0;
		istringstream iss(line);
		if (!(iss >> type) || type < 0 || type > CT_MOVING) {
			throw runtime_error("invalid cell type in " + fileName + ": '" + line + "'");
		}
		cellType[ind] = (unsigned char)type;
		boundary[ind] = type != CT_FLUID;
	}
}

double Block::lesTau(double tau0, const double* fneq, double rho) const {
	double pxx = 0, pyy = 0, pzz = 0, pxy = 0, pxz = 0, pyz = 0;
	for (int p = 0; p < nf; p++) {
		pxx += ex[p] * ex[p] * fneq[p];
		pyy += ey[p] * ey[p] * fneq[p];
		pzz += ez[p] * ez[p] * fneq[p];
		pxy += ex[p] * ey[p] * fneq[p];
		pxz += ex[p] * ez[p] * fneq[p];
		pyz += ey[p] * ez[p] * fneq[p];
	}
	double Q = sqrt(pxx * pxx + pyy * pyy + pzz * pzz + 2 * (pxy * pxy + pxz * pxz + pyz * pyz));
	double c2 = smagorinsky * smagorinsky;
	return 0.5 * (tau0 + sqrt(tau0 * tau0 + 18.0 * sqrt(2.0) * c2 * Q / rho));
}

void Block::equilibrium(double* feq, double Rho, double dRho, double U, double V, double W) {
	double u2 = U * U + V * V + W * W;
	for (int p = 0; p < nf; p++) {
		double x = ex[p] * U + ey[p] * V + ez[p] * W;
		feq[p] = w[p] * (Rho * (3.0 * x + 9.0 / 2.0 * x * x - 3.0 / 2.0 * u2) + dRho);
	}
}

void Block::equilibrium(double* feq, double Rho, double dRho, double U, double V, double W, double epsF) {
	double u2 = U * U + V * V + W * W;
	for (int p = 0; p < nf; p++) {
		double x = ex[p] * U + ey[p] * V + ez[p] * W;
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

void Block::outputVelocity(int ind, double* u) {
	u[0] = U[ind];
	u[1] = V[ind];
	u[2] = W[ind];
}

bool Block::hasPar(const string& quantity, const ParamMap& par) {
	auto it = par.find(quantity);
	return it != par.end() && !it->second.empty();
}

static double parseDouble(const string& quantity, const ParamMap& par, bool mandatory, double defaultValue, bool& found) {
	found = false;
	auto it = par.find(quantity);
	if (it != par.end() && !it->second.empty()) {
		try {
			size_t pos = 0;
			double val = stod(it->second, &pos);
			found = true;
			return val;
		}
		catch (const std::exception&) {
			cout << "Error: variable " << quantity << " has invalid value '" << it->second << "'" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (mandatory) {
		cout << "Error: variable " << quantity << " is not specified!" << endl;
		exit(EXIT_FAILURE);
	}
	return defaultValue;
}

int Block::getParInt(const string& quantity, const ParamMap& par, bool mandatory, int defaultValue) {
	bool found;
	double val = parseDouble(quantity, par, mandatory, defaultValue, found);
	if (val > INT_MAX) {
		return INT_MAX;
	}
	if (val < INT_MIN) {
		return INT_MIN;
	}
	return (int)llround(val);
}

double Block::getParDouble(const string& quantity, const ParamMap& par, bool mandatory, double defaultValue) {
	bool found;
	return parseDouble(quantity, par, mandatory, defaultValue, found);
}

double* Block::getParDoubleArray(const string& quantity, const ParamMap& par) {
	auto it = par.find(quantity);
	if (it == par.end() || it->second.empty()) {
		cout << "Error: variable " << quantity << " is not specified!" << endl;
		exit(EXIT_FAILURE);
	}
	string var = it->second;
	for (char& c : var) {
		if (c == ',' || c == ';') {
			c = ' ';
		}
	}
	vector<double> values;
	istringstream iss(var);
	string token;
	while (iss >> token) {
		try {
			values.push_back(stod(token));
		}
		catch (const std::exception&) {
			cout << "Error: variable " << quantity << " has invalid value '" << it->second << "'" << endl;
			exit(EXIT_FAILURE);
		}
	}
	double* out = new double[values.size() + 1];
	for (size_t i = 0; i < values.size(); i++) {
		out[i] = values[i];
	}
	out[values.size()] = 0;
	return out;
}
