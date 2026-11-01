/*
 * File:   BlockFreeSurface.cpp
 * Author: obublik
 *
 * Created on 2. února 2018, 12:27
 */

#include "BlockFreeSurface.h"
#include <stdexcept>
#include <limits>

using namespace std;

#define FSB(b) ((BlockFreeSurface*)blocks[b])

BlockFreeSurface::BlockFreeSurface() {
	mass = vof = vofS1 = vofS = normal = curvature = nullptr;
	flagMap = surfaceLevels = filled = emptied = airToFs = fluidToFs = artefact1 = artefact2 = nullptr;
	fsType = nullptr;
	viscosity = shearRate = fluidModelParam = tauCell = nullptr;
	massErr = massBlock = 0;
	nFilled = nEmptied = nAirToFs = nFluidToFs = nArtefacts1 = nArtefacts2 = 0;
	numFsCells = 0;
	viscosityModel = 0;
	curvatureModel = 0;
	fsTolerance = 0.05;
	nAirToFsWarnings = 0;
	surfTension = wetAngleEq = wetAngleAdv = wetAngleRec = 0;
}

BlockFreeSurface::BlockFreeSurface(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par)
	: BlockFreeSurface() {
	initCommon(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, coarseLevel, indAll, indIn, indBound, collisionModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
	initFreeSurface(isFirst, viscosityModel, par);
}

void BlockFreeSurface::initFreeSurface(bool isFirst, int viscosityModel, const ParamMap& par) {
	this->viscosityModel = viscosityModel;
	massErr = 0.0;
	surfTension = getParDouble("surfaceTension", par, false, 0) / Cst;
	curvatureModel = getParInt("curvatureModel", par, false, 0);
	fsTolerance = getParDouble("fsTolerance", par, false, 0.05);

	wetAngleEq = getParDouble("wetAngleEq", par, false, 0) * PI / 180.0; // wetting angle - Equlibrium [rad]
	wetAngleAdv = getParDouble("wetAngleAdv", par, false, 0) * PI / 180.0; // wetting angle - Advancing  [rad]
	wetAngleRec = getParDouble("wetAngleRec", par, false, 0) * PI / 180.0; // wetting angle - Receding   [rad]

	mass = new double[ntot]();
	vof = new double[ntot]();
	vofS = new double[ntot]();
	vofS1 = new double[ntot]();
	normal = new double[ntot * 3]();
	curvature = new double[ntot]();
	flagMap = new int[ntot]();
	filled = new int[ntot];
	emptied = new int[ntot];
	airToFs = new int[ntot];
	fluidToFs = new int[ntot];
	artefact1 = new int[ntot];
	artefact2 = new int[ntot];
	surfaceLevels = new int[ntot]();
	fsType = new unsigned char[ntot]();
	tauCell = new double[ntot];
	for (int i = 0; i < ntot; i++) {
		tauCell[i] = tau;
	}

	if (viscosityModel > 0) {
		viscosity = new double[ntot]();
		shearRate = new double[ntot]();
		fluidModelParam = getParDoubleArray("fluidModelParameters", par);
		if (isFirst) {
			switch (viscosityModel) {
			case 1: // power-law
				cout << "Power-Law model: v_kinematic = b[0]*gamma^(b[1]-1)" << endl;
				break;
			case 2: // bingham
				cout << "Bingham model: v_kinematic = b[0]+b[1]*gamma^(b[2]-1)" << endl;
				break;
			case 3: // cross
				cout << "Cross model: v_kinematic = b[0]+(b[1]-b[0])/(1+b[2]*gamma^b[3])" << endl;
				break;
			case 4: // Poly
				cout << "Polynomial model: v_kinematic = b[0]+b[1]*gamma+b[2]*gamma^2+b[3]*gamma^3)" << endl;
				break;
			case 5: // Modified power-Law
				cout << "Modified power-Law: v_kinematic = b[0]*gamma^(b[1]-1 + b[2]*log(gamma))" << endl;
				break;
			case 6:
				cout << "Modified power-Law 2: v_kinematic = b[0]*gamma^(b[1]-1 + b[2]*log(gamma) + b[3]*log(gamma)^2)" << endl;
				break;
			default:
				throw invalid_argument("unknown viscosity model");
			}
			double gammas[] = { 0, 0.1, 1, 10, 1e4 };
			for (double g : gammas) {
				cout << "tau(gamma = " << g << ") = " << 1.0 / 2 + 3 * nonNewtonianViscosityModel(g) / Cvis << endl;
			}
			cout << "tau is limited to <" << tauMin << ", " << tauMax << ">" << endl;
		}
	}
}

BlockFreeSurface::~BlockFreeSurface() {
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
	delete[] fsType;
	delete[] tauCell;
	delete[] viscosity;
	delete[] shearRate;
	delete[] fluidModelParam;
}

void BlockFreeSurface::cellForce(int ind, double* acc, double& epsF) {
	acc[0] = G[0];
	acc[1] = G[1];
	acc[2] = G[2];
	epsF = 1.0;
}

double BlockFreeSurface::localTau(int ind) {
	if (viscosityModel > 0) {
		shearRate[ind] = computeShearRate(ind) / Ct;
		viscosity[ind] = nonNewtonianViscosityModel(shearRate[ind]);
		double tauLoc = 0.5 + 3 * viscosity[ind] / Cvis;
		return min(max(tauLoc, tauMin), tauMax);
	}
	return tau;
}

void BlockFreeSurface::collideCell(int ind, const double* acc, double epsF, double tauLoc) {
	double feq[27], feqNew[27], fneq[27], delta[27];
	if (epsF == 1.0) {
		equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind]);
		equilibrium(feqNew, Rho[ind], dRho[ind], U[ind] + acc[0], V[ind] + acc[1], W[ind] + acc[2]);
	}
	else {
		equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind], epsF);
		equilibrium(feqNew, Rho[ind], dRho[ind], U[ind] + acc[0], V[ind] + acc[1], W[ind] + acc[2], epsF);
	}
	const pop_t* fp = &f[nf * ind];
	for (int p = 0; p < nf; p++) {
		fneq[p] = fp[p] - feq[p];
	}
	coll.relax(fneq, tauLoc, delta);
	pop_t* fcp = &fc[nf * ind];
	for (int p = 0; p < nf; p++) {
		fcp[p] = (pop_t)(fp[p] - delta[p] + feqNew[p] - feq[p]);
	}
	tauCell[ind] = tauLoc;
}

void BlockFreeSurface::collision() {
	double acc[3];
	double epsF;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			cellForce(ind, acc, epsF);
			collideCell(ind, acc, epsF, localTau(ind));
		}
		else if (viscosityModel > 0) {
			shearRate[ind] = 0;
			viscosity[ind] = 0;
		}
	}
}

void BlockFreeSurface::outputVelocity(int ind, double* u) {
	if (flagMap[ind] == 1 || flagMap[ind] == 3) {
		double acc[3], epsF;
		cellForce(ind, acc, epsF);
		u[0] = U[ind] + 0.5 * acc[0];
		u[1] = V[ind] + 0.5 * acc[1];
		u[2] = W[ind] + 0.5 * acc[2];
	}
	else {
		u[0] = u[1] = u[2] = 0;
	}
}

// reconstruction of the population coming from the gas phase (Koerner et al. 2005)
void BlockFreeSurface::koernerPopulation(int ind, int p) {
	double u2 = U[ind] * U[ind] + V[ind] * V[ind] + W[ind] * W[ind];
	double x = ex[p] * U[ind] + ey[p] * V[ind] + ez[p] * W[ind];
	double xOpp = -x;
	double RhoS = 1.0 - 6.0 * surfTension * curvature[ind];
	double feq = w[p] * (RhoS * (3.0 * x + 9.0 / 2.0 * x * x - 3.0 / 2.0 * u2) + (RhoS - rho0));
	double feqOpp = w[opp[p]] * (RhoS * (3.0 * xOpp + 9.0 / 2.0 * xOpp * xOpp - 3.0 / 2.0 * u2) + (RhoS - rho0));
	f[nf * ind + opp[p]] = (pop_t)(feq + feqOpp - fc[nf * ind + p]);
}

void BlockFreeSurface::propagation() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] != 1 && flagMap[ind] != 3) {
			continue;
		}
		bool isInterface = flagMap[ind] == 3;
		for (int p = 0; p < nf; p++) {
			int bl, loc;
			int ct;
			if (nbIndex[ind] < 0) {
				bl = indBlock;
				loc = ind + offs[p];
				ct = cellType[loc];
			}
			else {
				ct = neighbourType(ind, p, bl, loc);
			}
			if (ct == CT_FLUID) {
				if (FSB(bl)->flagMap[loc] == 2) {
					koernerPopulation(ind, p);
				}
				else { // standart propagation
					f[nf * ind + opp[p]] = blocks[bl]->fc[nf * loc + opp[p]];
				}
			}
			else {
				double weight = normal[3 * ind] * ex[p] + normal[3 * ind + 1] * ey[p] + normal[3 * ind + 2] * ez[p];
				if (isInterface && weight > 0 && (ct == CT_WALL || ct == CT_MOVING)) {
					koernerPopulation(ind, p);
				}
				else {
					f[nf * ind + opp[p]] = (pop_t)boundaryIncoming(ind, p, ct);
				}
			}
		}
	}
}

void BlockFreeSurface::makro() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 1 || flagMap[ind] == 3) {
			double d = 0, jx = 0, jy = 0, jz = 0;
			const pop_t* fp = &f[nf * ind];
			for (int p = 0; p < nf; p++) {
				double fVal = fp[p];
				d += fVal;
				jx += ex[p] * fVal;
				jy += ey[p] * fVal;
				jz += ez[p] * fVal;
			}
			dRho[ind] = d;
			Rho[ind] = d + rho0;
			U[ind] = jx / Rho[ind];
			V[ind] = jy / Rho[ind];
			W[ind] = jz / Rho[ind];
		}
		else {
			Rho[ind] = 0;
			dRho[ind] = 0;
			U[ind] = 0;
			V[ind] = 0;
			W[ind] = 0;
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
	string fileName = simulPath + "init/blockInit" + std::to_string(id) + ".txt";
	ifstream myfile(fileName);
	if (!myfile.is_open()) {
		cout << "warning: initial condition " << fileName << " not found, block is empty" << endl;
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (!boundary[ind]) {
				flagMap[ind] = 2;
			}
		}
		return;
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
			if (tokens.size() < 5) {
				throw runtime_error("expected 5 values on a line");
			}
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
	catch (const std::exception& e) {
		throw runtime_error("Error in " + fileName + ": " + e.what());
	}
}

void BlockFreeSurface::writeOutCondition(int id, string simulPath) {
	string fileName = simulPath + "out/blockInit" + std::to_string(id) + ".txt";
	ofstream myfile(fileName);
	if (!myfile.is_open()) {
		cout << "Error: cannot write " << fileName << endl;
		return;
	}
	myfile.precision(12);
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		myfile << Rho[ind] << " " << U[ind] * Cu << " " << V[ind] * Cu << " " << W[ind] * Cu << " " << flagMap[ind] << "\n";
	}
}

void BlockFreeSurface::computeInitDistribution() {
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

void BlockFreeSurface::initFlagMap() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 1) {
			for (int dir = 1; dir < nf; dir++) {
				int bl = nbBlock(ind, dir);
				int loc = nbLocal(ind, dir);
				if (bl >= 0 && FSB(bl)->flagMap[loc] == 2) {
					flagMap[ind] = 3;
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
		vofS[ind] = vof[ind];
		vofS1[ind] = vof[ind];
	}
}

void BlockFreeSurface::computeInterfaceTypes() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] != 3) {
			continue;
		}
		bool hasFluid = false;
		bool hasGas = false;
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			if (bl >= 0 && !blocks[bl]->boundary[loc]) {
				int fl = FSB(bl)->flagMap[loc];
				hasFluid = hasFluid || fl == 1;
				hasGas = hasGas || fl == 2;
			}
		}
		if (!hasFluid && hasGas) {
			fsType[ind] = FS_NO_FLUID;
		}
		else if (!hasGas && hasFluid) {
			fsType[ind] = FS_NO_EMPTY;
		}
		else {
			fsType[ind] = FS_STANDARD;
		}
	}
}

void BlockFreeSurface::massExchange() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] != 3) {
			continue;
		}
		double dm = 0;
		for (int p = 1; p < nf; p++) {
			int bl, loc;
			int ct = neighbourType(ind, p, bl, loc);
			double fout = fc[nf * ind + p];
			if (ct != CT_FLUID) {
				if (ct != CT_WALL) {
					dm += boundaryIncoming(ind, p, ct) - fout; // inlet, outlet and moving walls
				}
				continue;
			}
			BlockFreeSurface* nb = FSB(bl);
			int nflag = nb->flagMap[loc];
			double fin = nb->fc[nf * loc + opp[p]];
			if (nflag == 1) {
				dm += fin - fout;
			}
			else if (nflag == 3) {
				// one-sided terms need the unshifted populations
				double finT = fin + w[p] * rho0;
				double foutT = fout + w[p] * rho0;
				double ex = 0;
				int my = fsType[ind];
				int other = nb->fsType[loc];
				if (my == other) {
					ex = fin - fout;
				}
				else if (my == FS_STANDARD) {
					ex = (other == FS_NO_FLUID) ? finT : -foutT;
				}
				else if (my == FS_NO_FLUID) {
					ex = -foutT;
				}
				else { // FS_NO_EMPTY
					ex = finT;
				}
				dm += ex * (vof[ind] + nb->vof[loc]) / 2;
			}
		}
		mass[ind] += dm;
	}
}

void BlockFreeSurface::findSurfaceLevels(int level) {
	if (level == 1) {
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			surfaceLevels[ind] = (flagMap[ind] == 3) ? 1 : 0;
		}
	}
	else {
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			if (surfaceLevels[ind] == level - 1) {
				for (int dir = 1; dir < nf; dir++) {
					int bl = nbBlock(ind, dir);
					int loc = nbLocal(ind, dir);
					if (bl >= 0 && !blocks[bl]->boundary[loc]) {
						int* sl = FSB(bl)->surfaceLevels;
						if (sl[loc] == 0) {
							// several threads may write the same value here
#pragma omp atomic write
							sl[loc] = level;
						}
					}
				}
			}
		}
	}
}

void BlockFreeSurface::findFilledEmptied(int iter) {
	int fsStatus = iter % 2;
	nFilled = 0;
	nEmptied = 0;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3) { // free surface cell
			double rho = rho0;
			for (int p = 0; p < nf; p++) {
				rho += fc[nf * ind + p];
			}
			vof[ind] = mass[ind] / rho;
			if (fsStatus == 0) { //filling
				if (vof[ind] > 1.0 + fsTolerance) {
					filled[nFilled++] = ind;
				}
			}
			else { //emptying
				if (vof[ind] < 0.0 - fsTolerance) {
					emptied[nEmptied++] = ind;
				}
			}
		}
	}
}

void BlockFreeSurface::vicinityFilledEmptied() {
	for (int cell = 0; cell < nFilled; cell++) {
		int ind = filled[cell];
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			if (bl >= 0 && !blocks[bl]->boundary[loc]) {
				BlockFreeSurface* nb = FSB(bl);
				if (nb->flagMap[loc] == 2) {
#pragma omp critical(CriticalSection)
					{
						if (nb->nAirToFs < nb->ntot) {
							nb->airToFs[nb->nAirToFs++] = loc; // cell is marked as airToFs
						}
					}
				}
			}
		}
	}

	for (int cell = 0; cell < nEmptied; cell++) {
		int ind = emptied[cell];
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			if (bl >= 0 && !blocks[bl]->boundary[loc]) {
				BlockFreeSurface* nb = FSB(bl);
				if (nb->flagMap[loc] == 1) {
#pragma omp critical(CriticalSection)
					{
						if (nb->nFluidToFs < nb->ntot) {
							nb->fluidToFs[nb->nFluidToFs++] = loc; // cell is marked as fluidToFs
						}
					}
				}
			}
		}
	}
}

void BlockFreeSurface::calculateNormalVectors() {
	double normTol = 1e-5;
	double auxVof[27];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (surfaceLevels[ind] > 0 && surfaceLevels[ind] < 3) {
			for (int dir = 1; dir < nf; dir++) {
				int bl = nbBlock(ind, dir);
				int loc = nbLocal(ind, dir);
				if (bl >= 0 && !blocks[bl]->boundary[loc]) {
					auxVof[dir] = FSB(bl)->vofS[loc];
					continue;
				}
				bl = nbBlock(ind, opp[dir]);
				loc = nbLocal(ind, opp[dir]);
				if (bl >= 0 && !blocks[bl]->boundary[loc]) {
					auxVof[dir] = 2 * vofS[ind] - FSB(bl)->vofS[loc];
				}
				else {
					auxVof[dir] = vofS[ind];
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
}

void BlockFreeSurface::setWetAngle() {
	if (wetAngleEq == 0.0 || surfTension == 0.0) {
		return;
	}
	double wetAngleDyn;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (boundary[ind] || surfaceLevels[ind] <= 0) {
			continue;
		}
		int tmpWallNormal[3] = { 0, 0, 0 };
		int numWallCells = 0;
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			if (bl < 0 || blocks[bl]->boundary[loc]) {
				tmpWallNormal[0] += ex[dir];
				tmpWallNormal[1] += ey[dir];
				tmpWallNormal[2] += ez[dir];
				numWallCells++;
			}
		}

		double wallNormal[3];
		bool isContactCell = false;
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
					double weighTot = 0.0;
					for (int dir = 1; dir < nf; dir++) {
						int bl = nbBlock(ind, dir);
						int loc = nbLocal(ind, dir);
						if (bl >= 0) {
							int fl = FSB(bl)->flagMap[loc];
							if (fl == 1 || fl == 3) {
								double weight = max(FSB(bl)->vof[loc], 0.0);
								normalVelMag += (blocks[bl]->U[loc] * wallTang[0] + blocks[bl]->V[loc] * wallTang[1] + blocks[bl]->W[loc] * wallTang[2]) * weight;
								weighTot += weight;
							}
						}
					}
					normalVelMag = (weighTot != 0.0) ? normalVelMag / weighTot : 0.0;
				}
				else {
					normalVelMag = U[ind] * wallTang[0] + V[ind] * wallTang[1] + W[ind] * wallTang[2];
				}

				double Ca = vLB * normalVelMag / surfTension;
				double ka = 9e-9;
				double kr = 9e-8;
				if (normalVelMag >= 0.0) {
					wetAngleDyn = min(wetAngleEq + cbrt(Ca / ka) * (PI / 180.0), wetAngleAdv);
				}
				else {
					wetAngleDyn = max(wetAngleEq + cbrt(Ca / kr) * (PI / 180.0), wetAngleRec);
				}
				normal[3 * ind] = cos(wetAngleDyn) * wallNormal[0] + sin(wetAngleDyn) * wallTang[0];
				normal[3 * ind + 1] = cos(wetAngleDyn) * wallNormal[1] + sin(wetAngleDyn) * wallTang[1];
				normal[3 * ind + 2] = cos(wetAngleDyn) * wallNormal[2] + sin(wetAngleDyn) * wallTang[2];
			}
		}
	}
}

void BlockFreeSurface::smoothVOF(int smoothStep, int nSmoothStep) {
	switch (nSmoothStep) {
	case 0:
		if (smoothStep == 0) {
			for (int i = 0; i < nAll; i++) {
				vofS[indAll[i]] = vof[indAll[i]];
			}
		}
		break;
	case 1:
		if (smoothStep == 1) {
			smooth3D(vof, vofS, 1);
		}
		break;
	case 2:
		if (smoothStep == 1) {
			smooth3D(vof, vofS1, 1);
		}
		else if (smoothStep == 2) {
			smooth3D(vofS1, vofS, 2);
		}
		break;
	}
}

void BlockFreeSurface::smooth3D(double* arrayin, double* arrayout, int kod) {
	double auxVof[27];
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (surfaceLevels[ind] > 0) {
			for (int dir = 0; dir < nf; dir++) {
				int bl = nbBlock(ind, dir);
				int loc = nbLocal(ind, dir);
				if (bl >= 0 && !blocks[bl]->boundary[loc]) {
					auxVof[dir] = (kod == 1) ? FSB(bl)->vof[loc] : FSB(bl)->vofS1[loc];
					continue;
				}
				bl = nbBlock(ind, opp[dir]);
				loc = nbLocal(ind, opp[dir]);
				if (bl >= 0 && !blocks[bl]->boundary[loc]) {
					auxVof[dir] = 2 * arrayin[ind] - ((kod == 1) ? FSB(bl)->vof[loc] : FSB(bl)->vofS1[loc]);
				}
				else {
					auxVof[dir] = arrayin[ind];
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

// volume fraction of the cell shifted by (dx, dy, dz), NaN for walls and cells outside of the domain
double BlockFreeSurface::vofAt(int ind, int dx, int dy, int dz) {
	int bl, loc;
	local2globalShift(indBlock, ind, bl, loc, dx, dy, dz);
	if (bl < 0 || bl >= nBlockMax || !blocks[bl] || blocks[bl]->boundary[loc]) {
		return numeric_limits<double>::quiet_NaN();
	}
	BlockFreeSurface* b = FSB(bl);
	switch (b->flagMap[loc]) {
	case 1:
		return 1.0;
	case 2:
		return 0.0;
	case 3:
		return min(max(b->vof[loc], 0.0), 1.0);
	default:
		return numeric_limits<double>::quiet_NaN();
	}
}

// curvature from the divergence of the interface normals
double BlockFreeSurface::curvatureDivergence(int ind) {
	double auxNx[27], auxNy[27], auxNz[27];
	for (int dir = 1; dir < nf; dir++) {
		int bl = nbBlock(ind, dir);
		int loc = nbLocal(ind, dir);
		if (bl >= 0 && !blocks[bl]->boundary[loc]) {
			auxNx[dir] = FSB(bl)->normal[3 * loc];
			auxNy[dir] = FSB(bl)->normal[3 * loc + 1];
			auxNz[dir] = FSB(bl)->normal[3 * loc + 2];
			continue;
		}
		bl = nbBlock(ind, opp[dir]);
		loc = nbLocal(ind, opp[dir]);
		if (bl >= 0 && !blocks[bl]->boundary[loc]) {
			auxNx[dir] = 2 * normal[3 * ind] - FSB(bl)->normal[3 * loc];
			auxNy[dir] = 2 * normal[3 * ind + 1] - FSB(bl)->normal[3 * loc + 1];
			auxNz[dir] = 2 * normal[3 * ind + 2] - FSB(bl)->normal[3 * loc + 2];
			double norm = sqrt(auxNx[dir] * auxNx[dir] + auxNy[dir] * auxNy[dir] + auxNz[dir] * auxNz[dir]);
			if (norm > 1e-5) {
				auxNx[dir] /= norm;
				auxNy[dir] /= norm;
				auxNz[dir] /= norm;
			}
		}
		else {
			auxNx[dir] = normal[3 * ind];
			auxNy[dir] = normal[3 * ind + 1];
			auxNz[dir] = normal[3 * ind + 2];
		}
	}

	double curv = 0.0;
	for (int dir = 1; dir < nf; dir++) {
		curv += w[dir] * (ex[dir] * auxNx[dir] + ey[dir] * auxNy[dir] + ez[dir] * auxNz[dir]);
	}
	curv *= -3.0;
	curv = min(max(curv, -4.0), 4.0);
	// stored curvature is minus the mean curvature: Laplace pressure 2*sigma*H = -6*sigma*curvature/3
	return curv / 2.0;
}

// height function curvature (Cummins et al. 2005), returns false if the stencil is not usable
bool BlockFreeSurface::curvatureHeightFunction(int ind, double& curv) {
	double n[3] = { normal[3 * ind], normal[3 * ind + 1], normal[3 * ind + 2] };
	int a = 0;
	for (int d = 1; d < 3; d++) {
		if (fabs(n[d]) > fabs(n[a])) {
			a = d;
		}
	}
	if (fabs(n[a]) < 1e-6) {
		return false;
	}
	if (nf == 9 && a == 2) {
		return false;
	}
	int b = (a + 1) % 3;
	int c = (a + 2) % 3;
	bool twoD = (nf == 9);
	if (twoD) { // in 2D the lateral direction is the remaining in-plane axis
		b = (a == 0) ? 1 : 0;
	}
	int s = n[a] > 0 ? 1 : -1; // direction towards the gas
	const int L = 3;
	double H[3][3];
	for (int db = -1; db <= 1; db++) {
		for (int dc = -1; dc <= 1; dc++) {
			if (twoD && dc != 0) {
				H[db + 1][dc + 1] = 0;
				continue;
			}
			double h = 0;
			for (int t = -L; t <= L; t++) {
				int d[3] = { 0, 0, 0 };
				d[a] = t;
				d[b] += db;
				d[c] += dc;
				double v = vofAt(ind, d[0], d[1], d[2]);
				if (std::isnan(v)) {
					return false;
				}
				if (db == 0 && dc == 0) {
					// the central column has to contain the whole interface
					if (t == -L * s && v < 0.5) {
						return false;
					}
					if (t == L * s && v > 0.5) {
						return false;
					}
				}
				h += v;
			}
			H[db + 1][dc + 1] = h;
		}
	}
	double hb = (H[2][1] - H[0][1]) / 2;
	double hbb = H[2][1] - 2 * H[1][1] + H[0][1];
	double hc = 0, hcc = 0, hbc = 0;
	if (!twoD) {
		hc = (H[1][2] - H[1][0]) / 2;
		hcc = H[1][2] - 2 * H[1][1] + H[1][0];
		hbc = (H[2][2] - H[2][0] - H[0][2] + H[0][0]) / 4;
	}
	double denom = pow(1 + hb * hb + hc * hc, 1.5);
	double kappa = (hbb * (1 + hc * hc) + hcc * (1 + hb * hb) - 2 * hbc * hb * hc) / denom;
	curv = min(max(kappa / 2.0, -2.0), 2.0);
	return true;
}

void BlockFreeSurface::calculateCurvature() {
	if (surfTension == 0.0) {
		return;
	}
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3) {
			double curv;
			if (curvatureModel == 1 && curvatureHeightFunction(ind, curv)) {
				curvature[ind] = curv;
			}
			else {
				curvature[ind] = curvatureDivergence(ind);
			}
		}
	}
}

void BlockFreeSurface::newFsGeometryPrepare() {
	newCellState.assign(4 * nAirToFs, 0.0);
	for (int cell = 0; cell < nAirToFs; cell++) {
		int ind = airToFs[cell];
		double count = 0.0;
		double* st = &newCellState[4 * cell];
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			if (bl < 0) {
				continue;
			}
			BlockFreeSurface* nb = FSB(bl);
			if (nb->flagMap[loc] == 1 || nb->flagMap[loc] == 3) {
				double weight = w[dir] * min(max(nb->vof[loc], 0.0), 1.0);
				count += weight;
				st[0] += nb->Rho[loc] * weight;
				st[1] += nb->U[loc] * weight;
				st[2] += nb->V[loc] * weight;
				st[3] += nb->W[loc] * weight;
			}
		}
		if (count <= 0) {
			nAirToFsWarnings++;
			st[0] = 1.0;
			st[1] = st[2] = st[3] = 0;
		}
		else {
			for (int k = 0; k < 4; k++) {
				st[k] /= count;
			}
		}
	}
}

void BlockFreeSurface::newFsGeometryApply() {
	double feq[27];
	for (int cell = 0; cell < nAirToFs; cell++) {
		int ind = airToFs[cell];
		const double* st = &newCellState[4 * cell];
		flagMap[ind] = 3;
		vof[ind] = 0.0;
		mass[ind] = 0.0;
		Rho[ind] = st[0];
		dRho[ind] = Rho[ind] - rho0;
		U[ind] = st[1];
		V[ind] = st[2];
		W[ind] = st[3];
		equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind]);
		for (int p = 0; p < nf; p++) { // initializace equlibriem
			fc[nf * ind + p] = (pop_t)feq[p];
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
		curvature[ind] = 0.0;
	}

	for (int cell = 0; cell < nFluidToFs; cell++) {
		int ind = fluidToFs[cell];
		if (flagMap[ind] == 1) {
			flagMap[ind] = 3;
			mass[ind] = Rho[ind];
			vof[ind] = 1.0;
		}
	}
}

void BlockFreeSurface::newFsGeometry() {
	newFsGeometryPrepare();
	newFsGeometryApply();
}

void BlockFreeSurface::findArtefactsEmptied() {
	nArtefacts1 = 0;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3) {
			bool isArtefact = true;
			int nAir = 0;
			for (int dir = 1; dir < nf; dir++) {
				int bl = nbBlock(ind, dir);
				int loc = nbLocal(ind, dir);
				if (bl >= 0) {
					int fl = FSB(bl)->flagMap[loc];
					if (fl == 1 || blocks[bl]->cellType[loc] == CT_INLET) {
						isArtefact = false;
					}
					else if (fl == 2) {
						nAir++;
					}
				}
			}
			if (isArtefact) {
				// neither a filled artefact nor next to one (it would leave fluid touching gas)
				for (int dir = 0; dir < nf && isArtefact; dir++) {
					int bl = nbBlock(ind, dir);
					int loc = nbLocal(ind, dir);
					if (bl >= 0 && FSB(bl)->fsType[loc] == FS_ARTEFACT_FILLED) {
						isArtefact = false;
					}
				}
			}
			if (isArtefact && nAir > 0) {
				artefact1[nArtefacts1++] = ind; // cell is marked as artefact
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
				int bl = nbBlock(ind, dir);
				int loc = nbLocal(ind, dir);
				if (bl >= 0) {
					if (FSB(bl)->flagMap[loc] == 2) {
						isArtefact = false;
					}
					else if (blocks[bl]->boundary[loc] && weight > 0.6) {
						isArtefact = false;
					}
				}
				else if (weight > 0.6) {
					isArtefact = false;
				}
			}
			if (isArtefact) {
				artefact2[nArtefacts2++] = ind; // cell is marked as artefact
				fsType[ind] = FS_ARTEFACT_FILLED;
			}
		}
	}
}

void BlockFreeSurface::destroyArtefacts() {
	for (int cell = 0; cell < nArtefacts2; cell++) {
		int ind = artefact2[cell];
		fsType[ind] = FS_STANDARD;
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

// safety net: fluid cells must not touch gas, such cells become interface cells
void BlockFreeSurface::fixFluidGasContacts() {
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] != 1) {
			continue;
		}
		for (int dir = 1; dir < nf; dir++) {
			int bl = nbBlock(ind, dir);
			int loc = nbLocal(ind, dir);
			if (bl >= 0 && FSB(bl)->flagMap[loc] == 2) {
				flagMap[ind] = 3;
				vof[ind] = 1.0;
				mass[ind] = Rho[ind];
				fsType[ind] = FS_STANDARD;
				break;
			}
		}
	}
}

void BlockFreeSurface::massRedistribution() {
	for (int cell = 0; cell < nEmptied; cell++) {
		int ind = emptied[cell];
		redistribute(ind, mass[ind], -1);
		mass[ind] = 0;
	}
	for (int cell = 0; cell < nFilled; cell++) {
		int ind = filled[cell];
		redistribute(ind, mass[ind] - Rho[ind], 1);
		mass[ind] = Rho[ind];
	}
	for (int cell = 0; cell < nArtefacts1; cell++) {
		int ind = artefact1[cell];
		redistribute(ind, mass[ind], 0);
		mass[ind] = 0;
	}
	for (int cell = 0; cell < nArtefacts2; cell++) {
		int ind = artefact2[cell];
		redistribute(ind, mass[ind] - Rho[ind], 0);
		mass[ind] = Rho[ind];
	}
}

// distributes the mass of a converted cell to the neighbouring interface cells,
// direction 1 - along the normal (towards the gas), -1 - against the normal, 0 - uniformly
void BlockFreeSurface::redistribute(int ind, double massRedist, int direction) {
	int indNeigh[27];
	int indBlockNeigh[27];
	double listWi[27];
	int countNeigh = 0;
	double wtot = 0.0;
	double n[3] = { normal[3 * ind], normal[3 * ind + 1], normal[3 * ind + 2] };

	for (int dx = -1; dx < 2; dx++) {
		for (int dy = -1; dy < 2; dy++) {
			for (int dz = -1; dz < 2; dz++) {
				if (dx == 0 && dy == 0 && dz == 0) {
					continue;
				}
				int bl, loc;
				local2globalShift(indBlock, ind, bl, loc, dx, dy, dz);
				if (bl > -1 && bl < nBlockMax && blocks[bl] && FSB(bl)->flagMap[loc] == 3) {
					double wi = direction * (n[0] * dx + n[1] * dy + n[2] * dz);
					indNeigh[countNeigh] = loc;
					indBlockNeigh[countNeigh] = bl;
					listWi[countNeigh] = max(wi, 0.0);
					countNeigh++;
					wtot += max(wi, 0.0);
				}
			}
		}
	}
	if (countNeigh == 0) {
		massErr += massRedist;
		return;
	}
	if (wtot <= 1e-12) { // uniform distribution
		for (int ng = 0; ng < countNeigh; ng++) {
			listWi[ng] = 1.0;
		}
		wtot = countNeigh;
	}
	for (int ng = 0; ng < countNeigh; ng++) {
		double dm = (listWi[ng] / wtot) * massRedist;
		double* m = &FSB(indBlockNeigh[ng])->mass[indNeigh[ng]];
#pragma omp atomic
		*m += dm;
	}
}

// mass that could not be redistributed is spread over the interface cells of the block
void BlockFreeSurface::redistributeMassErr() {
	if (massErr == 0.0) {
		return;
	}
	int count = 0;
	for (int i = 0; i < nAll; i++) {
		if (flagMap[indAll[i]] == 3) {
			count++;
		}
	}
	if (count == 0) {
		return;
	}
	double dm = massErr / count;
	for (int i = 0; i < nAll; i++) {
		int ind = indAll[i];
		if (flagMap[ind] == 3) {
			mass[ind] += dm;
		}
	}
	massErr = 0.0;
}

void BlockFreeSurface::calculMassTotal() {
	massBlock = massErr;
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
				int bl = nbBlock(ind, dir);
				if (bl >= 0 && bl != indBlock) {
#pragma omp atomic write
					blocks[bl]->isBlockTouched = true;
				}
			}
		}
	}
}

// shear rate from the non-equilibrium stress, S = -3/(2 rho tau) sum e e fneq
double BlockFreeSurface::computeShearRate(int ind) {
	if (flagMap[ind] != 1 && flagMap[ind] != 3) {
		return 0;
	}
	double feq[27];
	equilibrium(feq, Rho[ind], dRho[ind], U[ind], V[ind], W[ind]);
	double pxx = 0, pyy = 0, pzz = 0, pxy = 0, pxz = 0, pyz = 0;
	for (int p = 0; p < nf; p++) {
		double fneq = f[nf * ind + p] - feq[p];
		pxx += ex[p] * ex[p] * fneq;
		pyy += ey[p] * ey[p] * fneq;
		pzz += ez[p] * ez[p] * fneq;
		pxy += ex[p] * ey[p] * fneq;
		pxz += ex[p] * ez[p] * fneq;
		pyz += ey[p] * ez[p] * fneq;
	}
	double c = -3.0 / (2.0 * Rho[ind] * tauCell[ind]);
	double sxx = c * pxx, syy = c * pyy, szz = c * pzz, sxy = c * pxy, sxz = c * pxz, syz = c * pyz;
	return sqrt(2.0 * (sxx * sxx + syy * syy + szz * szz + 2.0 * (sxy * sxy + sxz * sxz + syz * syz)));
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
