/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Mesh.h"
#include "MarchingCubes.h"
#include "Block.h"

using namespace std;

Mesh::Mesh(string geomPath, string simulPath, unordered_map<string, string> par) {
	this->geomPath = geomPath;
	this->simulPath = simulPath;
	Block b;
	ni = b.getParInt("ni", par, true, 0);
	nj = b.getParInt("nj", par, true, 0);
	nk = b.getParInt("nk", par, true, 0);
	problemID = b.getParInt("problemID", par, true, 0);
	viscosityModel = b.getParInt("viscosityModel", par, false, 0);
	initHydroPress = (b.getParInt("initializeWithHydroPress", par, false, 0) ? true : false);
	nSmoothStep = 2; // maximum je 2
	nIntIterParticle = b.getParInt("nIntIter", par, false, 1);

	initIndexes();
	nBlock = loadStructure(geomPath);
	nBlockMax = nBi * nBj * nBk;
	cout << "Number of blocks in x: " << nBi << endl;
	cout << "Number of blocks in y: " << nBj << endl;
	cout << "Number of blocks in z: " << nBk << endl;
	cout << "Total number of blocks: " << nBlock << endl;
	int nTot = 0;
	for (int i = 0; i < nBlock; i++) {
		nTot += ni * nj * nk;
	}
	cout << "Total number of elements: " << nTot << endl;
	model(b.getParInt("model", par, true, 0), &ex, &ey, &ez, &w, &opp);

	// blocks initialization
	bool isMRTmodel = (stoi(par["MRT"]) ? true : false);
	blocks = new Block * [nBlockMax];
	for (int i = 0; i < nBlockMax; i++) {
		blocks[i] = 0;
	}
	bool isFirst = true;
	int indBlock;
	for (int i = 0; i < nBlock; i++) {
		indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
		switch (problemID) {
		case 0:
			blocks[indBlock] = new BlockBasic(isFirst,blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, isMRTmodel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			break;
		case 1:
			blocks[indBlock] = new BlockFreeSurface(isFirst,blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, isMRTmodel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			if (nSmoothStep == 0) {
				delete[]((BlockFreeSurface*)blocks[indBlock])->vofS;
				delete[]((BlockFreeSurface*)blocks[indBlock])->vofS1;
			}
			break;
		case 2:
			blocks[indBlock] = new BlockFreeSurfaceLagrangeParticles(isFirst,blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, isMRTmodel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			if (nSmoothStep == 0) {
				delete[]((BlockFreeSurfaceLagrangeParticles*)blocks[indBlock])->vofS;
				delete[]((BlockFreeSurfaceLagrangeParticles*)blocks[indBlock])->vofS1;
			}
			break;
		case 21:
			blocks[indBlock] = new BlockFreeSurfaceLagrangeFibers(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, isMRTmodel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			if (nSmoothStep == 0) {
				delete[]((BlockFreeSurfaceLagrangeFibers*)blocks[indBlock])->vofS;
				delete[]((BlockFreeSurfaceLagrangeFibers*)blocks[indBlock])->vofS1;
			}
			break;
		case 3:
			blocks[indBlock] = new BlockFreeSurfaceEulerParticlesFVM(isFirst,blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, isMRTmodel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			break;
		case 4:
			blocks[indBlock] = new BlockFreeSurfaceEulerParticlesLBM(isFirst,blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, isMRTmodel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			break;
		}
		isFirst = false;
	}

	switch (problemID) {
	case 2:
		particles = new LagrangeParticles(blocks, nBlock, nBlockMax, nBi, nBj, nBk, ni, nj, nk, par);
		particles->readInitParticlesDistribution(simulPath, blocks[indBlock]->Ch, blocks[indBlock]->Cu);
		break;
	case 21:
		fibers = new LagrangeFibers(blocks, nBlock, nBlockMax, nBi, nBj, nBk, ni, nj, nk, par);
		fibers->readInitFibersDistribution(simulPath, blocks[indBlock]->Ch, blocks[indBlock]->Cu);
		break;
	}
}

Mesh::~Mesh() {
	for (int i = 0; i < nBlock; i++) {
		int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
		delete blocks[indBlock];
	}
	delete blocks;

	delete[] ex;
	delete[] ey;
	delete[] ez;
	delete[] w;
	delete[] indAll;
	delete[] indIn;
	delete[] indBound;
}

void Mesh::initIndexes() {
	indAll = new int[ni * nj * nk];
	indIn = new int[(ni - 2) * (nj - 2) * (nk - 2)];
	indBound = new int[ni * nj * nk - (ni - 2) * (nj - 2) * (nk - 2)];
	int sAll = 0;
	int sIn = 0;
	int sB = 0;
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			for (int k = 0; k < nk; k++) {
				int ind = (ni * nj) * k + ni * j + i;
				indAll[sAll] = ind;
				sAll++;
				if (i > 0 && i < ni - 1 && j > 0 && j < nj - 1 && k > 0 && k < nk - 1) {
					indIn[sIn] = ind;
					sIn++;
				}
				else {
					indBound[sB] = ind;
					sB++;
				}
			}
		}
	}
}

void Mesh::initBlocks() {
	for (int i = 0; i < nBlock; i++) { // init blocks boundary
		int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
		blocks[indBlock]->initBoundary(indBlock, geomPath); // set boundary
		blocks[indBlock]->initZeroCondition();
		blocks[indBlock]->readInitCondition(indBlock, simulPath);
	}
	if (initHydroPress) initHydrostaticPress();

	for (int i = 0; i < nBlock; i++) { // init blocks boundary
		int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
		blocks[indBlock]->computeInitDistribution();
		blocks[indBlock]->computeShiftIndexes();
	}

	// wall distance
	generateDistanceFunction();

	if (problemID == 1 || problemID == 2 || problemID == 21 || problemID == 3 || problemID == 4) {
		for (int i = 0; i < nBlock; i++) { // init blocks boundary
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			((BlockFreeSurface*)blocks[indBlock])->initFlagMap();
		}
		for (int i = 0; i < nBlock; i++) { // init blocks boundary
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			((BlockFreeSurface*)blocks[indBlock])->iniVofAndMass();
		}
	}

	if (problemID == 3) {
		for (int i = 0; i < nBlock; i++) { // init blocks boundary
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlock])->readInitConcentration(indBlock, simulPath);
		}
	} else if(problemID == 4) {
		for (int i = 0; i < nBlock; i++) { // init blocks boundary
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->readInitConcentration(indBlock, simulPath);
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->computeInitDistributionScalar();
		}
	}
}

void Mesh::evalBlocks(int iter, int step) { // dodelat verzi se zahustenim

	double t = 0;
	double omega = 0.89;
	double* G0;

	switch (step) {
	case 0:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				blocks[indBlock]->propagation();
			}
		}
		break;

	case 1:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				blocks[indBlock]->makro();
				blocks[indBlock]->collision();
			}
		}
		break;

	case 2:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->massExchange();

			}
		}
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->findFilledEmptied(iter);

			}
		}

		//int surfaceLevels = 4;
		for (int j = 1; j <= 4; j++) {
#pragma omp parallel for schedule(guided)
			for (int i = 0; i < nBlock; i++) {
				int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
				if (!blocks[indBlock]->isBlockEmpty) {
					((BlockFreeSurface*)blocks[indBlock])->findSurfaceLevels(j);
				}
			}
		}
		break;

	case 3:
#pragma omp parallel for schedule(guided)            
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->nAirToFs = 0;
				((BlockFreeSurface*)blocks[indBlock])->nFluidToFs = 0;
			}
		}
		break;

	case 4:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->vicinityFilledEmptied();
			}
		}
		break;

	case 5:
		for (int s = 0; s < nSmoothStep + 1; s++) {
#pragma omp parallel for schedule(guided)
			for (int i = 0; i < nBlock; i++) {
				int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
				if (!blocks[indBlock]->isBlockEmpty) {
					((BlockFreeSurface*)blocks[indBlock])->smoothVOF(s, nSmoothStep);
				}
			}
		}
		break;

	case 6:
#pragma omp parallel for schedule(guided)           
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->calculateNormalVectors();
				((BlockFreeSurface*)blocks[indBlock])->setWetAngle();
			}
		}
		break;

	case 7:
#pragma omp parallel for schedule(guided)          
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->calculateCurvature();
				//((BlockFreeSurface*) blocks[indBlock])->calculateCurvatureConvolution();
			}
		}
		break;

	case 8:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				blocks[indBlock]->makro();
				blocks[indBlock]->collision();
				((BlockFreeSurface*)blocks[indBlock])->newFsGeometry();
			}
		}
		break;

	case 9:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->findArtefactsFilled();
			}
		}
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->findArtefactsEmptied();
			}
		}
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->destroyArtefacts();
			}
		}
		break;

	case 10:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->massRedistribution();
			}
		}
		break;

	case 105:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->calculMassTotal();
			}
		}

		massTot = 0;
		activeCellsNumTot = 0;
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				massTot += ((BlockFreeSurface*)blocks[indBlock])->massBlock;
				activeCellsNumTot += ((BlockFreeSurface*)blocks[indBlock])->activeCellsNumBlock;
			}
		}
		break;

	case 11:
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			((BlockFreeSurface*)blocks[indBlock])->isBlockTouched = false;
		}
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurface*)blocks[indBlock])->checkIfBlockEmpty();
				((BlockFreeSurface*)blocks[indBlock])->touchNeighBlocks();
			}
		}
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (blocks[indBlock]->isBlockTouched) {
				blocks[indBlock]->isBlockEmpty = false;
			}
		}
		break;


		// Lagrangian particles
	case 20:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			((BlockFreeSurfaceLagrangeParticles*)blocks[indBlock])->computeGradP();
		}
		particles->boxParticles();
		particles->interpolateFluidVelocity();
		for (int s = 0; s < nIntIterParticle; s++) {
			particles->particleParticleInteraction();
			particles->particleWallInteraction();
			particles->particleBarrierForce(iter);
			particles->newParticlePosition(1.0 / nIntIterParticle);
		}
		particles->applyForce();
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			((BlockFreeSurfaceLagrangeParticles*)blocks[indBlock])->computeConcentrationAndFluidFraction();
		}
		break;

		// Lagrangian fibers
	case 21:
		//fibers->boxFibers();
		fibers->interpolateFluidVelocity();
		for (int s = 0; s < nIntIterParticle; s++) {
			//fibers->fiberWallInteraction();
			fibers->newFiberPosition(1.0 / nIntIterParticle);
		}
		fibers->applyForce();
		//for (int i = 0; i < nBlock; i++) {
		//	int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
		//	((BlockFreeSurfaceLagrangeFibers*)blocks[indBlock])->computeConcentrationAndFluidFraction();
		//}
		break;

		// Eulerian particles with FVM
	case 30:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlock])->computeGradP();
				((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlock])->interpolateFluidVelocity();
				((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlock])->computeParticlePressure();
			}
		}
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlock])->computeParticlePressureGradient();
			}
		}
		break;

	case 31:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlock])->computeNewParticleVelocityField();
			}
		}
		break;

	case 32:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlock])->smoothForce();
				((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlock])->upwind();
			}
		}
		break;

		// Eulerian particles with LBM
	case 40:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			//if (!blocks[indBlock]->isBlockEmpty) {
				//((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->computeGradP();
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->computeGradU();
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->interpolateFluidVelocity();
			//}
		}
		break;

	case 41:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			//if (!blocks[indBlock]->isBlockEmpty) {
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->computeGradVOFp();
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->computeNewParticleVelocityField(blocks[indBlock]->wallDistance, blocks[indBlock]->wallDistanceGradient);
			//}
		}
		break;

	case 42:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			//if (!blocks[indBlock]->isBlockEmpty) {
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->smoothForce();
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->propagationScalar();
			//}
		}
		break;
	case 43:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->setTau();
		}

#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			//if (!blocks[indBlock]->isBlockEmpty) {
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->calculConcentration();
			//}
		}

		nParticles = 0;
		maxParticlesInCell = 0;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->computeGradC();

			nParticles += int(((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->concentrationSum);
			if (int(((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->concentrationMax) > maxParticlesInCell)
				maxParticlesInCell = int(((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->concentrationMax);

		}
		break;

	case 45:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			//if (!blocks[indBlock]->isBlockEmpty) {
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->collisionScalar();
			//}
		}
		break;

	case 202:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				t = blocks[indBlock]->Ct * iter;
				G0 = ((BlockFreeSurfaceLagrangeParticles*)blocks[indBlock])->G0;
				((BlockFreeSurfaceLagrangeParticles*)blocks[indBlock])->G[0] = G0[0] * cos(omega* t) - G0[1] * sin(omega* t);
				((BlockFreeSurfaceLagrangeParticles*)blocks[indBlock])->G[1] = G0[0] * sin(omega* t) + G0[1] * cos(omega* t);
				((BlockFreeSurfaceLagrangeParticles*)blocks[indBlock])->G[2] = 0;
			}
		}

		G0 = particles->G0;
		particles->G[0] = G0[0] * cos(omega * t) - G0[1] * sin(omega * t);
		particles->G[1] = G0[0] * sin(omega * t) + G0[1] * cos(omega * t);
		particles->G[2] = 0;
		break;

	case 402:
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < nBlock; i++) {
			int indBlock = nBi * nBj * blockIndex[i][2] + nBi * blockIndex[i][1] + blockIndex[i][0];
			if (!blocks[indBlock]->isBlockEmpty) {
				t = blocks[indBlock]->Ct * iter;
				G0 = ((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->G0;
				((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->G[0] = G0[0] * cos(omega * t) - G0[1] * sin(omega * t);
				((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->G[1] = G0[0] * sin(omega * t) + G0[1] * cos(omega * t);
				((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->G[2] = 0;
			}
		}
		break;
	}
}

double Mesh::getUmax() {
	double uMax = 0;

	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int ni = blocks[bl]->ni;
		int nj = blocks[bl]->nj;
		int nk = blocks[bl]->nk;
		double* U = blocks[bl]->U;
		double* V = blocks[bl]->V;
		double* W = blocks[bl]->W;
		double Cu = blocks[bl]->Cu;
		for (int k = 0; k < nk; k++) {
			for (int j = 0; j < nj; j++) {
				for (int i = 0; i < ni; i++) {
					int ind = (ni * nj) * k + ni * j + i;
					double uAux = Cu * sqrt(U[ind] * U[ind] + V[ind] * V[ind] + W[ind] * W[ind]);
					if (uAux > uMax) {
						uMax = uAux;
					}
				}
			}
		}
	}
	return uMax;
}

void Mesh::model(int velocityModel, int** ex, int** ey, int** ez, double** w, int** opp) {
	double w0, w1, w2;
	switch (velocityModel) {
	case 29:
		nf = 9;
		*ex = new int[nf] {
			0, 1, 0, -1, 0, 1, -1, -1, 1
		};
		*ey = new int[nf] {
			0, 0, 1, 0, -1, 1, 1, -1, -1
		};
		*ez = new int[nf] {
			0, 0, 0, 0, 0, 0, 0, 0, 0
		};
		w0 = 16.0 / 36;
		w1 = 4.0 / 36;
		w2 = 1.0 / 36;
		*w = new double[nf] {
			w0, w1, w1, w1, w1, w2, w2, w2, w2
		};
		*opp = new int[nf] {
			0, 3, 4, 1, 2, 7, 8, 5, 6
		};
		break;

	case 315:
		nf = 15;
		*ex = new int[nf] {
			0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1
		};
		*ey = new int[nf] {
			0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1
		};
		*ez = new int[nf] {
			0, 0, 0, 0, 0, 1, -1, 1, 1, 1, 1, -1, -1, -1, -1
		};
		w0 = 16.0 / 72;
		w1 = 8.0 / 72;
		w2 = 1.0 / 72;
		*w = new double[nf] {
			w0, w1, w1, w1, w1, w1, w1, w2, w2, w2, w2, w2, w2, w2, w2
		};
		*opp = new int[nf] {
			0, 2, 1, 4, 3, 6, 5, 14, 13, 12, 11, 10, 9, 8, 7
		};
		break;

	case 319:
		nf = 19;
		*ex = new int[nf] {
			0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0
		};
		*ey = new int[nf] {
			0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1
		};
		*ez = new int[nf] {
			0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1
		};
		w0 = 1.0 / 3;
		w1 = 1.0 / 18;
		w2 = 1.0 / 36;
		*w = new double[nf] {
			w0, w1, w1, w1, w1, w1, w1, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2
		};
		*opp = new int[nf] {
			0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15
		};
		break;
	}
}

int Mesh::loadStructure(string geomPath) {
	int nBlocks = 0;
	string line;
	stringstream name;
	name << geomPath << "mesh/" << "blocks.txt";
	ifstream* myfile = new ifstream(name.str());

	nBi = 0;
	nBj = 0;
	nBk = 0;
	myfile->seekg(0, ios::beg);
	if (myfile->is_open()) {
		while (getline(*myfile, line)) {
			nBlocks++;
			istringstream iss(line);
			vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
			if (stoi(tokens[0]) + 1 > nBi) {
				nBi = stoi(tokens[0]) + 1;
			}
			if (stoi(tokens[1]) + 1 > nBj) {
				nBj = stoi(tokens[1]) + 1;
			}
			if (stoi(tokens[2]) + 1 > nBk) {
				nBk = stoi(tokens[2]) + 1;
			}
		}
	}
	blockMap = new int[nBi * nBj * nBk];
	for (int i = 0; i < nBi * nBj * nBk; i++) {
		blockMap[i] = -1;
	}

	blockIndex = new int* [nBlocks];
	myfile->clear();
	myfile->seekg(0, ios::beg);
	int s = 0;
	if (myfile->is_open()) {
		while (getline(*myfile, line)) {
			blockIndex[s] = new int[3];
			istringstream iss(line);
			vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
			for (int i = 0; i < 3; i++) {
				blockIndex[s][i] = stoi(tokens[i]);
			}
			int ind = (nBi * nBj) * blockIndex[s][2] + nBi * blockIndex[s][1] + blockIndex[s][0];
			blockMap[ind] = s;
			s++;
		}
		myfile->close();
	}
	delete myfile;

	return nBlocks;
}

void Mesh::writeResultsToOutFile() {
	for (int b = 0; b < nBlock; b++) {
		switch (problemID) {
			case 0:
				break;
			case 1:
				((BlockFreeSurface*)blocks[b])->writeOutCondition(b,simulPath);
				break;
			case 2:
				break;
			case 21:
				break;
			case 3:
				break;
			case 4:
				break;
		}
	}
}

void Mesh::writeResultsToVTK(int iter) {
	switch (problemID) {
	case 0:
		break;
	case 1:
		break;
	case 2:
		particles->saveParticlesPositionsToVTK(simulPath, iter);
		break;
	case 21:
		fibers->saveFibersPositionsToVTK(simulPath, iter);
		break;
	case 3:
		break;
	case 4:
		break;
	}

	Block** blocksField = blocks;

	int nPoints = 0;
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		nPoints += (blocksField[bl]->ni + 1) * (blocksField[bl]->nj + 1) * (blocksField[bl]->nk + 1);
	}

	ofstream myfile;
	int b0 = nBi * nBj * blockIndex[0][2] + nBi * blockIndex[0][1] + blockIndex[0][0];
	myfile.open(simulPath + "result/output_" + std::to_string((int)(iter * blocks[b0]->Ct * 1000000)) + "us.vtk", ios::out);
	myfile << "# vtk DataFile Version 2.6\noutput file\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	myfile << "POINTS " << nPoints << " float\n";
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		double cl = (double)blocksField[bl]->coarseLevel;
		int ni = blocksField[bl]->ni;
		int nj = blocksField[bl]->nj;
		int nk = blocksField[bl]->nk;
		for (int k = 0; k < nk + 1; k++) {
			for (int j = 0; j < nj + 1; j++) {
				for (int i = 0; i < ni + 1; i++) {
					//myfile << (i+0.5) / cl + (ni)*(globalBlockGridPosition[bl][0] - 1.0) / cl << " " << (j+0.5) / cl + (nj)*(globalBlockGridPosition[bl][1] - 1.0) / cl << " " << (k+0.5)/ cl + (nk)*(globalBlockGridPosition[bl][2] - 1) / cl << "\n";
					myfile << i / cl + ni * (blockIndex[b][0]) / cl << " " << j / cl + nj * (blockIndex[b][1]) / cl << " " << k / cl + nk * (blockIndex[b][2]) / cl << "\n";
				}
			}
		}
	}

	myfile << "\n";
	int* ei = new int[8]{
		0, 1, 1, 0, 0, 1, 1, 0
	};
	int* ej = new int[8]{
		0, 0, 1, 1, 0, 0, 1, 1
	};
	int* ek = new int[8]{
		0, 0, 0, 0, 1, 1, 1, 1
	};

	int nc = 0;
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int ni = blocksField[bl]->ni;
		int nj = blocksField[bl]->nj;
		int nk = blocksField[bl]->nk;
		bool* boundary = blocksField[bl]->boundary;
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nj; j++) {
				for (int k = 0; k < nk; k++) {
					int ind = (ni * nj) * k + ni * j + i;
					if (boundary[ind] == false) {
						nc++;
					}
				}
			}
		}
	}

	myfile << "CELLS " << nc << " " << 8 * nc + nc << "\n";
	int blInd = 0;
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int ni = blocksField[bl]->ni;
		int nj = blocksField[bl]->nj;
		int nk = blocksField[bl]->nk;
		bool* boundary = blocksField[bl]->boundary;
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nj; j++) {
				for (int k = 0; k < nk; k++) {
					int ind = (ni * nj) * k + ni * j + i;
					if (boundary[ind] == false) {
						myfile << "8 ";
						for (int p = 0; p < 8; p++) {
							int ind = ((ni + 1) * (nj + 1)) * (k + ek[p]) + (ni + 1) * (j + ej[p]) + i + ei[p] + blInd;
							myfile << ind << " ";
						}
						myfile << "\n";
					}
				}
			}
		}
		blInd += (ni + 1) * (nj + 1) * (nk + 1);
	}
	myfile << "\n";

	myfile << "CELL_TYPES " << nc << "\n";
	for (int i = 0; i < nc; i++) {
		myfile << 12 << "\n";
	}
	myfile << "\n";

	myfile << "CELL_DATA " << nc << "\n";
	myfile << "VECTORS rychlost float\n";
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int ni = blocksField[bl]->ni;
		int nj = blocksField[bl]->nj;
		int nk = blocksField[bl]->nk;
		bool* boundary = blocksField[bl]->boundary;
		double* U = blocks[bl]->U;
		double* V = blocksField[bl]->V;
		double* W = blocksField[bl]->W;
		double Cu = blocksField[bl]->Cu;
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nj; j++) {
				for (int k = 0; k < nk; k++) {
					int ind = (ni * nj) * k + ni * j + i;
					if (boundary[ind] == false) {
						myfile << (Cu * U[ind]) << " " << (Cu * V[ind]) << " " << (Cu * W[ind]) << "\n";
					}
				}
			}
		}
	}
	myfile << "\n";

	myfile << "SCALARS relative_pressure float 1\n";
	myfile << "LOOKUP_TABLE default\n";
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int ni = blocksField[bl]->ni;
		int nj = blocksField[bl]->nj;
		int nk = blocksField[bl]->nk;
		bool* boundary = blocksField[bl]->boundary;
		double* Rho = blocksField[bl]->Rho;
		double Cp = blocksField[bl]->Cp;
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nj; j++) {
				for (int k = 0; k < nk; k++) {
					int ind = (ni * nj) * k + ni * j + i;
					if (boundary[ind] == false) {
						if (Rho[ind] == 0) {
							myfile << 0 << "\n";
						}
						else {
							myfile << Cp * (Rho[ind] - 1) << "\n";
						}
					}
				}
			}
		}
	}
	myfile << "\n";

	if (iter == 0) { // wall distance is written only to first file
		myfile << "SCALARS wallDistance float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* wallDistance = blocksField[bl]->wallDistance;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << wallDistance[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "VECTORS wallDistanceGradient float\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* wallDistanceGradient = ((BlockFreeSurface*)blocksField[bl])->wallDistanceGradient;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << wallDistanceGradient[3 * ind] << " " << wallDistanceGradient[3 * ind + 1] << " " << wallDistanceGradient[3 * ind + 2] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";
	}

	if (problemID == 2) {
		particles->boxParticles();
		myfile << "SCALARS concentration float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* concentration = ((BlockFreeSurfaceLagrangeParticles*)blocksField[bl])->concentration;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << concentration[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "SCALARS particle_volume_fraction float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* particleFraction = ((BlockFreeSurfaceLagrangeParticles*)blocksField[bl])->particleFraction;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << particleFraction[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";
	}

	if (problemID == 3) {
		myfile << "SCALARS concentration float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* concentration = ((BlockFreeSurfaceEulerParticlesFVM*)blocksField[bl])->concentration;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << concentration[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "SCALARS particle_pressure float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* PP = ((BlockFreeSurfaceEulerParticlesFVM*)blocksField[bl])->PP;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << PP[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "VECTORS particle_pressure_gradient float\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* dPPx = ((BlockFreeSurfaceEulerParticlesFVM*)blocksField[bl])->dPPx;
			double* dPPy = ((BlockFreeSurfaceEulerParticlesFVM*)blocksField[bl])->dPPy;
			double* dPPz = ((BlockFreeSurfaceEulerParticlesFVM*)blocksField[bl])->dPPz;
			double Cu = blocksField[bl]->Cu;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << dPPx[ind] << " " << dPPy[ind] << " " << dPPz[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "VECTORS particle_velocity float\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* UP = ((BlockFreeSurfaceEulerParticlesFVM*)blocksField[bl])->UP;
			double* VP = ((BlockFreeSurfaceEulerParticlesFVM*)blocksField[bl])->VP;
			double* WP = ((BlockFreeSurfaceEulerParticlesFVM*)blocksField[bl])->WP;
			double Cu = blocksField[bl]->Cu;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << Cu * UP[ind] << " " << Cu * VP[ind] << " " << Cu * WP[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";
	}

	if (problemID == 4) {
		myfile << "SCALARS concentration float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* concentration = ((BlockFreeSurfaceEulerParticlesLBM*)blocksField[bl])->concentration;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << concentration[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "SCALARS user_def float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* user_def = ((BlockFreeSurfaceEulerParticlesLBM*)blocksField[bl])->user_def;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << user_def[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "VECTORS particle_velocity float\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* UP = ((BlockFreeSurfaceEulerParticlesLBM*)blocksField[bl])->UP;
			double* VP = ((BlockFreeSurfaceEulerParticlesLBM*)blocksField[bl])->VP;
			double* WP = ((BlockFreeSurfaceEulerParticlesLBM*)blocksField[bl])->WP;
			double Cu = blocksField[bl]->Cu;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << Cu * UP[ind] << " " << Cu * VP[ind] << " " << Cu * WP[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";
	}

	if (problemID == 1 || problemID == 2 || problemID == 21 || problemID == 3 || problemID == 4) {
		myfile << "VECTORS normal float\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* normal = ((BlockFreeSurface*)blocksField[bl])->normal;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << normal[3 * ind] << " " << normal[3 * ind + 1] << " " << normal[3 * ind + 2] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "SCALARS VOF float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* vof = ((BlockFreeSurface*)blocksField[bl])->vof;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << vof[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "SCALARS VOFprint float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			double* vof = ((BlockFreeSurface*)blocksField[bl])->vof;
			bool* boundary = blocksField[bl]->boundary;
			for (int i = 0; i < ni * nj * nk; i++) {
				double coeff = 1;
				int ind = indAll[i];
				for (int p = 0; p < nf; p++) {
					int nAll = blocksField[bl]->nAll;
					int indBlockNeigh = blocksField[bl]->bShift[nAll * p + ind];
					int indNeigh = blocksField[bl]->shift[nAll * p + ind];
					if ((!(indBlockNeigh > -1 && indBlockNeigh < nBlockMax && blocks[indBlockNeigh]) || blocks[bl]->boundary[ind])) {
						coeff = 0.5;
					}
				}
				if (boundary[ind] == false) {
					myfile << (coeff * vof[ind]) << "\n";
				}
			}
		}
		myfile << "\n";

		myfile << "SCALARS curvature float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			double* curvature = ((BlockFreeSurface*)blocksField[bl])->curvature;
			double Ch = blocksField[bl]->Ch;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << curvature[ind] / Ch << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "SCALARS FLAGMAP float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			int* flagMap = ((BlockFreeSurface*)blocksField[bl])->flagMap;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << flagMap[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		myfile << "SCALARS surfaceLevel float 1\n";
		myfile << "LOOKUP_TABLE default\n";
		for (int b = 0; b < nBlock; b++) {
			int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			int ni = blocksField[bl]->ni;
			int nj = blocksField[bl]->nj;
			int nk = blocksField[bl]->nk;
			bool* boundary = blocksField[bl]->boundary;
			int* surfaceLevels = ((BlockFreeSurface*)blocksField[bl])->surfaceLevels;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (boundary[ind] == false) {
							myfile << surfaceLevels[ind] << "\n";
						}
					}
				}
			}
		}
		myfile << "\n";

		if (viscosityModel) {
			myfile << "SCALARS shearRate float 1\n";
			myfile << "LOOKUP_TABLE default\n";
			for (int b = 0; b < nBlock; b++) {
				int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
				int ni = blocksField[bl]->ni;
				int nj = blocksField[bl]->nj;
				int nk = blocksField[bl]->nk;
				bool* boundary = blocksField[bl]->boundary;
				double* shearRate = ((BlockFreeSurface*)blocksField[bl])->shearRate;
				for (int i = 0; i < ni; i++) {
					for (int j = 0; j < nj; j++) {
						for (int k = 0; k < nk; k++) {
							int ind = (ni * nj) * k + ni * j + i;
							if (boundary[ind] == false) {
								myfile << shearRate[ind] << "\n";
							}
						}
					}
				}
			}
			myfile << "\n";

			myfile << "SCALARS viscosity float 1\n";
			myfile << "LOOKUP_TABLE default\n";
			for (int b = 0; b < nBlock; b++) {
				int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
				int ni = blocksField[bl]->ni;
				int nj = blocksField[bl]->nj;
				int nk = blocksField[bl]->nk;
				bool* boundary = blocksField[bl]->boundary;
				double* viscosity = ((BlockFreeSurface*)blocksField[bl])->viscosity;
				for (int i = 0; i < ni; i++) {
					for (int j = 0; j < nj; j++) {
						for (int k = 0; k < nk; k++) {
							int ind = (ni * nj) * k + ni * j + i;
							if (boundary[ind] == false) {
								myfile << viscosity[ind] << "\n";
							}
						}
					}
				}
			}
			myfile << "\n";
		}
	}
	myfile.close();

	// ===========================================================================
	ofstream myfileAux;
	myfileAux.open(simulPath + "result/init.txt", ios::out);
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int ni = blocksField[bl]->ni;
		int nj = blocksField[bl]->nj;
		int nk = blocksField[bl]->nk;
		double* rho = ((BlockFreeSurface*)blocksField[bl])->Rho;

		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nj; j++) {
				for (int k = 0; k < nk; k++) {
					int ind = (ni * nj) * k + ni * j + i;
					myfileAux << rho[ind] << " 0 0 0 1" << "\n";
				}
			}
		}
	}
}

void Mesh::writeBoundaryToVTK() {

	int nPoints = 0;
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		nPoints += (blocks[bl]->ni + 1) * (blocks[bl]->nj + 1) * (blocks[bl]->nk + 1);
	}

	ofstream myfile;
	myfile.open(simulPath + "result/boundary.vtk", ios::out);
	myfile << "# vtk DataFile Version 2.6\noutput file\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	myfile << "POINTS " << nPoints << " float\n";
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		double cl = (double)blocks[bl]->coarseLevel;
		int ni = blocks[bl]->ni;
		int nj = blocks[bl]->nj;
		int nk = blocks[bl]->nk;
		for (int k = 0; k < nk + 1; k++) {
			for (int j = 0; j < nj + 1; j++) {
				for (int i = 0; i < ni + 1; i++) {
					//myfile << (i+0.5) / cl + (ni)*(globalBlockGridPosition[bl][0] - 1.0) / cl << " " << (j+0.5) / cl + (nj)*(globalBlockGridPosition[bl][1] - 1.0) / cl << " " << (k+0.5)/ cl + (nk)*(globalBlockGridPosition[bl][2] - 1) / cl << "\n";
					//myfile << i / cl + ni * (blockIndex[b][0] - 1.0) / cl << " " << j / cl + nj * (blockIndex[b][1] - 1.0) / cl << " " << k / cl + nk * (blockIndex[b][2] - 1) / cl << "\n";
					myfile << i / cl + ni * (blockIndex[b][0]) / cl << " " << j / cl + nj * (blockIndex[b][1]) / cl << " " << k / cl + nk * (blockIndex[b][2]) / cl << "\n";
				}
			}
		}
	}
	myfile << "\n";
	int* ei = new int[8]{
		0, 1, 1, 0, 0, 1, 1, 0
	};
	int* ej = new int[8]{
		0, 0, 1, 1, 0, 0, 1, 1
	};
	int* ek = new int[8]{
		0, 0, 0, 0, 1, 1, 1, 1
	};

	int nc = 0;
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int ni = blocks[bl]->ni;
		int nj = blocks[bl]->nj;
		int nk = blocks[bl]->nk;
		bool* boundary = blocks[bl]->boundary;
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nj; j++) {
				for (int k = 0; k < nk; k++) {
					int ind = (ni * nj) * k + ni * j + i;
					if (boundary[ind] == true) {
						nc++;
					}
				}
			}
		}
	}

	myfile << "CELLS " << nc << " " << 8 * nc + nc << "\n";
	int blInd = 0;
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int ni = blocks[bl]->ni;
		int nj = blocks[bl]->nj;
		int nk = blocks[bl]->nk;
		bool* boundary = blocks[bl]->boundary;
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nj; j++) {
				for (int k = 0; k < nk; k++) {
					int ind = (ni * nj) * k + ni * j + i;
					if (boundary[ind] == true) {
						myfile << "8 ";
						for (int p = 0; p < 8; p++) {
							int ind = ((ni + 1) * (nj + 1)) * (k + ek[p]) + (ni + 1) * (j + ej[p]) + i + ei[p] + blInd;
							myfile << ind << " ";
						}
						myfile << "\n";
					}
				}
			}
		}
		blInd += (ni + 1) * (nj + 1) * (nk + 1);
	}
	myfile << "\n";

	myfile << "CELL_TYPES " << nc << "\n";
	for (int i = 0; i < nc; i++) {
		myfile << 12 << "\n";
	}
	myfile << "\n";
	myfile.close();
}

void Mesh::writeFreeSurfaceToVTK(int iter) {

	int n;
	double isolevel = 0.5;
	GRIDCELL grid;
	TRIANGLE triangles[10];
	//TRIANGLE tri[15000];
	TRIANGLE* tri = NULL;
	int ntri = 0;
	double minCoordX = 1e5;
	double maxCoordX = 0;

	int* ei = new int[8]{
		0, 1, 1, 0, 0, 1, 1, 0
	};
	int* ej = new int[8]{
		0, 0, 1, 1, 0, 0, 1, 1
	};
	int* ek = new int[8]{
		0, 0, 0, 0, 1, 1, 1, 1
	};

	Block** blocksField = blocks;
	int indBlockNeigh, indLocalNeigh;
	//cout << "Polygonising data ...\n" << endl;
	for (int b = 0; b < nBlock; b++) {
		int bl = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		double cl = (double)blocksField[bl]->coarseLevel;
		int ni = blocksField[bl]->ni;
		int nj = blocksField[bl]->nj;
		int nk = blocksField[bl]->nk;
		double* vof = ((BlockFreeSurface*)blocksField[bl])->vof;
		int* flagMap = ((BlockFreeSurface*)blocksField[bl])->flagMap;
		for (int k = 0; k < nk; k++) {
			int kg = k + blockIndex[b][2] * nk;
			for (int j = 0; j < nj; j++) {
				int jg = j + blockIndex[b][1] * nj;
				for (int i = 0; i < ni; i++) {
					int ig = i + blockIndex[b][0] * ni;
					int ind = (ni * nj) * k + ni * j + i;
					bool isSurface = false;
					for (int v = 0; v < 8; v++) {
						grid.val[v] = vof[ind];
						grid.p[v].x = ig + ei[v];
						grid.p[v].y = jg + ej[v];
						grid.p[v].z = kg + ek[v];

						double vofAvg = 0.0;
						for (int p = 0; p < 8; p++) {
							blocksField[bl]->local2globalShift(blocksField[bl]->indBlock, ind, indBlockNeigh, indLocalNeigh, ei[v] - 1 + ei[p], ej[v] - 1 + ej[p], ek[v] - 1 + ek[p]);
							if (indBlockNeigh > -1 && indBlockNeigh < nBlockMax && blocks[indBlockNeigh]) {
								if (!blocksField[indBlockNeigh]->boundary[indLocalNeigh]) {
									vofAvg += ((BlockFreeSurface*)blocksField[indBlockNeigh])->vof[indLocalNeigh];
								}
								else {
									vofAvg = isolevel - 0.01;
									break;
								}
							}
							else {
								vofAvg = isolevel - 0.01;
								break;
							}
						}
						grid.val[v] = vofAvg / 8;

						if (grid.val[v] > 0 || grid.val[v] < 1) {
							isSurface = true;
						}
					}
					if (isSurface) {
						n = Polygonise(grid, isolevel, triangles);
						tri = (TRIANGLE*)realloc(tri, (ntri + n) * sizeof(TRIANGLE));
						for (int l = 0; l < n; l++) {
							tri[ntri + l] = triangles[l];

							/*for (int ver = 0; ver < 3; ver++) {
								if (tri[ntri + l].p[ver].x < minCoordX) {
									minCoordX = tri[ntri + l].p[0].x;
								}
								if (tri[ntri + l].p[ver].x > maxCoordX) {
									maxCoordX = tri[ntri + l].p[0].x;
								}
							}*/

						}
						ntri += n;
					}
				}
			}
		}
	}
	//cout << "Total of triangles\n" << ntri << endl;

	ofstream myfile;
	int b0 = nBi * nBj * blockIndex[0][2] + nBi * blockIndex[0][1] + blockIndex[0][0];
	myfile.open(simulPath + "result/outputFS_" + std::to_string((int)(iter * blocks[b0]->Ct * 1000000)) + "us.vtk", ios::out);
	myfile << "# vtk DataFile Version 2.6\noutput file\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	myfile << "POINTS " << ntri * 3 << " float\n";
	for (int t = 0; t < ntri; t++) {
		myfile << tri[t].p[0].x << " " << tri[t].p[0].y << " " << tri[t].p[0].z << endl;
		myfile << tri[t].p[1].x << " " << tri[t].p[1].y << " " << tri[t].p[1].z << endl;
		myfile << tri[t].p[2].x << " " << tri[t].p[2].y << " " << tri[t].p[2].z << endl;
	}
	myfile << "\n ";
	myfile << "CELLS " << ntri << " " << (3 + 1) * ntri << endl;
	for (int p = 0; p < 3 * ntri; p += 3) {
		myfile << 3 << " " << p << " " << p + 1 << " " << p + 2 << endl;
	}
	myfile << "\n ";
	myfile << "CELL_TYPES " << ntri << endl;
	for (int t = 0; t < ntri; t++) {
		myfile << 5 << endl;
	}
	myfile.close();

	/*ofstream myfile2;
	myfile2.open(simulPath + "result/diameter", ios::out | ios::app);
	myfile2 << iter << " " << (maxCoordX - minCoordX) * 3.0303e-5 << "\n";
	myfile2.close();*/

	delete[] tri;
	delete[] ei;
	delete[] ej;
	delete[] ek;
}

void Mesh::generateDistanceFunction()
{
	// distance function
	cout << "generating distance function...";
	for (int b = 0; b < nBlock; b++) {
		int indBlock = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		double* wallDistance = blocks[indBlock]->wallDistance;
		for (int i = 0; i < blocks[indBlock]->nAll; i++) {
			int ind = indAll[i];
			wallDistance[ind] = 1e10;
			if (blocks[indBlock]->boundary[ind]) {
				wallDistance[ind] = -0.5;
			}
		}
	}

	bool change = true;
	int indBlockShift = 0;
	int indShift = 0;
	while (change) {
		change = false;
		for (int b = 0; b < nBlock; b++) {
			int indBlock = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
			for (int i = 0; i < blocks[indBlock]->nAll; i++) {
				int ind = indAll[i];
				if (!blocks[indBlock]->boundary[ind]) {
					for (int dx = -1; dx < 2; dx++) {
						for (int dy = -1; dy < 2; dy++) {
							for (int dz = -1; dz < 2; dz++) {
								blocks[indBlock]->local2globalShift(indBlock, ind, indBlockShift, indShift, dx, dy, dz);
								if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
									double L = sqrt(dx * dx + dy * dy + dz * dz);
									if (blocks[indBlockShift]->wallDistance[indShift] + L < blocks[indBlock]->wallDistance[ind]) {
										blocks[indBlock]->wallDistance[ind] = blocks[indBlockShift]->wallDistance[indShift] + L;
										change = true;
									}
								}
								else {
									double L = sqrt(dx * dx + dy * dy + dz * dz) / 2;
									if (L < blocks[indBlock]->wallDistance[ind]) {
										blocks[indBlock]->wallDistance[ind] = L;
										change = true;
									}
								}
							}
						}
					}
				}

			}
		}
	}

	// gradient of distance function
	double wallDist[27];
	for (int b = 0; b < nBlock; b++) {
		int indBlock = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int* bShift = blocks[indBlock]->bShift;
		int* shift = blocks[indBlock]->shift;
		int nAll = blocks[indBlock]->nAll;
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			for (int dir = 1; dir < nf; dir++) {
				indBlockShift = bShift[nAll * dir + ind];
				indShift = shift[nAll * dir + ind];
				bool dirExist = false;
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (!(blocks[indBlockShift]->boundary[indShift])) {
						wallDist[dir] = ((BlockFreeSurface*)blocks[indBlockShift])->wallDistance[indShift];
						dirExist = true;
					}
				}
				if (!dirExist) {
					indBlockShift = bShift[nAll * opp[dir] + ind];
					indShift = shift[nAll * opp[dir] + ind];
					bool oposDirExist = false;
					if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
						if (!(blocks[indBlockShift]->boundary[indShift])) {
							wallDist[dir] = -0.5; // 2 * blocks[indBlock]->wallDistance[ind] - ((BlockFreeSurface*)blocks[indBlockShift])->wallDistance[indShift];
							oposDirExist = true;
						}
					}
					if (!oposDirExist) {
						wallDist[dir] = -0.5; //blocks[indBlock]->wallDistance[ind];
					}
				}
			}

			double* wallDistanceGradient = blocks[indBlock]->wallDistanceGradient;
			wallDistanceGradient[3 * ind] = 0;
			wallDistanceGradient[3 * ind + 1] = 0;
			wallDistanceGradient[3 * ind + 2] = 0;
			double sumW = 3;
			for (int dir = 1; dir < nf; dir++) {
				wallDistanceGradient[3 * ind] += sumW * w[dir] * ex[dir] * wallDist[dir];
				wallDistanceGradient[3 * ind + 1] += sumW * w[dir] * ey[dir] * wallDist[dir];
				wallDistanceGradient[3 * ind + 2] += sumW * w[dir] * ez[dir] * wallDist[dir];
			}
		}

	}

	// wall distance polynome
	double* Ai = new double[64]
	{ 1, 0, 0, 0, 0, 0, 0, 0,
	 -1, 1, 0, 0, 0, 0, 0, 0,
	 -1, 0, 0, 1, 0, 0, 0, 0,
	 -1, 0, 0, 0, 1, 0, 0, 0,
	  1,-1, 1,-1, 0, 0, 0, 0,
	  1, 0, 0,-1,-1, 0, 0, 1,
	  1,-1, 0, 0,-1, 1, 0, 0,
	 -1, 1,-1, 1, 1,-1, 1,-1 };

	double dist[8];
	int* vx = new int[8]{ -1, 1, 1,-1,-1, 1, 1,-1 };
	int* vy = new int[8]{ -1,-1, 1, 1,-1,-1, 1, 1 };
	int* vz = new int[8]{ -1,-1,-1,-1, 1, 1, 1, 1 };
	for (int b = 0; b < nBlock; b++) {
		int indBlock = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		int nAll = blocks[indBlock]->nAll;
		for (int i = 0; i < nAll; i++) {
			int ind = indAll[i];
			for (int dir = 0; dir < 8; dir++) {
				dist[dir] = 0;
				blocks[indBlock]->local2globalShift(indBlock, ind, indBlockShift, indShift, vx[dir], vy[dir], vz[dir]);
				if (indBlockShift > -1 && indBlockShift < nBlockMax && blocks[indBlockShift]) {
					if (!(blocks[indBlockShift]->boundary[indShift])) {
						double wallDistance = blocks[indBlock]->wallDistance[ind];
						double wallDistanceShift = blocks[indBlockShift]->wallDistance[indShift];
						dist[dir] = (wallDistance + wallDistanceShift) / 2;
					}
				}
			}
			for (int p = 0; p < 8; p++) {
				double sum = 0;
				for (int q = 0; q < 8; q++) {
					sum += Ai[8 * p + q] * dist[q];
				}
				blocks[indBlock]->wallDistancePolynome[8 * ind + p] = sum;
			}
		}
	}

	delete[] Ai;
	delete[] vx;
	delete[] vy;
	delete[] vz;
	
	cout << " done" << endl;
}

void Mesh::initHydrostaticPress() {
	int ijkg[3];
	int ijkgmin[3];
	double scalarmin = 1e10;
	for (int b = 0; b < nBlock; b++) {
		int indBlock = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		if (!blocks[indBlock]->isBlockEmpty) {
			int ni = blocks[indBlock]->ni;
			int nj = blocks[indBlock]->nj;
			int nk = blocks[indBlock]->nk;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						bool lo = !blocks[indBlock]->boundary[ind];
						if (problemID == 1 || problemID == 2 || problemID == 3 || problemID == 4) {
							bool lo = ((BlockFreeSurface*)blocks[indBlock])->flagMap[ind] == 1;
						}
						if (lo) {
							ijkg[0] = blockIndex[b][0] * ni + i;
							ijkg[1] = blockIndex[b][1] * nj + j;
							ijkg[2] = blockIndex[b][2] * nk + k;
							double scalar = ijkg[0] * blocks[indBlock]->G[0] + ijkg[1] * blocks[indBlock]->G[1] + ijkg[2] * blocks[indBlock]->G[2];
							if (scalar < scalarmin) {
								scalarmin = scalar;
								ijkgmin[0] = ijkg[0];
								ijkgmin[1] = ijkg[1];
								ijkgmin[2] = ijkg[2];
							}
						}
					}
				}
			}
		}
	}
	for (int b = 0; b < nBlock; b++) {
		int indBlock = nBi * nBj * blockIndex[b][2] + nBi * blockIndex[b][1] + blockIndex[b][0];
		if (!blocks[indBlock]->isBlockEmpty) {
			int ni = blocks[indBlock]->ni;
			int nj = blocks[indBlock]->nj;
			int nk = blocks[indBlock]->nk;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						bool lo = !blocks[indBlock]->boundary[ind];
						if (problemID == 1 || problemID == 2 || problemID == 3 || problemID == 4) {
							bool lo = ((BlockFreeSurface*)blocks[indBlock])->flagMap[ind] == 1;
						}
						if (lo) {
							ijkg[0] = blockIndex[b][0] * ni + i;
							ijkg[1] = blockIndex[b][1] * nj + j;
							ijkg[2] = blockIndex[b][2] * nk + k;
							double absG = sqrt(pow(blocks[indBlock]->G[0], 2) + pow(blocks[indBlock]->G[1], 2) + pow(blocks[indBlock]->G[2], 2));
							double h = ((ijkg[0] - ijkgmin[0]) * blocks[indBlock]->G[0] + (ijkg[1] - ijkgmin[1]) * blocks[indBlock]->G[1] + (ijkg[2] - ijkgmin[2]) * blocks[indBlock]->G[2]) / absG;
							blocks[indBlock]->Rho[ind] += h * absG * 3.0;
							//cout << h * absG * 3.0 << endl;
							//cout << h * blocks[indBlock]->Ch * 9.81 * 1000 / blocks[indBlock]->Cp << endl;
							blocks[indBlock]->dRho[ind] = blocks[indBlock]->Rho[ind] - blocks[indBlock]->rho0;
						}
					}
				}
			}
		}
	}
}

double Mesh::getCt() {
	int b0 = nBi * nBj * blockIndex[0][2] + nBi * blockIndex[0][1] + blockIndex[0][0];

	return blocks[b0]->Ct;
}