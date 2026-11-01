/*
 * File:   Mesh.cpp
 * Author: obublik
 */

#include "Mesh.h"
#include "MarchingCubes.h"
#include "Block.h"
#include <queue>
#include <array>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <limits>

using namespace std;

#define FSB(b) ((BlockFreeSurface*)(b))

Mesh::Mesh(string geomPath, string simulPath, const ParamMap& par) {
	this->geomPath = geomPath;
	this->simulPath = simulPath;
	particles = nullptr;
	fibers = nullptr;
	massTot = 0;
	nParticles = 0;
	maxParticlesInCell = 0;
	activeCellsNumTot = 0;
	ni = Block::getParInt("ni", par, true, 0);
	nj = Block::getParInt("nj", par, true, 0);
	nk = Block::getParInt("nk", par, true, 0);
	problemID = Block::getParInt("problemID", par, true, 0);
	viscosityModel = Block::getParInt("viscosityModel", par, false, 0);
	initHydroPress = Block::getParInt("initializeWithHydroPress", par, false, 0) != 0;
	nSmoothStep = min(max(Block::getParInt("nSmoothStep", par, false, 2), 0), 2);
	nIntIterParticle = Block::getParInt("nIntIter", par, false, 1);
	gravityRotationOmega = Block::getParDouble("gravityRotationOmega", par, false, 0.89);
	vtkBinary = Block::getParInt("vtkBinary", par, false, 1) != 0;
	writeFields = Block::getParInt("writeFields", par, false, 1) != 0;
	if (Block::hasPar("collisionModel", par)) {
		collisionModel = Block::getParInt("collisionModel", par, false, 0);
	}
	else {
		collisionModel = Block::getParInt("MRT", par, false, 0) ? COLLISION_MRT : COLLISION_BGK;
	}
	if (collisionModel < 0 || collisionModel > 2) {
		throw invalid_argument("collisionModel has to be 0 (BGK), 1 (MRT) or 2 (TRT)");
	}

	model(Block::getParInt("model", par, true, 0), &ex, &ey, &ez, &w, &opp);
	if (nf == 9 && nk != 1) {
		throw invalid_argument("D2Q9 model requires nk = 1");
	}
	if (ni < 1 || nj < 1 || nk < 1) {
		throw invalid_argument("block dimensions ni, nj, nk have to be positive");
	}

	initIndexes();
	nBlock = loadStructure(geomPath);
	nBlockMax = nBi * nBj * nBk;
	cout << "Number of blocks in x: " << nBi << endl;
	cout << "Number of blocks in y: " << nBj << endl;
	cout << "Number of blocks in z: " << nBk << endl;
	cout << "Total number of blocks: " << nBlock << endl;
	cout << "Total number of elements: " << (long long)nBlock * ni * nj * nk << endl;

	// blocks initialization
	blocks = new Block * [nBlockMax];
	for (int i = 0; i < nBlockMax; i++) {
		blocks[i] = nullptr;
	}
	bool isFirst = true;
	for (int i = 0; i < nBlock; i++) {
		int indBlock = blockId[i];
		switch (problemID) {
		case 0:
			blocks[indBlock] = new BlockBasic(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, collisionModel, viscosityModel != 0, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			break;
		case 1:
			blocks[indBlock] = new BlockFreeSurface(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			break;
		case 2:
			blocks[indBlock] = new BlockFreeSurfaceLagrangeParticles(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			break;
		case 21:
			blocks[indBlock] = new BlockFreeSurfaceLagrangeFibers(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			break;
		case 3:
			blocks[indBlock] = new BlockFreeSurfaceEulerParticlesFVM(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			break;
		case 4:
			blocks[indBlock] = new BlockFreeSurfaceEulerParticlesLBM(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, 1, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
			break;
		default:
			throw invalid_argument("unknown problemID " + to_string(problemID));
		}
		isFirst = false;
	}

	Block* b0 = blocks[blockId[0]];
	switch (problemID) {
	case 2:
		particles = new LagrangeParticles(blocks, nBlock, nBlockMax, nBi, nBj, nBk, ni, nj, nk, par);
		particles->readInitParticlesDistribution(simulPath, b0->Ch, b0->Cu);
		break;
	case 21:
		fibers = new LagrangeFibers(blocks, nBlock, nBlockMax, nBi, nBj, nBk, ni, nj, nk, par);
		fibers->readInitFibersDistribution(simulPath, b0->Ch, b0->Cu);
		break;
	}
}

Mesh::~Mesh() {
	for (int i = 0; i < nBlock; i++) {
		delete blocks[blockId[i]];
	}
	delete[] blocks;
	for (int i = 0; i < nBlock; i++) {
		delete[] blockIndex[i];
	}
	delete[] blockIndex;
	delete[] blockId;
	delete[] blockMap;
	delete particles;
	delete fibers;

	delete[] ex;
	delete[] ey;
	delete[] ez;
	delete[] w;
	delete[] opp;
	delete[] indAll;
	delete[] indIn;
	delete[] indBound;
}

bool Mesh::isFreeSurface() const {
	return problemID == 1 || problemID == 2 || problemID == 21 || problemID == 3 || problemID == 4;
}

void Mesh::initIndexes() {
	Block::countCells(ni, nj, nk, nIn, nBound);
	indAll = new int[ni * nj * nk];
	indIn = new int[max(nIn, 1)];
	indBound = new int[max(nBound, 1)];
	int sAll = 0;
	int sIn = 0;
	int sB = 0;
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			for (int k = 0; k < nk; k++) {
				int ind = (ni * nj) * k + ni * j + i;
				indAll[sAll++] = ind;
				if (Block::isInterior(i, j, k, ni, nj, nk)) {
					indIn[sIn++] = ind;
				}
				else {
					indBound[sB++] = ind;
				}
			}
		}
	}
}

void Mesh::initBlocks() {
	for (int i = 0; i < nBlock; i++) {
		int indBlock = blockId[i];
		blocks[indBlock]->initBoundary(indBlock, geomPath); // set boundary
		blocks[indBlock]->initZeroCondition();
		blocks[indBlock]->readInitCondition(indBlock, simulPath);
	}
	for (int i = 0; i < nBlock; i++) {
		blocks[blockId[i]]->computeShiftIndexes();
	}

	if (isFreeSurface()) {
		for (int i = 0; i < nBlock; i++) {
			FSB(blocks[blockId[i]])->initFlagMap();
		}
	}
	if (initHydroPress) {
		initHydrostaticPress();
	}
	if (isFreeSurface()) {
		for (int i = 0; i < nBlock; i++) {
			FSB(blocks[blockId[i]])->iniVofAndMass();
		}
	}
	for (int i = 0; i < nBlock; i++) {
		blocks[blockId[i]]->computeInitDistribution();
	}

	// wall distance
	generateDistanceFunction();

	if (problemID == 3) {
		for (int i = 0; i < nBlock; i++) {
			int indBlock = blockId[i];
			((BlockFreeSurfaceEulerParticlesFVM*)blocks[indBlock])->readInitConcentration(indBlock, simulPath);
		}
	}
	else if (problemID == 4) {
		for (int i = 0; i < nBlock; i++) {
			int indBlock = blockId[i];
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->readInitConcentration(indBlock, simulPath);
			((BlockFreeSurfaceEulerParticlesLBM*)blocks[indBlock])->computeInitDistributionScalar();
		}
	}
}

void Mesh::evalBlocks(int iter, int step) {
	switch (step) {
	case 0:
		forBlocks([](Block* b) { b->propagation(); });
		break;

	case 1: // classic LBM: moments and collision are fused in BlockBasic::collision
		forBlocks([](Block* b) { b->collision(); });
		break;

	case 2: // mass exchange, find filled/emptied cells
		forBlocks([](Block* b) { FSB(b)->computeInterfaceTypes(); });
		forBlocks([](Block* b) { FSB(b)->massExchange(); });
		forBlocks([iter](Block* b) { FSB(b)->findFilledEmptied(iter); });
		for (int j = 1; j <= 4; j++) {
			forBlocks([j](Block* b) { FSB(b)->findSurfaceLevels(j); });
		}
		break;

	case 3:
		forBlocks([](Block* b) {
			FSB(b)->nAirToFs = 0;
			FSB(b)->nFluidToFs = 0;
		});
		break;

	case 4:
		forBlocks([](Block* b) { FSB(b)->vicinityFilledEmptied(); });
		break;

	case 5:
		for (int s = 0; s < nSmoothStep + 1; s++) {
			int ns = nSmoothStep;
			forBlocks([s, ns](Block* b) { FSB(b)->smoothVOF(s, ns); });
		}
		break;

	case 6:
		forBlocks([](Block* b) {
			FSB(b)->calculateNormalVectors();
			FSB(b)->setWetAngle();
		});
		break;

	case 7:
		forBlocks([](Block* b) { FSB(b)->calculateCurvature(); });
		break;

	case 8: // macro, collision, new free surface geometry (separate passes, neighbours read these fields)
		forBlocks([](Block* b) { b->makro(); });
		forBlocks([](Block* b) {
			b->collision();
			FSB(b)->newFsGeometryPrepare();
		});
		forBlocks([](Block* b) { FSB(b)->newFsGeometryApply(); });
		break;

	case 9:
		forBlocks([](Block* b) { FSB(b)->findArtefactsFilled(); });
		forBlocks([](Block* b) { FSB(b)->findArtefactsEmptied(); });
		forBlocks([](Block* b) { FSB(b)->destroyArtefacts(); });
		forBlocks([](Block* b) { FSB(b)->fixFluidGasContacts(); });
		break;

	case 10:
		forBlocks([](Block* b) { FSB(b)->massRedistribution(); });
		forBlocks([](Block* b) { FSB(b)->redistributeMassErr(); });
		break;

	case 105: {
		forBlocks([](Block* b) { FSB(b)->calculMassTotal(); });
		massTot = 0;
		activeCellsNumTot = 0;
		for (int i = 0; i < nBlock; i++) {
			Block* b = blocks[blockId[i]];
			if (!b->isBlockEmpty) {
				massTot += FSB(b)->massBlock;
				activeCellsNumTot += FSB(b)->activeCellsNumBlock;
			}
			else {
				massTot += FSB(b)->massErr * b->Crho * b->Ch * b->Ch * b->Ch;
			}
		}
		break;
	}

	case 11:
		for (int i = 0; i < nBlock; i++) {
			blocks[blockId[i]]->isBlockTouched = false;
		}
		forBlocks([](Block* b) {
			FSB(b)->checkIfBlockEmpty();
			FSB(b)->touchNeighBlocks();
		});
		for (int i = 0; i < nBlock; i++) {
			Block* b = blocks[blockId[i]];
			if (b->isBlockTouched) {
				b->isBlockEmpty = false;
			}
		}
		break;

		// Lagrangian particles
	case 20:
		forBlocks([](Block* b) { ((BlockFreeSurfaceLagrangeParticles*)b)->computeGradP(); }, false);
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
			((BlockFreeSurfaceLagrangeParticles*)blocks[blockId[i]])->computeConcentrationAndFluidFraction();
		}
		break;

		// Lagrangian fibers
	case 21:
		fibers->interpolateFluidVelocity();
		for (int s = 0; s < nIntIterParticle; s++) {
			fibers->newFiberPosition(1.0 / nIntIterParticle);
		}
		fibers->applyForce();
		break;

		// Eulerian particles with FVM
	case 30:
		forBlocks([](Block* b) {
			BlockFreeSurfaceEulerParticlesFVM* e = (BlockFreeSurfaceEulerParticlesFVM*)b;
			e->computeGradP();
			e->interpolateFluidVelocity();
			e->computeParticlePressure();
		});
		forBlocks([](Block* b) { ((BlockFreeSurfaceEulerParticlesFVM*)b)->computeParticlePressureGradient(); });
		break;

	case 31:
		forBlocks([](Block* b) { ((BlockFreeSurfaceEulerParticlesFVM*)b)->computeNewParticleVelocityField(); });
		break;

	case 32: // fluxes read the neighbour concentration, so they are applied in a separate pass
		forBlocks([](Block* b) {
			BlockFreeSurfaceEulerParticlesFVM* e = (BlockFreeSurfaceEulerParticlesFVM*)b;
			e->smoothForce();
			e->computeFlux();
		});
		forBlocks([](Block* b) { ((BlockFreeSurfaceEulerParticlesFVM*)b)->applyFlux(); });
		break;

		// Eulerian particles with LBM
	case 40:
		forBlocks([](Block* b) {
			((BlockFreeSurfaceEulerParticlesLBM*)b)->computeGradU();
			((BlockFreeSurfaceEulerParticlesLBM*)b)->interpolateFluidVelocity();
		}, false);
		break;

	case 41:
		forBlocks([](Block* b) {
			((BlockFreeSurfaceEulerParticlesLBM*)b)->computeGradVOFp();
			((BlockFreeSurfaceEulerParticlesLBM*)b)->computeNewParticleVelocityField(b->wallDistance, b->wallDistanceGradient);
		}, false);
		break;

	case 42:
		forBlocks([](Block* b) {
			((BlockFreeSurfaceEulerParticlesLBM*)b)->smoothForce();
			((BlockFreeSurfaceEulerParticlesLBM*)b)->propagationScalar();
		}, false);
		break;

	case 43: {
		forBlocks([](Block* b) { ((BlockFreeSurfaceEulerParticlesLBM*)b)->setTau(); }, false);
		forBlocks([](Block* b) { ((BlockFreeSurfaceEulerParticlesLBM*)b)->calculConcentration(); }, false);
		forBlocks([](Block* b) { ((BlockFreeSurfaceEulerParticlesLBM*)b)->computeGradC(); }, false);
		double sum = 0;
		double cMax = 0;
		for (int i = 0; i < nBlock; i++) {
			BlockFreeSurfaceEulerParticlesLBM* b = (BlockFreeSurfaceEulerParticlesLBM*)blocks[blockId[i]];
			sum += b->concentrationSum;
			cMax = max(cMax, b->concentrationMax);
		}
		nParticles = (int)sum;
		maxParticlesInCell = (int)cMax;
		break;
	}

	case 45:
		forBlocks([](Block* b) { ((BlockFreeSurfaceEulerParticlesLBM*)b)->collisionScalar(); }, false);
		break;

	case 202:
	case 402: { // rotating gravity
		double t = blocks[blockId[0]]->Ct * iter;
		double c = cos(gravityRotationOmega * t);
		double s = sin(gravityRotationOmega * t);
		forBlocks([c, s](Block* b) {
			b->G[0] = b->G0[0] * c - b->G0[1] * s;
			b->G[1] = b->G0[0] * s + b->G0[1] * c;
			b->G[2] = 0;
		});
		if (step == 202) {
			double* G0 = particles->G0;
			particles->G[0] = G0[0] * c - G0[1] * s;
			particles->G[1] = G0[0] * s + G0[1] * c;
			particles->G[2] = 0;
		}
		break;
	}
	}
}

// one time step of the selected problem
void Mesh::timeStep(int iter) {
	switch (problemID) {
	case 0: // classic LBM
		evalBlocks(iter, 0);
		evalBlocks(iter, 1);
		break;
	case 1: // free surface LBM
		evalBlocks(iter, 2); // mass exchange, find filled/emptied cells
		evalBlocks(iter, 3);
		evalBlocks(iter, 4);
		evalBlocks(iter, 5); // smooth
		evalBlocks(iter, 6); // normals and wet angle
		evalBlocks(iter, 7); // curvature
		evalBlocks(iter, 0); // propagation
		evalBlocks(iter, 8); // macro, collision, new fs
		evalBlocks(iter, 9); // destroy artefacts
		evalBlocks(iter, 10); // mass redistribution
		evalBlocks(iter, 11); // check empty block flag
		break;
	case 2: // free surface LBM with particles
		evalBlocks(iter, 202);
		evalBlocks(iter, 2);
		evalBlocks(iter, 3);
		evalBlocks(iter, 4);
		evalBlocks(iter, 5);
		evalBlocks(iter, 6); // normals and wet angle
		evalBlocks(iter, 7); // curvature
		evalBlocks(iter, 0); // propagation
		evalBlocks(iter, 20); // particles
		evalBlocks(iter, 8); // macro, collision, new fs
		evalBlocks(iter, 9); // destroy artefacts
		evalBlocks(iter, 10);
		evalBlocks(iter, 11);
		break;
	case 21: // free surface LBM with fibers
		evalBlocks(iter, 2);
		evalBlocks(iter, 3);
		evalBlocks(iter, 4);
		evalBlocks(iter, 5);
		evalBlocks(iter, 6); // normals and wet angle
		evalBlocks(iter, 7); // curvature
		evalBlocks(iter, 0); // propagation
		evalBlocks(iter, 21); // fibers dynamic
		evalBlocks(iter, 8); // macro, collision, new fs
		evalBlocks(iter, 9); // destroy artefacts
		evalBlocks(iter, 10);
		evalBlocks(iter, 11);
		break;
	case 3: // free surface LBM with euler particles - solve transport eq. using FVM
		evalBlocks(iter, 2);  // mass exchange, find filled/emptied cells
		evalBlocks(iter, 3);
		evalBlocks(iter, 4);
		evalBlocks(iter, 5);  // smooth
		evalBlocks(iter, 6);  // normals and wet angle
		evalBlocks(iter, 7);  // curvature
		evalBlocks(iter, 0);  // propagation
		evalBlocks(iter, 30); // compute GradP + interpolate fluid velocity
		evalBlocks(iter, 31); // compute new particle velocity
		evalBlocks(iter, 32); // smooth force + upwind
		evalBlocks(iter, 8);  // macro, collision, new fs
		evalBlocks(iter, 9);  // destroy artefacts
		evalBlocks(iter, 10); // mass redistribution
		evalBlocks(iter, 11); // check empty block flag
		break;
	case 4: // free surface LBM with euler particles - solve transport eq. using LBM
		evalBlocks(iter, 402);
		evalBlocks(iter, 2);  // mass exchange, find filled/emptied cells
		evalBlocks(iter, 3);
		evalBlocks(iter, 4);
		evalBlocks(iter, 5);  // smooth
		evalBlocks(iter, 6);  // normals and wet angle
		evalBlocks(iter, 7);  // curvature
		evalBlocks(iter, 0);  // propagation
		evalBlocks(iter, 40); // compute GradP + interpolate fluid velocity
		evalBlocks(iter, 41); // compute new particle velocity
		evalBlocks(iter, 42); // smooth force + propagation of scalar
		evalBlocks(iter, 43); // compute concentration
		evalBlocks(iter, 45); // collision of scalar
		evalBlocks(iter, 8);  // macro, collision, new fs
		evalBlocks(iter, 9);  // destroy artefacts
		evalBlocks(iter, 10); // mass redistribution
		evalBlocks(iter, 11); // check empty block flag
		break;
	}
}

double Mesh::getUmax() {
	vector<double> uMaxBlock(nBlock, 0.0);
	vector<char> nanBlock(nBlock, 0);
#pragma omp parallel for schedule(dynamic)
	for (int b = 0; b < nBlock; b++) {
		Block* bl = blocks[blockId[b]];
		double u[3];
		for (int i = 0; i < bl->nAll; i++) {
			int ind = bl->indAll[i];
			if (bl->boundary[ind]) {
				continue;
			}
			bl->outputVelocity(ind, u);
			double uAux = bl->Cu * sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
			if (std::isnan(uAux)) {
				nanBlock[b] = 1;
			}
			else if (uAux > uMaxBlock[b]) {
				uMaxBlock[b] = uAux;
			}
		}
	}
	double uMax = 0;
	for (int b = 0; b < nBlock; b++) {
		if (nanBlock[b]) {
			return numeric_limits<double>::quiet_NaN();
		}
		uMax = max(uMax, uMaxBlock[b]);
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

	default:
		throw invalid_argument("unknown velocity model " + to_string(velocityModel) + " (use 29, 315 or 319)");
	}
}

int Mesh::loadStructure(string geomPath) {
	string fileName = geomPath + "mesh/blocks.txt";
	ifstream myfile(fileName);
	if (!myfile.is_open()) {
		throw runtime_error("cannot open " + fileName);
	}
	vector<array<int, 3> > list;
	string line;
	nBi = nBj = nBk = 0;
	while (getline(myfile, line)) {
		istringstream iss(line);
		array<int, 3> idx;
		if (!(iss >> idx[0] >> idx[1] >> idx[2])) {
			continue; // empty line
		}
		if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) {
			throw runtime_error("negative block index in " + fileName);
		}
		nBi = max(nBi, idx[0] + 1);
		nBj = max(nBj, idx[1] + 1);
		nBk = max(nBk, idx[2] + 1);
		list.push_back(idx);
	}
	if (list.empty()) {
		throw runtime_error("no blocks in " + fileName);
	}
	int nBlocks = (int)list.size();
	blockMap = new int[nBi * nBj * nBk];
	for (int i = 0; i < nBi * nBj * nBk; i++) {
		blockMap[i] = -1;
	}
	blockIndex = new int* [nBlocks];
	blockId = new int[nBlocks];
	for (int s = 0; s < nBlocks; s++) {
		blockIndex[s] = new int[3]{ list[s][0], list[s][1], list[s][2] };
		int ind = (nBi * nBj) * list[s][2] + nBi * list[s][1] + list[s][0];
		if (blockMap[ind] != -1) {
			throw runtime_error("duplicate block in " + fileName);
		}
		blockMap[ind] = s;
		blockId[s] = ind;
	}
	return nBlocks;
}

void Mesh::writeResultsToOutFile() {
	if (!isFreeSurface()) {
		return;
	}
	for (int b = 0; b < nBlock; b++) {
		FSB(blocks[blockId[b]])->writeOutCondition(blockId[b], simulPath);
	}
}

string Mesh::outputName(const string& prefix, int iter, const string& ext) {
	long long us = llround(iter * blocks[blockId[0]]->Ct * 1000000);
	return prefix + std::to_string(us) + "us" + ext;
}

void Mesh::writePvd(const string& fileName, const vector<pair<double, string> >& files) {
	ofstream pvd(simulPath + "result/" + fileName);
	pvd << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
	for (auto& f : files) {
		pvd << "<DataSet timestep=\"" << f.first << "\" file=\"" << f.second << "\"/>\n";
	}
	pvd << "</Collection>\n</VTKFile>\n";
}

namespace {
	// legacy VTK writer, binary data are big endian
	class VtkStream {
	public:
		ofstream os;
		bool binary;
		VtkStream(const string& name, bool binary) : os(name, ios::out | ios::binary), binary(binary) {
			os << "# vtk DataFile Version 2.6\noutput file\n" << (binary ? "BINARY" : "ASCII") << "\n";
		}
		static void swap4(char* p) {
			swap(p[0], p[3]);
			swap(p[1], p[2]);
		}
		void floats(const vector<float>& v, int perLine) {
			if (binary) {
				vector<float> t(v);
				for (float& x : t) {
					swap4((char*)&x);
				}
				os.write((const char*)t.data(), t.size() * sizeof(float));
				os << "\n";
			}
			else {
				for (size_t i = 0; i < v.size(); i++) {
					os << v[i] << (((i + 1) % perLine == 0) ? "\n" : " ");
				}
			}
		}
		void ints(const vector<int32_t>& v, int perLine) {
			if (binary) {
				vector<int32_t> t(v);
				for (int32_t& x : t) {
					swap4((char*)&x);
				}
				os.write((const char*)t.data(), t.size() * sizeof(int32_t));
				os << "\n";
			}
			else {
				for (size_t i = 0; i < v.size(); i++) {
					os << v[i] << (((i + 1) % perLine == 0) ? "\n" : " ");
				}
			}
		}
	};

	// hexahedral cells of the selected cells of all blocks
	void writeHexaGrid(VtkStream& vtk, Block** blocks, int* blockId, int** blockIndex, int nBlock, bool selectBoundary, long long& nc) {
		vector<float> pts;
		for (int b = 0; b < nBlock; b++) {
			Block* bl = blocks[blockId[b]];
			for (int k = 0; k < bl->nk + 1; k++) {
				for (int j = 0; j < bl->nj + 1; j++) {
					for (int i = 0; i < bl->ni + 1; i++) {
						pts.push_back((float)(i + bl->ni * blockIndex[b][0]));
						pts.push_back((float)(j + bl->nj * blockIndex[b][1]));
						pts.push_back((float)(k + bl->nk * blockIndex[b][2]));
					}
				}
			}
		}
		vtk.os << "DATASET UNSTRUCTURED_GRID\nPOINTS " << pts.size() / 3 << " float\n";
		vtk.floats(pts, 3);
		vtk.os << "\n";

		int ei[8] = { 0, 1, 1, 0, 0, 1, 1, 0 };
		int ej[8] = { 0, 0, 1, 1, 0, 0, 1, 1 };
		int ek[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
		vector<int32_t> cells;
		nc = 0;
		long long blInd = 0;
		for (int b = 0; b < nBlock; b++) {
			Block* bl = blocks[blockId[b]];
			int ni = bl->ni, nj = bl->nj, nk = bl->nk;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (bl->boundary[ind] == selectBoundary) {
							cells.push_back(8);
							for (int p = 0; p < 8; p++) {
								cells.push_back((int32_t)(((ni + 1) * (nj + 1)) * (k + ek[p]) + (ni + 1) * (j + ej[p]) + i + ei[p] + blInd));
							}
							nc++;
						}
					}
				}
			}
			blInd += (ni + 1) * (nj + 1) * (nk + 1);
		}
		vtk.os << "CELLS " << nc << " " << 9 * nc << "\n";
		vtk.ints(cells, 9);
		vtk.os << "\nCELL_TYPES " << nc << "\n";
		vtk.ints(vector<int32_t>(nc, 12), 1);
		vtk.os << "\n";
	}

	// iterates over the selected cells in the order of writeHexaGrid
	template <class F> void forCells(Block** blocks, int* blockId, int nBlock, bool selectBoundary, F fn) {
		for (int b = 0; b < nBlock; b++) {
			Block* bl = blocks[blockId[b]];
			int ni = bl->ni, nj = bl->nj, nk = bl->nk;
			for (int i = 0; i < ni; i++) {
				for (int j = 0; j < nj; j++) {
					for (int k = 0; k < nk; k++) {
						int ind = (ni * nj) * k + ni * j + i;
						if (bl->boundary[ind] == selectBoundary) {
							fn(bl, ind);
						}
					}
				}
			}
		}
	}
}

void Mesh::writeResultsToVTK(int iter) {
	switch (problemID) {
	case 2:
		particles->saveParticlesPositionsToVTK(simulPath, iter);
		break;
	case 21:
		fibers->saveFibersPositionsToVTK(simulPath, iter);
		break;
	}

	string name = outputName("output_", iter, ".vtk");
	VtkStream vtk(simulPath + "result/" + name, vtkBinary);
	long long nc;
	writeHexaGrid(vtk, blocks, blockId, blockIndex, nBlock, false, nc);
	vtk.os << "CELL_DATA " << nc << "\n";

	auto scalar = [&](const string& fieldName, auto fn) {
		vector<float> v;
		v.reserve(nc);
		forCells(blocks, blockId, nBlock, false, [&](Block* b, int ind) { v.push_back((float)fn(b, ind)); });
		vtk.os << "SCALARS " << fieldName << " float 1\nLOOKUP_TABLE default\n";
		vtk.floats(v, 1);
		vtk.os << "\n";
	};
	auto vector3 = [&](const string& fieldName, auto fn) {
		vector<float> v;
		v.reserve(3 * nc);
		double a[3];
		forCells(blocks, blockId, nBlock, false, [&](Block* b, int ind) {
			fn(b, ind, a);
			v.push_back((float)a[0]);
			v.push_back((float)a[1]);
			v.push_back((float)a[2]);
		});
		vtk.os << "VECTORS " << fieldName << " float\n";
		vtk.floats(v, 3);
		vtk.os << "\n";
	};

	vector3("rychlost", [](Block* b, int ind, double* a) {
		b->outputVelocity(ind, a);
		for (int d = 0; d < 3; d++) {
			a[d] *= b->Cu;
		}
	});
	scalar("relative_pressure", [](Block* b, int ind) { return b->Rho[ind] == 0 ? 0.0 : b->Cp * (b->Rho[ind] - 1); });

	if (iter == 0) { // wall distance is written only to first file
		scalar("wallDistance", [](Block* b, int ind) { return b->wallDistance[ind]; });
		vector3("wallDistanceGradient", [](Block* b, int ind, double* a) {
			for (int d = 0; d < 3; d++) {
				a[d] = b->wallDistanceGradient[3 * ind + d];
			}
		});
	}

	if (problemID == 2) {
		particles->boxParticles();
		scalar("concentration", [](Block* b, int ind) { return ((BlockFreeSurfaceLagrangeParticles*)b)->concentration[ind]; });
		scalar("particle_volume_fraction", [](Block* b, int ind) { return ((BlockFreeSurfaceLagrangeParticles*)b)->particleFraction[ind]; });
	}

	if (problemID == 3) {
		typedef BlockFreeSurfaceEulerParticlesFVM E;
		scalar("concentration", [](Block* b, int ind) { return ((E*)b)->concentration[ind]; });
		scalar("particle_pressure", [](Block* b, int ind) { return ((E*)b)->PP[ind]; });
		vector3("particle_pressure_gradient", [](Block* b, int ind, double* a) {
			a[0] = ((E*)b)->dPPx[ind];
			a[1] = ((E*)b)->dPPy[ind];
			a[2] = ((E*)b)->dPPz[ind];
		});
		vector3("particle_velocity", [](Block* b, int ind, double* a) {
			a[0] = b->Cu * ((E*)b)->UP[ind];
			a[1] = b->Cu * ((E*)b)->VP[ind];
			a[2] = b->Cu * ((E*)b)->WP[ind];
		});
	}

	if (problemID == 4) {
		typedef BlockFreeSurfaceEulerParticlesLBM E;
		scalar("concentration", [](Block* b, int ind) { return ((E*)b)->concentration[ind]; });
		scalar("user_def", [](Block* b, int ind) { return ((E*)b)->user_def[ind]; });
		vector3("particle_velocity", [](Block* b, int ind, double* a) {
			a[0] = b->Cu * ((E*)b)->UP[ind];
			a[1] = b->Cu * ((E*)b)->VP[ind];
			a[2] = b->Cu * ((E*)b)->WP[ind];
		});
	}

	if (isFreeSurface()) {
		vector3("normal", [](Block* b, int ind, double* a) {
			for (int d = 0; d < 3; d++) {
				a[d] = FSB(b)->normal[3 * ind + d];
			}
		});
		scalar("VOF", [](Block* b, int ind) { return FSB(b)->vof[ind]; });
		int nBlockMax = this->nBlockMax;
		scalar("VOFprint", [nBlockMax](Block* b, int ind) {
			double coeff = 1;
			for (int p = 0; p < b->nf; p++) {
				int bl = b->nbBlock(ind, p);
				if (bl < 0 || bl >= nBlockMax || b->blocks[bl] == nullptr || b->boundary[ind]) {
					coeff = 0.5;
				}
			}
			return coeff * FSB(b)->vof[ind];
		});
		scalar("curvature", [](Block* b, int ind) { return FSB(b)->curvature[ind] / b->Ch; });
		scalar("FLAGMAP", [](Block* b, int ind) { return (double)FSB(b)->flagMap[ind]; });
		scalar("surfaceLevel", [](Block* b, int ind) { return (double)FSB(b)->surfaceLevels[ind]; });
		if (viscosityModel) {
			scalar("shearRate", [](Block* b, int ind) { return FSB(b)->shearRate[ind]; });
			scalar("viscosity", [](Block* b, int ind) { return FSB(b)->viscosity[ind]; });
		}
	}
	vtk.os.close();

	vtkFiles.push_back(make_pair(iter * blocks[blockId[0]]->Ct, name));
	writePvd("output.pvd", vtkFiles);
	if (writeFields) {
		writeFieldsBinary(iter);
	}
}

// compact output for the GUI: all cells of all blocks, float32, i fastest
void Mesh::writeFieldsBinary(int iter) {
	ofstream os(simulPath + "result/fields_" + std::to_string(iter) + ".bin", ios::out | ios::binary);
	vector<string> names = { "ux", "uy", "uz", "p", "type" };
	if (isFreeSurface()) {
		names.push_back("vof");
		names.push_back("flag");
	}
	Block* b0 = blocks[blockId[0]];
	auto i32 = [&](int32_t v) { os.write((const char*)&v, 4); };
	auto f64 = [&](double v) { os.write((const char*)&v, 8); };
	os.write("LBMF", 4);
	i32(1); // version
	i32(iter);
	f64(iter * b0->Ct);
	i32(problemID);
	i32(nBlock);
	i32(ni);
	i32(nj);
	i32(nk);
	f64(b0->Ch);
	i32((int32_t)names.size());
	for (auto& n : names) {
		char buf[16] = { 0 };
		strncpy(buf, n.c_str(), 15);
		os.write(buf, 16);
	}
	vector<float> data(ni * nj * nk);
	for (int b = 0; b < nBlock; b++) {
		Block* bl = blocks[blockId[b]];
		i32(blockIndex[b][0]);
		i32(blockIndex[b][1]);
		i32(blockIndex[b][2]);
		for (size_t fld = 0; fld < names.size(); fld++) {
			double u[3];
			for (int ind = 0; ind < bl->nAll; ind++) {
				double v = 0;
				switch (fld) {
				case 0:
				case 1:
				case 2:
					if (!bl->boundary[ind]) {
						bl->outputVelocity(ind, u);
						v = bl->Cu * u[fld];
					}
					break;
				case 3:
					v = (bl->boundary[ind] || bl->Rho[ind] == 0) ? 0.0 : bl->Cp * (bl->Rho[ind] - 1);
					break;
				case 4:
					v = bl->cellType[ind];
					break;
				case 5:
					v = FSB(bl)->vof[ind];
					break;
				case 6:
					v = FSB(bl)->flagMap[ind];
					break;
				}
				data[ind] = (float)v;
			}
			os.write((const char*)data.data(), data.size() * sizeof(float));
		}
	}
}

void Mesh::writeBoundaryToVTK() {
	VtkStream vtk(simulPath + "result/boundary.vtk", vtkBinary);
	long long nc;
	writeHexaGrid(vtk, blocks, blockId, blockIndex, nBlock, true, nc);
	vtk.os << "CELL_DATA " << nc << "\nSCALARS cell_type float 1\nLOOKUP_TABLE default\n";
	vector<float> v;
	forCells(blocks, blockId, nBlock, true, [&](Block* b, int ind) { v.push_back((float)b->cellType[ind]); });
	vtk.floats(v, 1);
	vtk.os << "\n";
}

void Mesh::writeFreeSurfaceToVTK(int iter) {
	double isolevel = 0.5;
	GRIDCELL grid;
	TRIANGLE triangles[16];
	vector<TRIANGLE> tri;
	int ei[8] = { 0, 1, 1, 0, 0, 1, 1, 0 };
	int ej[8] = { 0, 0, 1, 1, 0, 0, 1, 1 };
	int ek[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };

	for (int b = 0; b < nBlock; b++) {
		Block* bl = blocks[blockId[b]];
		if (bl->isBlockEmpty) {
			continue;
		}
		int ni = bl->ni;
		int nj = bl->nj;
		int nk = bl->nk;
		for (int k = 0; k < nk; k++) {
			int kg = k + blockIndex[b][2] * nk;
			for (int j = 0; j < nj; j++) {
				int jg = j + blockIndex[b][1] * nj;
				for (int i = 0; i < ni; i++) {
					int ig = i + blockIndex[b][0] * ni;
					int ind = (ni * nj) * k + ni * j + i;
					double vmin = 1e10, vmax = -1e10;
					for (int v = 0; v < 8; v++) {
						grid.p[v].x = ig + ei[v];
						grid.p[v].y = jg + ej[v];
						grid.p[v].z = kg + ek[v];

						// vertex value is the average of the 8 surrounding cells
						double vofAvg = 0.0;
						for (int p = 0; p < 8; p++) {
							int indBlockNeigh, indLocalNeigh;
							bl->local2globalShift(bl->indBlock, ind, indBlockNeigh, indLocalNeigh, ei[v] - 1 + ei[p], ej[v] - 1 + ej[p], ek[v] - 1 + ek[p]);
							if (indBlockNeigh > -1 && indBlockNeigh < nBlockMax && blocks[indBlockNeigh] && !blocks[indBlockNeigh]->boundary[indLocalNeigh]) {
								vofAvg += FSB(blocks[indBlockNeigh])->vof[indLocalNeigh];
							}
							else {
								vofAvg = 8 * (isolevel - 0.01);
								break;
							}
						}
						grid.val[v] = vofAvg / 8;
						vmin = min(vmin, grid.val[v]);
						vmax = max(vmax, grid.val[v]);
					}
					if (vmin < isolevel && vmax > isolevel) {
						int n = Polygonise(grid, isolevel, triangles);
						tri.insert(tri.end(), triangles, triangles + n);
					}
				}
			}
		}
	}

	string name = outputName("outputFS_", iter, ".vtk");
	VtkStream vtk(simulPath + "result/" + name, vtkBinary);
	long long ntri = (long long)tri.size();
	vector<float> pts;
	pts.reserve(9 * ntri);
	for (auto& t : tri) {
		for (int v = 0; v < 3; v++) {
			pts.push_back((float)t.p[v].x);
			pts.push_back((float)t.p[v].y);
			pts.push_back((float)t.p[v].z);
		}
	}
	vtk.os << "DATASET UNSTRUCTURED_GRID\nPOINTS " << ntri * 3 << " float\n";
	vtk.floats(pts, 3);
	vtk.os << "\nCELLS " << ntri << " " << 4 * ntri << "\n";
	vector<int32_t> cells;
	cells.reserve(4 * ntri);
	for (long long t = 0; t < ntri; t++) {
		cells.push_back(3);
		cells.push_back((int32_t)(3 * t));
		cells.push_back((int32_t)(3 * t + 1));
		cells.push_back((int32_t)(3 * t + 2));
	}
	vtk.ints(cells, 4);
	vtk.os << "\nCELL_TYPES " << ntri << "\n";
	vtk.ints(vector<int32_t>(ntri, 5), 1);
	vtk.os.close();

	vtkFsFiles.push_back(make_pair(iter * blocks[blockId[0]]->Ct, name));
	writePvd("outputFS.pvd", vtkFsFiles);
}

// distance from the walls (Dijkstra on the 26-neighbourhood graph, exact fixed point of the former iteration)
void Mesh::generateDistanceFunction() {
	cout << "generating distance function...";
	int nAllB = ni * nj * nk;
	typedef pair<double, long long> Item;
	priority_queue<Item, vector<Item>, greater<Item> > queue;
	vector<int> slot(nBlockMax, -1);
	for (int b = 0; b < nBlock; b++) {
		slot[blockId[b]] = b;
	}

	for (int b = 0; b < nBlock; b++) {
		Block* bl = blocks[blockId[b]];
		for (int ind = 0; ind < nAllB; ind++) {
			bl->wallDistance[ind] = bl->boundary[ind] ? -0.5 : 1e10;
		}
	}
	for (int b = 0; b < nBlock; b++) {
		int indBlock = blockId[b];
		Block* bl = blocks[indBlock];
		for (int ind = 0; ind < nAllB; ind++) {
			if (bl->boundary[ind]) {
				continue;
			}
			double d = 1e10;
			for (int dx = -1; dx < 2; dx++) {
				for (int dy = -1; dy < 2; dy++) {
					for (int dz = -1; dz < 2; dz++) {
						if (dx == 0 && dy == 0 && dz == 0) {
							continue;
						}
						int bs, is;
						bl->local2globalShift(indBlock, ind, bs, is, dx, dy, dz);
						double L = sqrt((double)(dx * dx + dy * dy + dz * dz));
						if (bs > -1 && bs < nBlockMax && blocks[bs]) {
							if (blocks[bs]->boundary[is]) {
								d = min(d, -0.5 + L);
							}
						}
						else {
							d = min(d, L / 2);
						}
					}
				}
			}
			if (d < 1e10) {
				bl->wallDistance[ind] = d;
				queue.push(Item(d, (long long)b * nAllB + ind));
			}
		}
	}
	while (!queue.empty()) {
		Item it = queue.top();
		queue.pop();
		int b = (int)(it.second / nAllB);
		int ind = (int)(it.second % nAllB);
		int indBlock = blockId[b];
		Block* bl = blocks[indBlock];
		if (it.first > bl->wallDistance[ind]) {
			continue;
		}
		for (int dx = -1; dx < 2; dx++) {
			for (int dy = -1; dy < 2; dy++) {
				for (int dz = -1; dz < 2; dz++) {
					if (dx == 0 && dy == 0 && dz == 0) {
						continue;
					}
					int bs, is;
					bl->local2globalShift(indBlock, ind, bs, is, dx, dy, dz);
					if (bs > -1 && bs < nBlockMax && blocks[bs] && !blocks[bs]->boundary[is]) {
						double nd = it.first + sqrt((double)(dx * dx + dy * dy + dz * dz));
						if (nd < blocks[bs]->wallDistance[is]) {
							blocks[bs]->wallDistance[is] = nd;
							queue.push(Item(nd, (long long)slot[bs] * nAllB + is));
						}
					}
				}
			}
		}
	}

	// gradient of distance function
	forBlocks([this](Block* bl) {
		double wallDist[27];
		for (int i = 0; i < bl->nAll; i++) {
			int ind = bl->indAll[i];
			for (int dir = 1; dir < nf; dir++) {
				int bs = bl->nbBlock(ind, dir);
				int is = bl->nbLocal(ind, dir);
				if (bs >= 0 && !blocks[bs]->boundary[is]) {
					wallDist[dir] = blocks[bs]->wallDistance[is];
				}
				else {
					wallDist[dir] = -0.5;
				}
			}
			double* g = &bl->wallDistanceGradient[3 * ind];
			g[0] = g[1] = g[2] = 0;
			for (int dir = 1; dir < nf; dir++) {
				g[0] += 3 * w[dir] * ex[dir] * wallDist[dir];
				g[1] += 3 * w[dir] * ey[dir] * wallDist[dir];
				g[2] += 3 * w[dir] * ez[dir] * wallDist[dir];
			}
		}
	}, false);

	// wall distance polynome
	static const double Ai[64] = {
		1, 0, 0, 0, 0, 0, 0, 0,
		-1, 1, 0, 0, 0, 0, 0, 0,
		-1, 0, 0, 1, 0, 0, 0, 0,
		-1, 0, 0, 0, 1, 0, 0, 0,
		1,-1, 1,-1, 0, 0, 0, 0,
		1, 0, 0,-1,-1, 0, 0, 1,
		1,-1, 0, 0,-1, 1, 0, 0,
		-1, 1,-1, 1, 1,-1, 1,-1 };
	static const int vx[8] = { -1, 1, 1,-1,-1, 1, 1,-1 };
	static const int vy[8] = { -1,-1, 1, 1,-1,-1, 1, 1 };
	static const int vz[8] = { -1,-1,-1,-1, 1, 1, 1, 1 };
	forBlocks([this](Block* bl) {
		double dist[8];
		for (int i = 0; i < bl->nAll; i++) {
			int ind = bl->indAll[i];
			for (int dir = 0; dir < 8; dir++) {
				dist[dir] = 0;
				int bs, is;
				bl->local2globalShift(bl->indBlock, ind, bs, is, vx[dir], vy[dir], vz[dir]);
				if (bs > -1 && bs < nBlockMax && blocks[bs] && !blocks[bs]->boundary[is]) {
					dist[dir] = (bl->wallDistance[ind] + blocks[bs]->wallDistance[is]) / 2;
				}
			}
			for (int p = 0; p < 8; p++) {
				double sum = 0;
				for (int q = 0; q < 8; q++) {
					sum += Ai[8 * p + q] * dist[q];
				}
				bl->wallDistancePolynome[8 * ind + p] = sum;
			}
		}
	}, false);

	cout << " done" << endl;
}

// hydrostatic density profile below the highest fluid cell (in the direction of gravity)
void Mesh::initHydrostaticPress() {
	bool fs = isFreeSurface();
	auto isLiquid = [fs](Block* b, int ind) {
		if (b->boundary[ind]) {
			return false;
		}
		return !fs || FSB(b)->flagMap[ind] == 1 || FSB(b)->flagMap[ind] == 3;
	};
	double* G = blocks[blockId[0]]->G;
	double absG = sqrt(G[0] * G[0] + G[1] * G[1] + G[2] * G[2]);
	if (absG == 0) {
		return;
	}
	double scalarmin = 1e300;
	for (int b = 0; b < nBlock; b++) {
		Block* bl = blocks[blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			if (isLiquid(bl, ind)) {
				int k = ind / (ni * nj), j = (ind % (ni * nj)) / ni, i = ind % ni;
				double scalar = (blockIndex[b][0] * ni + i) * G[0] + (blockIndex[b][1] * nj + j) * G[1] + (blockIndex[b][2] * nk + k) * G[2];
				scalarmin = min(scalarmin, scalar);
			}
		}
	}
	for (int b = 0; b < nBlock; b++) {
		Block* bl = blocks[blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			if (isLiquid(bl, ind)) {
				int k = ind / (ni * nj), j = (ind % (ni * nj)) / ni, i = ind % ni;
				double scalar = (blockIndex[b][0] * ni + i) * G[0] + (blockIndex[b][1] * nj + j) * G[1] + (blockIndex[b][2] * nk + k) * G[2];
				// depth below the top fluid cell is h = (scalar - scalarmin) / |g|, rho = 1 + 3 |g| h
				bl->Rho[ind] += 3.0 * (scalar - scalarmin);
				bl->dRho[ind] = bl->Rho[ind] - bl->rho0;
			}
		}
	}
}

Block* Mesh::cellAt(int gi, int gj, int gk, int& ind) {
	int bi = gi / ni, bj = gj / nj, bk = gk / nk;
	if (gi < 0 || gj < 0 || gk < 0 || bi >= nBi || bj >= nBj || bk >= nBk) {
		return nullptr;
	}
	ind = (ni * nj) * (gk % nk) + ni * (gj % nj) + gi % ni;
	return blocks[(nBi * nBj) * bk + nBi * bj + bi];
}

double Mesh::getCt() {
	return blocks[blockId[0]]->Ct;
}
