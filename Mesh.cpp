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
	curLevel = -1;
	nGhost = 0;
	grid.periodic[0] = Block::getParInt("xPeriodicity", par, false, 0) != 0;
	grid.periodic[1] = Block::getParInt("yPeriodicity", par, false, 0) != 0;
	grid.periodic[2] = Block::getParInt("zPeriodicity", par, false, 0) != 0;
	nBlock = loadStructure(geomPath);
	cout << "Number of blocks in x: " << nBi << endl;
	cout << "Number of blocks in y: " << nBj << endl;
	cout << "Number of blocks in z: " << nBk << endl;
	cout << "Total number of blocks: " << nBlock << endl;
	cout << "Total number of elements: " << (long long)nBlock * ni * nj * nk << endl;
	if (maxLevel > 0) {
		cout << "Refinement levels: " << maxLevel + 1 << " (blocks per level:";
		for (int l = 0; l <= maxLevel; l++) {
			cout << " " << levelBlocks[l].size();
		}
		cout << ", ghost blocks: " << nGhost << ")" << endl;
	}

	// blocks initialization
	blocks = new Block * [nBlockMax];
	for (int i = 0; i < nBlockMax; i++) {
		blocks[i] = nullptr;
	}
	for (int i = 0; i < nBlock; i++) {
		int slot = blockId[i];
		Block* b = createBlock(slot, blockLevel[i], i == 0, par);
		blocks[slot] = b;
		for (int d = 0; d < 3; d++) {
			b->bc[d] = blockIndex[i][d];
		}
		b->grid = &grid;
		b->fileTag = blockLevel[i] == 0 ? to_string(slot) : "L" + to_string(blockLevel[i]) + "_" + to_string(b->bc[0]) + "_" + to_string(b->bc[1]) + "_" + to_string(b->bc[2]);
	}
	for (int l = 0; l <= maxLevel; l++) {
		for (GhostBlock& g : levelGhosts[l]) {
			Block* b = createBlock(g.slot, l, false, par);
			blocks[g.slot] = b;
			for (int d = 0; d < 3; d++) {
				b->bc[d] = g.bc[d];
			}
			b->grid = &grid;
			b->ghost = g.source;
		}
	}
	initGhostMaps();
	initGhostNear();
	explodeStamp.assign(maxLevel + 1, 0);
	for (int l = 0; l <= maxLevel; l++) {
		for (GhostBlock& g : levelGhosts[l]) {
			for (int d = 0; d < 27; d++) {
				int slot = grid.find(l, g.bc[0] + d % 3 - 1, g.bc[1] + d / 3 % 3 - 1, g.bc[2] + d / 9 - 1);
				if (d != 13 && slot >= 0 && slot < nBlockMax && blocks[slot] && blocks[slot]->ghost == GHOST_NONE) {
					g.nbReal.push_back(slot);
				}
			}
		}
	}
	for (int l = 0; l <= maxLevel && maxLevel > 0; l++) {
		for (int slot : levelBlocks[l]) {
			blocks[slot]->massCorr = new double[ni * nj * nk]();
		}
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
	for (auto& level : levelGhosts) {
		for (GhostBlock& g : level) {
			delete blocks[g.slot];
		}
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
	if (maxLevel > 0) {
		syncGhostGeometry();
		for (int l = 0; l <= maxLevel; l++) {
			forGhosts(l, 0, [this](GhostBlock& g) {
				blocks[g.slot]->initZeroCondition();
				blocks[g.slot]->computeShiftIndexes();
			}, true);
		}
		syncGhosts();
	}

	if (isFreeSurface()) {
		for (int i = 0; i < nBlock; i++) {
			FSB(blocks[blockId[i]])->initFlagMap();
		}
		syncGhosts();
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
	syncGhosts();

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
		forBlocks([](Block* b) { FSB(b)->calculMassTotal(); }, false);
		massTot = 0;
		activeCellsNumTot = 0;
		for (int i = 0; i < nBlock; i++) {
			Block* b = blocks[blockId[i]];
			massTot += FSB(b)->massBlock;
			if (!b->isBlockEmpty) {
				activeCellsNumTot += (long long)FSB(b)->activeCellsNumBlock << b->level; // 2^level steps per iteration
			}
		}
		break;
	}

	case 11:
		forBlocks([](Block* b) { b->isBlockTouched = false; }, false);
		forBlocks([](Block* b) {
			FSB(b)->checkIfBlockEmpty();
			FSB(b)->touchNeighBlocks();
		});
		forBlocks([](Block* b) {
			if (b->isBlockTouched) {
				b->isBlockEmpty = false;
			}
		}, false);
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
		double t = getCt() * iter;
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

// one step of the level and (recursively) two steps of the finer levels, iterLevel counts the steps of the level
void Mesh::advanceLevel(int level, int iterLevel) {
	if (level < maxLevel) {
		computeLost(level);
		explodeToGhosts(level + 1);
		advanceLevel(level + 1, 2 * iterLevel);
		advanceLevel(level + 1, 2 * iterLevel + 1);
		coalesceToGhosts(level, true);
	}
	if (level > 0) {
		refreshGhosts(level);
	}
	curLevel = level;
	switch (problemID) {
	case 0:
		evalBlocks(iterLevel, 0);
		if (level > 0) {
			forGhosts(level, GHOST_FROM_COARSER, [this](GhostBlock& g) { blocks[g.slot]->propagateCells(g.near2); });
			accountInterface(level);
		}
		if (level < maxLevel) {
			applyMassCorr(level);
		}
		evalBlocks(iterLevel, 1);
		if (level > 0) {
			// the coarser populations stream through the ghost layer during the two fine steps without collision
			forGhosts(level, GHOST_FROM_COARSER, [this](GhostBlock& g) {
				Block* b = blocks[g.slot];
				for (int ind : g.near2) {
					for (int p = 0; p < nf; p++) {
						b->fc[nf * ind + p] = b->f[nf * ind + p];
					}
				}
			});
		}
		break;
	case 1:
		advanceFreeSurfaceLevel(level, iterLevel);
		break;
	}
	curLevel = -1;
}

// one free surface step of the level (Mesh::timeStep, case 1), the ghosts take part in the local operations
void Mesh::advanceFreeSurfaceLevel(int level, int iterLevel) {
	bool coarser = level > 0, finer = level < maxLevel;
	auto ghostsDo = [this, level](auto fn) { forGhosts(level, 0, [this, &fn](GhostBlock& g) { fn(FSB(blocks[g.slot])); }); };

	// mass exchange, find filled/emptied cells
	forBlocks([](Block* b) { FSB(b)->applyMassAcc(); });
	forBlocks([](Block* b) { FSB(b)->computeInterfaceTypes(); });
	ghostsDo([](BlockFreeSurface* b) { b->computeInterfaceTypes(); });
	forBlocks([](Block* b) { FSB(b)->massExchange(); });
	if (coarser) {
		accountFsExchange(level);
	}
	forBlocks([iterLevel](Block* b) { FSB(b)->findFilledEmptied(iterLevel); });
	ghostsDo([](BlockFreeSurface* b) { b->findSurfaceLevels(1); });
	forBlocks([](Block* b) { FSB(b)->findSurfaceLevels(1); });
	for (int j = 2; j <= 4; j++) {
		forBlocks([j](Block* b) { FSB(b)->findSurfaceLevels(j); });
		ghostsDo([j](BlockFreeSurface* b) { b->findSurfaceLevels(j); });
	}
	evalBlocks(iterLevel, 3);
	evalBlocks(iterLevel, 4);
	forBlocks([](Block* b) { FSB(b)->injectPending(); });
	for (int s = 0; s < nSmoothStep + 1; s++) {
		int ns = nSmoothStep;
		forBlocks([s, ns](Block* b) { FSB(b)->smoothVOF(s, ns); });
		ghostsDo([s, ns](BlockFreeSurface* b) { b->smoothVOF(s, ns); });
	}
	forBlocks([](Block* b) {
		FSB(b)->calculateNormalVectors();
		FSB(b)->setWetAngle();
	});
	ghostsDo([](BlockFreeSurface* b) {
		b->calculateNormalVectors();
		b->setWetAngle();
	});
	forBlocks([](Block* b) { FSB(b)->calculateCurvature(); });
	ghostsDo([](BlockFreeSurface* b) { b->calculateCurvature(); });

	// propagation, the populations of the coarser level stream through the ghost layer
	evalBlocks(iterLevel, 0);
	if (coarser) {
		forGhosts(level, GHOST_FROM_COARSER, [this](GhostBlock& g) { blocks[g.slot]->propagateCells(g.near2); });
		accountInterface(level);
	}
	applyMassCorr(level);

	// macro, collision, new free surface geometry (also of the ghost cells, the requests go to their level)
	forBlocks([](Block* b) { b->makro(); });
	forBlocks([](Block* b) {
		b->collision();
		FSB(b)->newFsGeometryPrepare();
	});
	ghostsDo([](BlockFreeSurface* b) { b->newFsGeometryPrepare(); });
	forBlocks([](Block* b) { FSB(b)->newFsGeometryApply(); });
	ghostsDo([](BlockFreeSurface* b) { b->newFsGeometryApply(); });
	forwardGhostRequests(level);
	if (coarser) {
		forGhosts(level, GHOST_FROM_COARSER, [this](GhostBlock& g) {
			BlockFreeSurface* b = FSB(blocks[g.slot]);
			for (int ind : g.near2) {
				if (b->flagMap[ind] == 1 || b->flagMap[ind] == 3) {
					for (int p = 0; p < nf; p++) {
						b->fc[nf * ind + p] = b->f[nf * ind + p];
					}
				}
			}
		});
	}

	evalBlocks(iterLevel, 9);  // destroy artefacts
	evalBlocks(iterLevel, 10); // mass redistribution
	forwardGhostMass(level);
	evalBlocks(iterLevel, 11); // check empty block flag
	activateFromGhosts(level);
}

// one time step of the selected problem
void Mesh::timeStep(int iter) {
	if (maxLevel > 0) {
		advanceLevel(0, iter);
		return;
	}
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

Block* Mesh::createBlock(int indBlock, int level, bool isFirst, const ParamMap& par) {
	switch (problemID) {
	case 0:
		return new BlockBasic(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel != 0, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
	case 1:
		return new BlockFreeSurface(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
	case 2:
		return new BlockFreeSurfaceLagrangeParticles(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
	case 21:
		return new BlockFreeSurfaceLagrangeFibers(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
	case 3:
		return new BlockFreeSurfaceEulerParticlesFVM(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
	case 4:
		return new BlockFreeSurfaceEulerParticlesLBM(isFirst, blocks, nBlock, nBlockMax, blockMap, nBi, nBj, nBk, indBlock, level, indAll, indIn, indBound, collisionModel, viscosityModel, nf, ex, ey, ez, w, opp, ni, nj, nk, par);
	default:
		throw invalid_argument("unknown problemID " + to_string(problemID));
	}
}

// blocks.txt: "bi bj bk" (level 0) or "bi bj bk level" with the block coordinates at its level
int Mesh::loadStructure(string geomPath) {
	string fileName = geomPath + "mesh/blocks.txt";
	ifstream myfile(fileName);
	if (!myfile.is_open()) {
		throw runtime_error("cannot open " + fileName);
	}
	vector<array<int, 4> > list;
	string line;
	nBi = nBj = nBk = 0;
	maxLevel = 0;
	while (getline(myfile, line)) {
		istringstream iss(line);
		array<int, 4> idx = { 0, 0, 0, 0 };
		if (!(iss >> idx[0] >> idx[1] >> idx[2])) {
			continue; // empty line
		}
		if (!(iss >> idx[3])) {
			idx[3] = 0;
		}
		if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0 || idx[3] < 0 || idx[3] > 15) {
			throw runtime_error("invalid block index in " + fileName + ": '" + line + "'");
		}
		maxLevel = max(maxLevel, idx[3]);
		nBi = max(nBi, (idx[0] >> idx[3]) + 1);
		nBj = max(nBj, (idx[1] >> idx[3]) + 1);
		nBk = max(nBk, (idx[2] >> idx[3]) + 1);
		list.push_back(idx);
	}
	if (list.empty()) {
		throw runtime_error("no blocks in " + fileName);
	}
	if (maxLevel > 0) {
		if (nf == 9) {
			throw invalid_argument("grid refinement is not available for the D2Q9 model");
		}
		if (ni % 2 || nj % 2 || nk % 2) {
			throw invalid_argument("grid refinement requires even block dimensions ni, nj, nk");
		}
		if (problemID != 0 && problemID != 1) {
			throw invalid_argument("grid refinement is available for problemID 0 and 1 only");
		}
	}
	if (nk == 1 && nBk == 1) {
		grid.periodic[2] = true; // 2D model, see Block::initCommon
	}
	grid.nBi = nBi;
	grid.nBj = nBj;
	grid.nBk = nBk;
	grid.maxLevel = maxLevel;
	int nLinear = nBi * nBj * nBk;
	int nBlocks = (int)list.size();
	blockMap = new int[nLinear];
	for (int i = 0; i < nLinear; i++) {
		blockMap[i] = -1;
	}
	blockIndex = new int* [nBlocks];
	blockId = new int[nBlocks];
	blockLevel.assign(nBlocks, 0);
	levelBlocks.assign(maxLevel + 1, vector<int>());
	int nextSlot = nLinear;
	for (int s = 0; s < nBlocks; s++) {
		int l = list[s][3];
		blockIndex[s] = new int[3]{ list[s][0], list[s][1], list[s][2] };
		blockLevel[s] = l;
		int slot;
		if (l == 0) {
			slot = (nBi * nBj) * list[s][2] + nBi * list[s][1] + list[s][0];
			if (blockMap[slot] != -1) {
				throw runtime_error("duplicate block in " + fileName);
			}
			blockMap[slot] = s;
		}
		else {
			unsigned long long key = BlockGrid::key(l, list[s][0], list[s][1], list[s][2]);
			if (grid.slots.count(key)) {
				throw runtime_error("duplicate block in " + fileName);
			}
			slot = nextSlot++;
			grid.slots[key] = slot;
		}
		blockId[s] = slot;
		levelBlocks[l].push_back(slot);
	}
	nBlockMax = nextSlot;
	createGhosts(list);
	return nBlocks;
}

// ghost blocks around the blocks of every level where the neighbourhood is held by another level
void Mesh::createGhosts(const vector<array<int, 4> >& list) {
	levelGhosts.assign(maxLevel + 1, vector<GhostBlock>());
	if (maxLevel == 0) {
		return;
	}
	auto exists = [this](int l, int i, int j, int k) {
		if (l < 0 || l > maxLevel) {
			return false;
		}
		int slot = grid.find(l, i, j, k);
		if (slot < 0) {
			return false;
		}
		return l == 0 ? blockMap[slot] >= 0 : true;
	};
	for (auto& b : list) {
		for (int l = 0; l < b[3]; l++) {
			int sh = b[3] - l;
			if (exists(l, b[0] >> sh, b[1] >> sh, b[2] >> sh)) {
				throw runtime_error("blocks.txt: block (" + to_string(b[0]) + ", " + to_string(b[1]) + ", " + to_string(b[2]) + ") of the level " + to_string(b[3]) + " overlaps a block of the level " + to_string(l));
			}
		}
	}
	unordered_map<unsigned long long, int> known;
	int nextSlot = nBlockMax;
	for (auto& b : list) {
		int l = b[3];
		for (int dk = -1; dk <= 1; dk++) {
			for (int dj = -1; dj <= 1; dj++) {
				for (int di = -1; di <= 1; di++) {
					int p[3] = { b[0] + di, b[1] + dj, b[2] + dk };
					if ((di == 0 && dj == 0 && dk == 0) || !grid.wrap(l, p[0], p[1], p[2]) || exists(l, p[0], p[1], p[2])) {
						continue;
					}
					unsigned long long key = BlockGrid::key(l, p[0], p[1], p[2]);
					if (known.count(key)) {
						continue;
					}
					int source = GHOST_NONE;
					if (exists(l - 1, p[0] >> 1, p[1] >> 1, p[2] >> 1)) {
						source = GHOST_FROM_COARSER;
					}
					else {
						for (int o = 0; o < 8; o++) {
							int c[3] = { 2 * p[0] + (o & 1), 2 * p[1] + ((o >> 1) & 1), 2 * p[2] + (o >> 2) };
							if (exists(l + 1, c[0], c[1], c[2])) {
								source = GHOST_FROM_FINER;
							}
							for (int g = 0; g < 8; g++) {
								if (exists(l + 2, 2 * c[0] + (g & 1), 2 * c[1] + ((g >> 1) & 1), 2 * c[2] + (g >> 2))) {
									// the half of the neighbour block next to this block has to be of the level l + 1
									bool near = true;
									int off[3] = { di, dj, dk };
									int co[3] = { o & 1, (o >> 1) & 1, o >> 2 };
									for (int d = 0; d < 3; d++) {
										near = near && (off[d] == 0 || co[d] == (off[d] > 0 ? 0 : 1));
									}
									if (near) {
										throw runtime_error("blocks.txt: the mesh is not 2:1 balanced, a block of the level " + to_string(l) + " touches a block of the level " + to_string(l + 2));
									}
								}
							}
						}
						if (source == GHOST_NONE) {
							for (int lc = l - 2; lc >= 0; lc--) {
								int sh = l - lc;
								if (exists(lc, p[0] >> sh, p[1] >> sh, p[2] >> sh)) {
									throw runtime_error("blocks.txt: the mesh is not 2:1 balanced, a block of the level " + to_string(l) + " touches a block of the level " + to_string(lc));
								}
							}
						}
					}
					if (source == GHOST_NONE) {
						continue; // exterior of the domain
					}
					known[key] = 1;
					GhostBlock gb;
					gb.slot = (l == 0) ? grid.find(0, p[0], p[1], p[2]) : -1; // slots of the finer levels are assigned below
					gb.source = source;
					for (int d = 0; d < 3; d++) {
						gb.bc[d] = p[d];
					}
					levelGhosts[l].push_back(gb);
					nGhost++;
				}
			}
		}
	}
	// the hash holds only the computed blocks during the search above
	for (int l = 1; l <= maxLevel; l++) {
		for (GhostBlock& g : levelGhosts[l]) {
			g.slot = nextSlot++;
			grid.slots[BlockGrid::key(l, g.bc[0], g.bc[1], g.bc[2])] = g.slot;
		}
	}
	nBlockMax = nextSlot;
}

// source cells of the ghost blocks
void Mesh::initGhostMaps() {
	int nAllB = ni * nj * nk;
	for (int l = 0; l <= maxLevel; l++) {
		for (GhostBlock& g : levelGhosts[l]) {
			if (g.source == GHOST_FROM_COARSER) {
				int slot = grid.find(l - 1, g.bc[0] >> 1, g.bc[1] >> 1, g.bc[2] >> 1);
				g.srcBlock.assign(nAllB, slot);
				g.srcLocal.assign(nAllB, 0);
				for (int k = 0; k < nk; k++) {
					for (int j = 0; j < nj; j++) {
						for (int i = 0; i < ni; i++) {
							int ci = ((g.bc[0] & 1) * ni + i) >> 1;
							int cj = ((g.bc[1] & 1) * nj + j) >> 1;
							int ck = ((g.bc[2] & 1) * nk + k) >> 1;
							g.srcLocal[(ni * nj) * k + ni * j + i] = (ni * nj) * ck + ni * cj + ci;
						}
					}
				}
			}
			else {
				g.srcBlock.assign(8 * nAllB, -1);
				g.srcLocal.assign(8 * nAllB, 0);
				for (int k = 0; k < nk; k++) {
					for (int j = 0; j < nj; j++) {
						for (int i = 0; i < ni; i++) {
							int ind = (ni * nj) * k + ni * j + i;
							for (int o = 0; o < 8; o++) {
								int fi = 2 * i + (o & 1), fj = 2 * j + ((o >> 1) & 1), fk = 2 * k + (o >> 2);
								int slot = grid.find(l + 1, 2 * g.bc[0] + fi / ni, 2 * g.bc[1] + fj / nj, 2 * g.bc[2] + fk / nk);
								if (slot >= 0 && blocks[slot] && blocks[slot]->ghost == GHOST_NONE) {
									g.srcBlock[8 * ind + o] = slot;
									g.srcLocal[8 * ind + o] = (ni * nj) * (fk % nk) + ni * (fj % nj) + fi % ni;
								}
							}
						}
					}
				}
			}
		}
	}
}

// distance (in cells of the level, 26-neighbourhood) of the ghost cells from the computed cells, 1 or 2
void Mesh::initGhostNear() {
	vector<int> ghostIndex(nBlockMax, -1);
	vector<vector<char> > dist(0);
	vector<GhostBlock*> all;
	for (int l = 0; l <= maxLevel; l++) {
		for (GhostBlock& g : levelGhosts[l]) {
			ghostIndex[g.slot] = (int)all.size();
			all.push_back(&g);
		}
	}
	// free surface stencils of the computed cells reach 3 cells (height function), the smoothed volume fraction,
	// normals and curvature of the ghost cells they read need 3 more
	const int GHOST_HALO = 6;
	const char FAR = 100;
	dist.assign(all.size(), vector<char>(ni * nj * nk, FAR));
	for (int pass = 1; pass <= GHOST_HALO; pass++) {
#pragma omp parallel for schedule(dynamic)
		for (int gi = 0; gi < (int)all.size(); gi++) {
			Block* b = blocks[all[gi]->slot];
			for (int ind = 0; ind < b->nAll; ind++) {
				if (dist[gi][ind] < FAR) {
					continue;
				}
				bool hit = false;
				for (int dz = -1; dz <= 1 && !hit; dz++) {
					for (int dy = -1; dy <= 1 && !hit; dy++) {
						for (int dx = -1; dx <= 1 && !hit; dx++) {
							int bl, loc;
							b->local2globalShift(b->indBlock, ind, bl, loc, dx, dy, dz);
							if (bl < 0 || !blocks[bl]) {
								continue;
							}
							if (pass == 1) {
								hit = blocks[bl]->ghost == GHOST_NONE;
							}
							else {
								hit = ghostIndex[bl] >= 0 && dist[ghostIndex[bl]][loc] == pass - 1;
							}
						}
					}
				}
				if (hit) {
					dist[gi][ind] = (char)(-pass); // marked, finalized below so that the pass reads only older marks
				}
			}
		}
		for (auto& d : dist) {
			for (char& c : d) {
				if (c < 0) {
					c = (char)-c;
				}
			}
		}
	}
	for (size_t gi = 0; gi < all.size(); gi++) {
		all[gi]->near1.clear();
		all[gi]->near2.clear();
		all[gi]->halo.clear();
		for (int ind = 0; ind < ni * nj * nk; ind++) {
			if (dist[gi][ind] <= GHOST_HALO) {
				all[gi]->halo.push_back(ind);
			}
			if (dist[gi][ind] == 1) {
				all[gi]->near1.push_back(ind);
			}
			if (dist[gi][ind] <= 2) {
				all[gi]->near2.push_back(ind);
			}
		}
	}
}

// cell types of the ghosts: copy of the coarser cell or the prevailing type of the 8 finer cells
void Mesh::syncGhostGeometry() {
	for (int l = 0; l <= maxLevel; l++) {
		forGhosts(l, 0, [this](GhostBlock& g) {
			Block* b = blocks[g.slot];
			for (int ind = 0; ind < b->nAll; ind++) {
				int type = CT_WALL;
				if (g.source == GHOST_FROM_COARSER) {
					type = blocks[g.srcBlock[ind]]->cellType[g.srcLocal[ind]];
				}
				else {
					int count[CT_MOVING + 1] = { 0 };
					for (int o = 0; o < 8; o++) {
						int sb = g.srcBlock[8 * ind + o];
						count[sb < 0 ? CT_WALL : blocks[sb]->cellType[g.srcLocal[8 * ind + o]]]++;
					}
					if (count[CT_FLUID] >= 4) {
						type = CT_FLUID;
					}
					else {
						for (int t = CT_WALL + 1; t <= CT_MOVING; t++) {
							if (count[t] > count[type]) {
								type = t;
							}
						}
					}
				}
				b->cellType[ind] = (unsigned char)type;
				b->boundary[ind] = type != CT_FLUID;
			}
		}, true);
	}
}

// state of all ghosts from the current state of the computed blocks (initialization)
void Mesh::syncGhosts() {
	for (int l = 1; l <= maxLevel; l++) {
		explodeToGhosts(l, true);
	}
	for (int l = maxLevel - 1; l >= 0; l--) {
		coalesceToGhosts(l, false, true);
	}
}

void Mesh::explodeGhost(GhostBlock& g) {
	Block* b = blocks[g.slot];
	b->resetGhost();
	for (int ind : g.halo) {
		b->prolongCell(ind, blocks[g.srcBlock[ind]], g.srcLocal[ind]);
	}
	g.stamp = explodeStamp[b->level];
}

void Mesh::explodeToGhosts(int level, bool all) {
	explodeStamp[level]++;
	forGhosts(level, GHOST_FROM_COARSER, [this](GhostBlock& g) { explodeGhost(g); }, all);
}

void Mesh::refreshGhosts(int level) {
	int stamp = explodeStamp[level];
	forGhosts(level, GHOST_FROM_COARSER, [this, stamp](GhostBlock& g) {
		if (g.stamp != stamp) {
			explodeGhost(g); // the coarser level does not change during the finer steps
		}
	});
}

void Mesh::coalesceToGhosts(int level, bool streamed, bool all) {
	forGhosts(level, GHOST_FROM_FINER, [this, level, streamed](GhostBlock& g) {
		Block* b = blocks[g.slot];
		Block* src[8];
		b->resetGhost();
		if (!streamed || isFreeSurface()) {
			for (int ind : g.halo) {
				for (int o = 0; o < 8; o++) {
					int sb = g.srcBlock[8 * ind + o];
					src[o] = sb < 0 ? nullptr : blocks[sb];
				}
				b->restrictCell(ind, src, &g.srcLocal[8 * ind]);
			}
		}
		if (!streamed) {
			return;
		}
		// populations entering the computed cells of this level: average of the populations that streamed
		// into the finer ghost cells covering the destination cell during the two fine steps (Rohde et al. 2006)
		for (int ind : g.near1) {
			if (b->massExchangeKind(ind) == 0) {
				continue; // wall or gas (the destination reconstructs the population)
			}
			for (int dir = 1; dir < nf; dir++) {
				int bl = b->nbBlock(ind, dir);
				int loc = b->nbLocal(ind, dir);
				if (bl < 0 || blocks[bl]->ghost != GHOST_NONE || blocks[bl]->massExchangeKind(loc) == 0) {
					continue;
				}
				Block* dest = blocks[bl];
				int gi, gj, gk;
				dest->globalCoords(loc, gi, gj, gk);
				double sum = 0;
				int count = 0;
				for (int o = 0; o < 8; o++) {
					int fl;
					Block* fb = cellAtLevel(level + 1, 2 * gi + (o & 1), 2 * gj + ((o >> 1) & 1), 2 * gk + (o >> 2), fl);
					if (fb && fb->massExchangeKind(fl) != 0) {
						sum += fb->f[nf * fl + dir];
						count++;
					}
				}
				if (count == 0) {
					continue;
				}
				double v = sum / count;
				b->fc[nf * ind + dir] = (pop_t)v;
				if (dest->massExchangeKind(loc) == 1) {
					double* mc = &dest->massCorr[loc];
#pragma omp atomic
					*mc -= v + w[dir] * dest->rho0; // received by the coarser level
				}
			}
		}
	}, all);
}

// mass of the post-collision populations of the level streaming towards the finer level
void Mesh::computeLost(int level) {
	int saved = curLevel;
	curLevel = level;
	forBlocks([this](Block* b) {
		for (int i = 0; i < b->nBound; i++) {
			int ind = b->indBound[i];
			if (b->massExchangeKind(ind) != 1) {
				continue;
			}
			for (int p = 1; p < nf; p++) {
				int bl = b->nbBlock(ind, p);
				int loc = b->nbLocal(ind, p);
				if (bl >= 0 && blocks[bl]->ghost == GHOST_FROM_FINER && !blocks[bl]->boundary[loc]) {
					b->massCorr[ind] += b->fc[nf * ind + p] + w[p] * b->rho0;
				}
			}
		}
	}, false);
	curLevel = saved;
}

// populations exchanged by the computed cells of the level and its ghosts from the coarser level in the last
// propagation, the balance goes to the coarser cell covering the ghost cell
void Mesh::accountInterface(int level) {
	forGhosts(level, GHOST_FROM_COARSER, [this](GhostBlock& g) {
		Block* b = blocks[g.slot];
		for (int ind : g.near1) {
			if (b->massExchangeKind(ind) != 1) {
				continue; // exchange with free surface cells is balanced in accountFsExchange
			}
			double d = 0;
			for (int i = 1; i < nf; i++) {
				int bl = b->nbBlock(ind, i);
				int loc = b->nbLocal(ind, i);
				if (bl < 0 || blocks[bl]->ghost != GHOST_NONE || blocks[bl]->massExchangeKind(loc) == 0) {
					continue;
				}
				d -= b->fc[nf * ind + i] + w[i] * b->rho0;       // received by the computed cell
				d += b->f[nf * ind + opp[i]] + w[i] * b->rho0;   // sent by the computed cell
			}
			if (d != 0) {
				blocks[g.srcBlock[ind]]->massCorr[g.srcLocal[ind]] += d / 8;
			}
		}
	});
}

// free surface: the mass exchanged by the computed cells with the interface ghost cells from the coarser level
// is taken from the coarser cell (the finer cells are authoritative, the coarser ones skip these links)
void Mesh::accountFsExchange(int level) {
	forGhosts(level, GHOST_FROM_COARSER, [this](GhostBlock& g) {
		BlockFreeSurface* b = FSB(blocks[g.slot]);
		for (int ind : g.near1) {
			if (b->massExchangeKind(ind) != 2) {
				continue;
			}
			double d = 0;
			for (int i = 1; i < nf; i++) {
				int bl = b->nbBlock(ind, i);
				int loc = b->nbLocal(ind, i);
				if (bl < 0 || blocks[bl]->ghost != GHOST_NONE) {
					continue;
				}
				BlockFreeSurface* c = FSB(blocks[bl]);
				int kind = c->massExchangeKind(loc);
				if (kind == 2) {
					d -= c->linkExchange(loc, opp[i], b, ind);
				}
				else if (kind == 1) {
					d -= b->fc[nf * ind + i] - c->fc[nf * loc + opp[i]];
				}
			}
			if (d != 0) {
				FSB(blocks[g.srcBlock[ind]])->massAcc[g.srcLocal[ind]] += d / 8;
			}
		}
	});
}

// conversions of ghost cells requested by the computed cells of the level are applied by the level holding them
void Mesh::forwardGhostRequests(int level) {
	for (GhostBlock& g : levelGhosts[level]) {
		BlockFreeSurface* b = FSB(blocks[g.slot]);
		for (int c = 0; c < b->nAirToFs + b->nFluidToFs; c++) {
			bool air = c < b->nAirToFs;
			int ind = air ? b->airToFs[c] : b->fluidToFs[c - b->nAirToFs];
			if (g.source == GHOST_FROM_COARSER) {
				BlockFreeSurface* t = FSB(blocks[g.srcBlock[ind]]);
				(air ? t->pendingAir : t->pendingFluid).push_back(g.srcLocal[ind]);
			}
			else {
				for (int o = 0; o < 8; o++) {
					int sb = g.srcBlock[8 * ind + o];
					if (sb >= 0) {
						BlockFreeSurface* t = FSB(blocks[sb]);
						(air ? t->pendingAir : t->pendingFluid).push_back(g.srcLocal[8 * ind + o]);
					}
				}
			}
		}
		b->nAirToFs = b->nFluidToFs = 0;
	}
}

// mass redistributed into ghost cells goes to the cells of the level holding them
void Mesh::forwardGhostMass(int level) {
	forGhosts(level, 0, [this](GhostBlock& g) {
		BlockFreeSurface* b = FSB(blocks[g.slot]);
		for (int ind = 0; ind < b->nAll; ind++) {
			double dm = b->massAcc[ind];
			if (dm == 0) {
				continue;
			}
			b->massAcc[ind] = 0;
			if (g.source == GHOST_FROM_COARSER) {
				FSB(blocks[g.srcBlock[ind]])->massAcc[g.srcLocal[ind]] += dm / 8;
				continue;
			}
			// 8 finer cells hold the mass of the coarser cell, interface cells take it first
			int target[8];
			int n = 0;
			for (int kind = 2; kind >= 1 && n == 0; kind--) {
				for (int o = 0; o < 8; o++) {
					int sb = g.srcBlock[8 * ind + o];
					if (sb >= 0 && blocks[sb]->massExchangeKind(g.srcLocal[8 * ind + o]) == kind) {
						target[n++] = o;
					}
				}
			}
			if (n == 0) {
				for (int o = 0; o < 8 && n == 0; o++) {
					if (g.srcBlock[8 * ind + o] >= 0) {
						target[n++] = o;
					}
				}
			}
			for (int t = 0; t < n; t++) {
				int o = target[t];
				FSB(blocks[g.srcBlock[8 * ind + o]])->massAcc[g.srcLocal[8 * ind + o]] += 8 * dm / n;
			}
		}
	});
}

// blocks touched through a ghost block by the surface of another level become active
void Mesh::activateFromGhosts(int level) {
	for (GhostBlock& g : levelGhosts[level]) {
		Block* b = blocks[g.slot];
		if (!b->isBlockTouched) {
			continue;
		}
		b->isBlockTouched = false;
		if (g.source == GHOST_FROM_COARSER) {
			blocks[g.srcBlock[0]]->isBlockEmpty = false;
		}
		else {
			for (size_t c = 0; c < g.srcBlock.size(); c++) {
				if (g.srcBlock[c] >= 0) {
					blocks[g.srcBlock[c]]->isBlockEmpty = false;
				}
			}
		}
	}
}

// corrections go only to fluid cells, empty blocks have none (and are skipped)
void Mesh::applyMassCorr(int level) {
	int saved = curLevel;
	curLevel = level;
	forBlocks([this](Block* b) {
		for (int ind = 0; ind < b->nAll; ind++) {
			if (b->massCorr[ind] != 0 && b->massExchangeKind(ind) != 1) {
				b->massCorr[ind] = 0; // the mass of the free surface cells is balanced separately
			}
			if (b->massCorr[ind] != 0) {
				// distributed like the density part of the equilibrium, the momentum stays
				double dm = b->massCorr[ind];
				for (int p = 0; p < nf; p++) {
					b->f[nf * ind + p] += (pop_t)(w[p] * dm);
				}
				b->massCorr[ind] = 0;
			}
		}
	});
	curLevel = saved;
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
	long long us = llround(iter * getCt() * 1000000);
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
	// coordinates are in the cells of the level 0
	void writeHexaGrid(VtkStream& vtk, Block** blocks, int* blockId, int nBlock, bool selectBoundary, long long& nc) {
		vector<float> pts;
		for (int b = 0; b < nBlock; b++) {
			Block* bl = blocks[blockId[b]];
			double h = 1.0 / (1 << bl->level);
			for (int k = 0; k < bl->nk + 1; k++) {
				for (int j = 0; j < bl->nj + 1; j++) {
					for (int i = 0; i < bl->ni + 1; i++) {
						pts.push_back((float)(h * (i + bl->ni * bl->bc[0])));
						pts.push_back((float)(h * (j + bl->nj * bl->bc[1])));
						pts.push_back((float)(h * (k + bl->nk * bl->bc[2])));
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
	writeHexaGrid(vtk, blocks, blockId, nBlock, false, nc);
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

	vtkFiles.push_back(make_pair(iter * getCt(), name));
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
	i32(maxLevel > 0 ? 2 : 1); // version 2 adds the level of the blocks
	i32(iter);
	f64(iter * b0->Ct);
	i32(problemID);
	i32(nBlock);
	i32(ni);
	i32(nj);
	i32(nk);
	f64(b0->Ch * (1 << b0->level)); // cell size of the level 0
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
		if (maxLevel > 0) {
			i32(blockLevel[b]);
		}
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
	writeHexaGrid(vtk, blocks, blockId, nBlock, true, nc);
	vtk.os << "CELL_DATA " << nc << "\nSCALARS cell_type float 1\nLOOKUP_TABLE default\n";
	vector<float> v;
	forCells(blocks, blockId, nBlock, true, [&](Block* b, int ind) { v.push_back((float)b->cellType[ind]); });
	vtk.floats(v, 1);
	vtk.os << "\n";
}

void Mesh::writeFreeSurfaceToVTK(int iter) {
	double isolevel = 0.5;
	GRIDCELL cube;
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
		double h = 1.0 / (1 << bl->level);
		for (int k = 0; k < nk; k++) {
			int kg = k + bl->bc[2] * nk;
			for (int j = 0; j < nj; j++) {
				int jg = j + bl->bc[1] * nj;
				for (int i = 0; i < ni; i++) {
					int ig = i + bl->bc[0] * ni;
					int ind = (ni * nj) * k + ni * j + i;
					double vmin = 1e10, vmax = -1e10;
					for (int v = 0; v < 8; v++) {
						cube.p[v].x = h * (ig + ei[v]);
						cube.p[v].y = h * (jg + ej[v]);
						cube.p[v].z = h * (kg + ek[v]);

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
						cube.val[v] = vofAvg / 8;
						vmin = min(vmin, cube.val[v]);
						vmax = max(vmax, cube.val[v]);
					}
					if (vmin < isolevel && vmax > isolevel) {
						int n = Polygonise(cube, isolevel, triangles);
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

	vtkFsFiles.push_back(make_pair(iter * getCt(), name));
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
							if (blocks[bs]->ghost != GHOST_NONE) {
								continue; // other level
							}
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
					if (bs > -1 && bs < nBlockMax && blocks[bs] && blocks[bs]->ghost == GHOST_NONE && !blocks[bs]->boundary[is]) {
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
	// rho = 1 + 3 (x . g) / Cu^2 relative to the highest fluid cell, x in metres
	auto potential = [](Block* b, int ind) {
		int gi, gj, gk;
		b->globalCoords(ind, gi, gj, gk);
		return ((gi + 0.5) * b->G[0] + (gj + 0.5) * b->G[1] + (gk + 0.5) * b->G[2]) * b->Ch * b->Cf;
	};
	Block* r = blocks[blockId[0]];
	if (r->G[0] == 0 && r->G[1] == 0 && r->G[2] == 0) {
		return;
	}
	double scalarmin = 1e300;
	for (int b = 0; b < nBlock; b++) {
		Block* bl = blocks[blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			if (isLiquid(bl, ind)) {
				scalarmin = min(scalarmin, potential(bl, ind));
			}
		}
	}
	for (int b = 0; b < nBlock; b++) {
		Block* bl = blocks[blockId[b]];
		for (int ind = 0; ind < bl->nAll; ind++) {
			if (isLiquid(bl, ind)) {
				bl->Rho[ind] += 3.0 * (potential(bl, ind) - scalarmin) / (bl->Cu * bl->Cu);
				bl->dRho[ind] = bl->Rho[ind] - bl->rho0;
			}
		}
	}
}

Block* Mesh::cellAt(int gi, int gj, int gk, int& ind) {
	Block* b = cellAtLevel(0, gi, gj, gk, ind);
	return (b && b->ghost == GHOST_NONE) ? b : nullptr;
}

Block* Mesh::cellAtLevel(int level, int gi, int gj, int gk, int& ind) {
	if (gi < 0 || gj < 0 || gk < 0) {
		return nullptr;
	}
	int slot = grid.find(level, gi / ni, gj / nj, gk / nk);
	if (slot < 0 || slot >= nBlockMax) {
		return nullptr;
	}
	ind = (ni * nj) * (gk % nk) + ni * (gj % nj) + gi % ni;
	return blocks[slot];
}

// cell updates of one iteration (the level L makes 2^L steps)
long long Mesh::cellUpdatesPerIteration() {
	long long n = 0;
	for (int b = 0; b < nBlock; b++) {
		n += (long long)ni * nj * nk << blockLevel[b];
	}
	return n;
}

// cell updates of one iteration in the non empty blocks
long long Mesh::activeBlockCellUpdates() {
	long long n = 0;
	for (int b = 0; b < nBlock; b++) {
		if (!blocks[blockId[b]]->isBlockEmpty) {
			n += (long long)ni * nj * nk << blockLevel[b];
		}
	}
	return n;
}

// time step of the level 0
double Mesh::getCt() {
	return blocks[blockId[0]]->Ct * (1 << blockLevel[0]);
}
