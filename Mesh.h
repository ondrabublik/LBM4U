/*
 * File:   Mesh.h
 * Author: obublik
 *
 * Created on 6. února 2018, 13:53
 */

#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <array>
#include "Block.h"
#include "BlockBasic.h"
#include "BlockFreeSurface.h"
#include "BlockFreeSurfaceLagrangeParticles.h"
#include "BlockFreeSurfaceLagrangeFibers.h"
#include "BlockFreeSurfaceEulerParticlesFVM.h"
#include "BlockFreeSurfaceEulerParticlesLBM.h"
#include "LagrangeParticles.h"
#include "LagrangeFibers.h"

using namespace std;

class Mesh {
public:
	string geomPath;
	string simulPath;
	int problemID;
	int viscosityModel;
	int collisionModel;
	bool initHydroPress;
	int nBlock, nBlockMax;
	int** blockIndex;
	int* blockId;     // global index of the i-th block
	int* blockMap;
	Block** blocks;
	int ni, nj, nk, nf;
	int nIn, nBound;
	int nBi, nBj, nBk;

	// octree refinement (static): blocks of the level L have the cell size H / ni / 2^L and make 2^L steps per iteration
	BlockGrid grid;
	int maxLevel;
	int nGhost;
	vector<vector<int> > levelBlocks;    // slots of the computed blocks of each level
	vector<int> blockLevel;              // level of the i-th block (blockIndex holds its coordinates at the level)
	// ghost block: copy of the neighbourhood of the blocks of its level held by blocks of another level
	struct GhostBlock {
		int slot;
		int source;                      // GHOST_FROM_COARSER or GHOST_FROM_FINER
		int bc[3];                       // block coordinates
		vector<int> srcBlock, srcLocal;  // per cell 1 (coarser) or 8 (finer) source cells, srcBlock -1 if missing
		vector<int> near1, near2;        // cells touching a computed cell of the level, at most 2 cells from it
		vector<int> halo;                // cells read by the stencils of the computed cells (at most GHOST_HALO cells from them)
		vector<int> nbReal;              // computed blocks of the level in the 26-neighbourhood
		int stamp = -1;                  // explosion of the coarser state the ghost holds
	};
	vector<vector<GhostBlock> > levelGhosts;
	int curLevel;                        // level processed by evalBlocks, -1 for all levels
	vector<int> explodeStamp;            // number of explosions to the ghosts of each level
	int* ex, * ey, * ez, * opp;
	double* w;
	int* indIn; // interior indexes
	int* indBound; // boundary indexes
	int* indAll;
	double massTot;
	int nParticles;
	int maxParticlesInCell;
	long long activeCellsNumTot;
	int nSmoothStep;
	int nIntIterParticle;
	double gravityRotationOmega;
	bool vtkBinary;
	bool compressResults; // zlib compressed VTK arrays and deflated fields_<iter>.bin
	bool writeFields;
	vector<pair<double, string> > vtkFiles, vtkFsFiles;

	// particles
	LagrangeParticles* particles;

	// fiberss
	LagrangeFibers* fibers;

	Mesh(string geomPath, string simulPath, const ParamMap& par);
	~Mesh();

	bool isFreeSurface() const;
	void model(int model, int** ex, int** ey, int** ez, double** w, int** opp);
	int loadStructure(string geomPath);
	void initIndexes();
	void initBlocks();
	void evalBlocks(int iter, int step);
	void timeStep(int iter);
	// block and local index of a cell given by global indexes of the level 0, nullptr if no computed block holds it
	Block* cellAt(int gi, int gj, int gk, int& ind);
	// block (computed or ghost) holding the cell given by global indexes of the level, nullptr if none
	Block* cellAtLevel(int level, int gi, int gj, int gk, int& ind);
	double getUmax();
	double getCt();
	long long cellUpdatesPerIteration();
	long long activeBlockCellUpdates();
	void writeResultsToVTK(int iter);
	void writeResultsToOutFile();
	void writeBoundaryToVTK();
	void writeFreeSurfaceToVTK(int iter);
	void writeFieldsBinary(int iter);
	void generateDistanceFunction();
	void initHydrostaticPress();

	// runs fn(block) for every (non empty) block of the current level (all levels if curLevel < 0) in parallel
	template <class F> void forBlocks(F fn, bool onlyActive = true) {
		if (curLevel < 0) {
#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < nBlock; i++) {
				Block* b = blocks[blockId[i]];
				if (!onlyActive || !b->isBlockEmpty) {
					fn(b);
				}
			}
			return;
		}
		const vector<int>& list = levelBlocks[curLevel];
		int n = (int)list.size();
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < n; i++) {
			Block* b = blocks[list[i]];
			if (!onlyActive || !b->isBlockEmpty) {
				fn(b);
			}
		}
	}
	// a ghost block is needed only next to a non empty computed block (empty free surface blocks are skipped)
	bool ghostActive(const GhostBlock& g) const {
		for (int slot : g.nbReal) {
			if (!blocks[slot]->isBlockEmpty) {
				return true;
			}
		}
		return false;
	}
	// runs fn(ghost) for the (active, all with all = true) ghost blocks of the level with the given source (0 - all) in parallel
	template <class F> void forGhosts(int level, int source, F fn, bool all = false) {
		vector<GhostBlock>& list = levelGhosts[level];
		int n = (int)list.size();
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < n; i++) {
			if ((source == 0 || list[i].source == source) && (all || ghostActive(list[i]))) {
				fn(list[i]);
			}
		}
	}

private:
	Block* createBlock(int indBlock, int level, bool isFirst, const ParamMap& par);
	void createGhosts(const vector<array<int, 4> >& list);
	void initGhostMaps();
	void initGhostNear();
	void syncGhostGeometry();
	// copies the state of the level - 1 to the ghosts of the level (populations: homogeneous explosion)
	void explodeToGhosts(int level, bool all = false);
	void explodeGhost(GhostBlock& g);
	// ghosts that became active after the explosion of the current coarser step
	void refreshGhosts(int level);
	// fills the ghosts of the level from the level + 1 (populations: coalescence of the streamed ones)
	void coalesceToGhosts(int level, bool streamed, bool all = false);
	void syncGhosts();
	void computeLost(int level);
	void accountFsExchange(int level);
	void forwardGhostRequests(int level);
	void forwardGhostMass(int level);
	void activateFromGhosts(int level);
	void advanceFreeSurfaceLevel(int level, int iterLevel);
	void accountInterface(int level);
	void applyMassCorr(int level);
	void advanceLevel(int level, int iterLevel);
	string outputName(const string& prefix, int iter, const string& ext);
	void writePvd(const string& fileName, const vector<pair<double, string> >& files);
};

#endif /* MESH_H */
