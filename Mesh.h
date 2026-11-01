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
	// block and local index of a cell given by global indexes, nullptr if the block does not exist
	Block* cellAt(int gi, int gj, int gk, int& ind);
	double getUmax();
	double getCt();
	void writeResultsToVTK(int iter);
	void writeResultsToOutFile();
	void writeBoundaryToVTK();
	void writeFreeSurfaceToVTK(int iter);
	void writeFieldsBinary(int iter);
	void generateDistanceFunction();
	void initHydrostaticPress();

	// runs fn(block) for every (non empty) block in parallel
	template <class F> void forBlocks(F fn, bool onlyActive = true) {
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nBlock; i++) {
			Block* b = blocks[blockId[i]];
			if (!onlyActive || !b->isBlockEmpty) {
				fn(b);
			}
		}
	}

private:
	string outputName(const string& prefix, int iter, const string& ext);
	void writePvd(const string& fileName, const vector<pair<double, string> >& files);
};

#endif /* MESH_H */
