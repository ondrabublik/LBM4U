/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Mesh.h
 * Author: obublik
 *
 * Created on 6. února 2018, 13:53
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include "math.h"
#include <windows.h>
#include <string>
#include <unordered_map>

#ifndef MESH_H
#define MESH_H
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
	bool initHydroPress;
    int nBlock, nBlockMax;
    int **blockIndex;
    int * blockMap;
    Block **blocks;
    int* coarseLevels;
    int maxCoarseLevel;
    int ni, nj, nk, nf;
    int nBi, nBj, nBk;
    int *ex, *ey, *ez, *opp;
    double *w;
    double *G;
    int *indIn; // inlet indexes
    int *indBound; // boundary indexes
    int *indAll;
    double massTot;
    int nParticles;
    int maxParticlesInCell;
    int activeCellsNumTot;
    int nSmoothStep;
	int nIntIterParticle;

    // particles
    LagrangeParticles *particles;

	// fiberss
	LagrangeFibers* fibers;

    Mesh(string geomPath, string simulPath, unordered_map<string, string> par);
    ~Mesh();

    void model(int model, int** ex, int** ey, int** ez, double** w, int** opp);
    int loadStructure(string geomPath);
    void initIndexes();
    void initBlocks();
    void setWall(Block *b, int direction);
    void evalBlocks(int iter, int step);
    void changeData(int level);
    void exchange(Block *b, Block *bs, int indexTab);
    void exchangeVOF(Block *b, Block *bs, int indexTab);
    double getUmax();
	double getCt();
    void writeResults();
    void writeResultsToVTK(int iter);
	void writeResultsToOutFile();
    void writeBoundaryToVTK();
    void writeFreeSurfaceToVTK(int iter);
    double average(int*ei, int*ej, int*ek, int indBlock, int indLocal, double vofCenter);
	void generateDistanceFunction();
	void initHydrostaticPress();
};

#endif /* MESH_H */

