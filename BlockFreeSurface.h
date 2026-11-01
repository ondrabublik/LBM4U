/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BlockFreeSurface.h
 * Author: obublik
 *
 * Created on 9. února 2018, 9:22
 */

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <unordered_map>
#include "math.h"
#include "Block.h"

#ifndef BLOCKFREESURFACE_H
#define BLOCKFREESURFACE_H

class BlockFreeSurface : public Block {
public:
    
    double surfTension, wetAngleEq, wetAngleAdv, wetAngleRec;
    double massErr;
    double massBlock;
    double *mass, *vof, *vofS1, *vofS, *normal, *curvature;
    int nFilled, nEmptied, nAirToFs, nFluidToFs, nArtefacts1, nArtefacts2;
    int *flagMap;  // 1 - fluid, 2 - gas, 3 - interface
    int *surfaceLevels;
    int *filled, *emptied, *airToFs, *fluidToFs, *artefact1, *artefact2;
    int numFsCells;
    double *T, *M, *iM, *S, *T0, *T1;
    
    // nonewtonian fluid
	int viscosityModel;
    double *viscosity, *shearRate;
    double *fluidModelParam;

    BlockFreeSurface(bool isFirst, Block **blocks, int nBlock, int nBlockMax, int *blockMap, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int *indAll, int *indIn, int *indBound, bool isMRTModel, int viscosityModel, int nf, int *ex, int *ey, int *ez, double *w, int* opp, int ni, int nj, int nk, unordered_map<string, string> par);
    
    ~BlockFreeSurface();

    void initZeroCondition();
    void readInitCondition(int id, string simulPath);
	void writeOutCondition(int id, string simulPath);
    void computeInitDistribution();
    void initFlagMap();
    void iniVofAndMass();
    void collision();
    void propagation();
    void makro();
    void writeBlock(int id);
    void massExchange();
    void findFilledEmptied(int iter);
    void vicinityFilledEmptied();
    void newFsGeometry();
    void findArtefactsEmptied();
    void findArtefactsFilled();
    void destroyArtefacts();
    void massRedistribution();
    void redistribute(int ind, double massRedist);
    void calculateNormalVectors();
    void smoothVOF(int smoothStep, int nSmoothStep);
    void smooth3D(double* arrayin,double* arrayout,int kod);
    void calculateCurvature();
    void checkIfBlockEmpty();
    void touchNeighBlocks();
    void calculMassTotal();
    void defineMRTmatrix();
    void setWetAngle();
    void findSurfaceLevels(int level);
    double computeShearRate(int ind);
    void setMRTdiagonal(double tau);
    double nonNewtonianViscosityModel(double gamma);
};

#endif /* BLOCKFREESURFACEBGK_H */

