/*
 * File:   BlockFreeSurface.h
 * Author: obublik
 *
 * Created on 9. února 2018, 9:22
 *
 * Free surface LBM block (Koerner et al. 2005, Thuerey 2007).
 */

#ifndef BLOCKFREESURFACE_H
#define BLOCKFREESURFACE_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <cmath>
#include "Block.h"

// interface cell types for the mass exchange (Thuerey 2007)
enum InterfaceType {
	FS_STANDARD = 0,
	FS_NO_FLUID = 1,  // no fluid neighbours
	FS_NO_EMPTY = 2,  // no gas neighbours
	FS_ARTEFACT_FILLED = 3 // marked in findArtefactsFilled, becomes fluid
};

class BlockFreeSurface : public Block {
public:

	double surfTension, wetAngleEq, wetAngleAdv, wetAngleRec;
	double massErr;       // mass that could not be redistributed yet
	double massBlock;
	double* mass, * vof, * vofS1, * vofS, * normal, * curvature;
	int nFilled, nEmptied, nAirToFs, nFluidToFs, nArtefacts1, nArtefacts2;
	int* flagMap;  // 1 - fluid, 2 - gas, 3 - interface, 4 - not set
	int* surfaceLevels;
	int* filled, * emptied, * airToFs, * fluidToFs, * artefact1, * artefact2;
	unsigned char* fsType;
	int numFsCells;
	int curvatureModel;   // 0 - divergence of smoothed normals (robust), 1 - height function (accurate for smooth VOF, noisier in dynamics)
	double fsTolerance;   // overfilling/underfilling tolerance for the cell conversion
	int nAirToFsWarnings;
	std::vector<double> newCellState; // averaged rho, u, v, w for the gas cells becoming interface

	// nonewtonian fluid
	int viscosityModel;
	double* viscosity, * shearRate;
	double* fluidModelParam;
	double* tauCell; // relaxation time used in the last collision

	BlockFreeSurface();
	BlockFreeSurface(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int viscosityModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par);

	~BlockFreeSurface();

	void initZeroCondition();
	void readInitCondition(int id, string simulPath);
	void writeOutCondition(int id, string simulPath);
	void computeInitDistribution();
	void initFlagMap();
	void iniVofAndMass();
	void collision();
	void propagation();
	void propagateCells(const vector<int>& cells);
	void propagateCell(int ind);
	void makro();
	void writeBlock(int id);
	void computeInterfaceTypes();
	void massExchange();
	void findFilledEmptied(int iter);
	void vicinityFilledEmptied();
	void newFsGeometry();
	void newFsGeometryPrepare();
	void newFsGeometryApply();
	void findArtefactsEmptied();
	void findArtefactsFilled();
	void destroyArtefacts();
	void fixFluidGasContacts();
	void massRedistribution();
	void redistributeMassErr();
	void redistribute(int ind, double massRedist, int direction);
	void calculateNormalVectors();
	void smoothVOF(int smoothStep, int nSmoothStep);
	void smooth3D(double* arrayin, double* arrayout, int kod);
	void calculateCurvature();
	void checkIfBlockEmpty();
	void touchNeighBlocks();
	void calculMassTotal();
	void setWetAngle();
	void findSurfaceLevels(int surfLevel);
	// shear rate in lattice units computed locally from the non-equilibrium part
	double computeShearRate(int ind);
	double nonNewtonianViscosityModel(double gamma);
	void outputVelocity(int ind, double* u);

	// grid refinement: computed blocks receive mass from the other levels in massAcc and conversions in pending*,
	// ghost blocks collect there the mass and the conversions their cells got from the computed cells
	double* massAcc;
	std::vector<int> pendingAir, pendingFluid;
	void prolongCell(int ind, const Block* src, int sInd);
	void restrictCell(int ind, Block* const* src, const int* sInd);
	int massExchangeKind(int ind) const;
	void resetGhost();
	// mass exchanged by the interface cell ind with the fluid or interface neighbour nb in the direction p
	double linkExchange(int ind, int p, const BlockFreeSurface* nb, int loc) const;
	void applyMassAcc();
	void injectPending();

	// acceleration acting on the cell and fluid fraction (porosity) of the cell
	virtual void cellForce(int ind, double* acc, double& epsF);
	// relaxation time of the cell (non-Newtonian model)
	virtual double localTau(int ind);
	// collision with the exact difference method forcing
	void collideCell(int ind, const double* acc, double epsF, double tauLoc);

protected:
	void initFreeSurface(bool isFirst, int viscosityModel, const ParamMap& par);
	double curvatureDivergence(int ind);
	bool curvatureHeightFunction(int ind, double& curv);
	double vofAt(int ind, int dx, int dy, int dz);
	void koernerPopulation(int ind, int p);
};

#endif /* BLOCKFREESURFACE_H */
