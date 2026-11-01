/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Block.h
 * Author: obublik
 *
 * Created on 8. února 2018, 15:53
 */

#ifndef BLOCK_H
#define BLOCK_H

using namespace std;

class Block {
public:
	double tau;
	double vLB;
	double rho0;
    bool mrtModel;
    Block **blocks;
    int nBlock;
    int nBlockMax;
    int *blockMap;
    int nBi, nBj, nBk;
    int indBlock;
    bool isBlockEmpty;
    bool isBlockTouched;
    int coarseLevel;
    int *indAll;
    int *indIn;
    int *indBound;
    int *shift;
    int *bShift;
    int velocityModel;
    double Ch, Cu, Ct, Cp, Cf, Cv, Cvis, Cst, Crho;
    int *ex, *ey, *ez, *opp;
    double *w;
    double PI = 3.14159265;
    int ni, nj, nk, nf, nAll, nIn, nBound, ntot, activeCellsNumBlock;
    
    bool *boundary;
	double* U, * V, * W, * Rho, * dRho;
    double *fc, *f, *fcOld;
	double* wallDistance;
	double* wallDistanceGradient;
	double* wallDistancePolynome;
    
    double* G, *G0; // gravity field
    double *velWall; // velocity at the wall
    bool xPeriodicity = 0;
    bool yPeriodicity = 0;
    bool zPeriodicity = 0;
    
    virtual void initZeroCondition();
    virtual void readInitCondition(int id, string simulPath);
    virtual void computeInitDistribution();
    void initBoundary(int id, string geomPath);
    virtual void writeBlock(int id);
    virtual void propagation();
    virtual void collision();
    virtual void makro();
    void local2globalShift(int indBlock, int ind, int& indBlockNeigh, int& indLocalNeigh, int ex, int ey , int ez);
    void computeShiftIndexes();
	int getParInt(string quantity, unordered_map<string, string> par, bool mandatory, int defaultValue);
	double getParDouble(string quantity, unordered_map<string, string> par, bool mandatory, double defaultValue);
	double* getParDoubleArray(string quantity, unordered_map<string, string> par);
	void equilibrium(double* feq, double Rho, double dRho, double U, double V, double W);
	void equilibrium(double* feq, double Rho, double dRho, double U, double V, double W, double eps);
};

#endif /* BLOCK_H */

