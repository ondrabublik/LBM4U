/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Block.h
 * Author: obublik
 *
 * Created on 2. února 2018, 12:27
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <unordered_map>
#include "math.h"
#include "Block.h"

#ifndef BLOCKBGK_H
#define BLOCKBGK_H

class BlockBasic : public Block{
public:

    double *T;

    BlockBasic(bool isFirst, Block **blocks, int nBlock, int nBlockMax, int *blockIndex, int nBi, int nBj, int nBk, int indBlock, int coarseLevel, int* indAll, int* indIn, int* indBound, bool isMRTModel, bool isTauVariable, int nf, int *ex, int *ey, int *ez, double *w, int* opp, int ni, int nj, int nk, unordered_map<string, string> par);

    ~BlockBasic();

    void initZeroCondition();
    void readInitCondition(int i, string simulPath);
    void computeInitDistribution();
    void eval(int iter);
    void collision();
    void propagation();
    void makro();
    void defineMRTmatrix();
    void writeBlock(int id);
};

#endif /* BLOCKBGK_H */

