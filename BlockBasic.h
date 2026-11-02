/*
 * File:   BlockBasic.h
 * Author: obublik
 *
 * Created on 2. února 2018, 12:27
 *
 * Single phase LBM block (BGK/TRT/MRT) with Guo forcing.
 */

#ifndef BLOCKBGK_H
#define BLOCKBGK_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <cmath>
#include "Block.h"

class BlockBasic : public Block {
public:

	BlockBasic(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, bool isTauVariable, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par);

	~BlockBasic();

	void initZeroCondition();
	void readInitCondition(int i, string simulPath);
	void computeInitDistribution();
	// fused computation of the moments and the collision
	void collision();
	void propagation();
	void propagateCells(const vector<int>& cells);
	void makro();
	void writeBlock(int id);

private:
	void moments(int ind, double& rho, double* u);
};

#endif /* BLOCKBGK_H */
