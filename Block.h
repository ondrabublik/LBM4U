/*
 * File:   Block.h
 * Author: obublik
 *
 * Created on 8. února 2018, 15:53
 */

#ifndef BLOCK_H
#define BLOCK_H

#include <string>
#include <unordered_map>
#include <vector>
#include "Collision.h"

using namespace std;

// precision of the distribution functions (compile with -DLBM_FP32 for single precision)
#ifdef LBM_FP32
typedef float pop_t;
#else
typedef double pop_t;
#endif

// cell types stored in the boundary files
enum CellType {
	CT_FLUID = 0,
	CT_WALL = 1,
	CT_INLET = 2,   // velocity inlet (velInlet)
	CT_OUTLET = 3,  // pressure outlet (rhoOutlet)
	CT_MOVING = 4   // moving wall (velWall)
};

typedef unordered_map<string, string> ParamMap;

// sources of the data of ghost blocks (copies of the neighbourhood at another refinement level)
enum GhostSource {
	GHOST_NONE = 0,
	GHOST_FROM_COARSER = -1, // region of a block of the level - 1
	GHOST_FROM_FINER = 1     // region of blocks of the level + 1
};

// octree of blocks: level 0 blocks have the slot (nBi * nBj) * bk + nBi * bj + bi, finer levels are hashed
struct BlockGrid {
	int nBi = 0, nBj = 0, nBk = 0; // number of blocks of the level 0
	bool periodic[3] = { false, false, false };
	int maxLevel = 0;
	unordered_map<unsigned long long, int> slots; // levels >= 1

	static unsigned long long key(int level, int bi, int bj, int bk) {
		return ((unsigned long long)level << 60) | ((unsigned long long)bi << 40) | ((unsigned long long)bk << 20) | (unsigned long long)bj;
	}
	// wraps periodic block coordinates, false outside of the domain
	bool wrap(int level, int& bi, int& bj, int& bk) const;
	// slot of the block (it may be empty for the level 0), -1 outside of the domain or not existing
	int find(int level, int bi, int bj, int bk) const;
};

class Block {
public:
	double tau = 0.6;
	double vLB = 0;
	double rho0 = 0;          // reference density subtracted from the stored populations (DDF shifting)
	int collisionModel = COLLISION_BGK;
	bool mrtModel = false;
	CollisionOperator coll;
	double tauMin = 0.5001, tauMax = 1.5; // limits for variable relaxation time (non-Newtonian fluids)
	double smagorinsky = 0;              // Smagorinsky constant of the LES model (0 - laminar)
	Block** blocks = nullptr;
	int nBlock = 0;
	int nBlockMax = 0;
	int* blockMap = nullptr;
	int nBi = 0, nBj = 0, nBk = 0;
	int indBlock = 0;
	bool isBlockEmpty = false;
	bool isBlockTouched = false;
	int level = 0;              // refinement level, cell size H / ni / 2^level
	int bc[3] = { 0, 0, 0 };    // block coordinates at its level
	int ghost = GHOST_NONE;     // GhostSource of ghost blocks (they are not computed, only filled)
	const BlockGrid* grid = nullptr;
	string fileTag;             // identification of the block in the file names
	int* indAll = nullptr;
	int* indIn = nullptr;
	int* indBound = nullptr;
	// neighbour tables: interior cells use ind + offs[dir], boundary cells use the tables
	int* nbIndex = nullptr;     // per cell: row in nbBlockTab/nbLocalTab or -1
	int* nbBlockTab = nullptr;  // [nBound * nf] neighbour block (-1 outside of the domain)
	int* nbLocalTab = nullptr;  // [nBound * nf] neighbour local index
	int offs[27];
	int velocityModel = 0;
	double Ch = 0, Cu = 0, Ct = 0, Cp = 0, Cf = 0, Cv = 0, Cvis = 0, Cst = 0, Crho = 0;
	int* ex = nullptr, * ey = nullptr, * ez = nullptr, * opp = nullptr;
	double* w = nullptr;
	double PI = 3.14159265358979323846;
	int ni = 0, nj = 0, nk = 0, nf = 0, nAll = 0, nIn = 0, nBound = 0, ntot = 0, activeCellsNumBlock = 0;

	bool* boundary = nullptr;          // true for every non-fluid cell (wall, inlet, outlet, ...)
	unsigned char* cellType = nullptr; // CellType of each cell
	double* U = nullptr, * V = nullptr, * W = nullptr, * Rho = nullptr, * dRho = nullptr;
	pop_t* fc = nullptr, * f = nullptr;
	double* wallDistance = nullptr;
	double* wallDistanceGradient = nullptr;
	double* wallDistancePolynome = nullptr;

	double* G = nullptr, * G0 = nullptr; // gravity field
	double velWall[3] = { 0, 0, 0 };     // velocity of moving walls
	double velInlet[3] = { 0, 0, 0 };    // inlet velocity
	double rhoOutlet = 1.0;              // outlet density
	bool movingLidTop = true;            // domain boundary at z max moves with velWall
	bool xPeriodicity = false;
	bool yPeriodicity = false;
	bool zPeriodicity = false;

	Block() {}
	virtual ~Block();

	// common part of the constructors of derived blocks
	void initCommon(bool isFirst, Block** blocks, int nBlock, int nBlockMax, int* blockMap, int nBi, int nBj, int nBk, int indBlock, int level, int* indAll, int* indIn, int* indBound, int collisionModel, int nf, int* ex, int* ey, int* ez, double* w, int* opp, int ni, int nj, int nk, const ParamMap& par);
	static void countCells(int ni, int nj, int nk, int& nIn, int& nBound);
	static bool isInterior(int i, int j, int k, int ni, int nj, int nk);

	virtual void initZeroCondition();
	virtual void readInitCondition(int id, string simulPath);
	virtual void computeInitDistribution();
	void initBoundary(int id, string geomPath);
	// global cell coordinates at the level of the block
	void globalCoords(int ind, int& gi, int& gj, int& gk) const {
		gk = bc[2] * nk + ind / (ni * nj);
		gj = bc[1] * nj + (ind % (ni * nj)) / ni;
		gi = bc[0] * ni + ind % ni;
	}
	virtual void writeBlock(int id);
	virtual void propagation();
	// propagation of the listed cells only (ghost layers)
	virtual void propagateCells(const vector<int>& cells);
	virtual void collision();
	virtual void makro();
	// velocity written to the outputs (includes the half force correction)
	virtual void outputVelocity(int ind, double* u);
	// ghost cell ind takes the state of the coarser cell sInd of the block src (homogeneous explosion of the populations)
	virtual void prolongCell(int ind, const Block* src, int sInd);
	// ghost cell ind takes the average state of the 8 finer cells (src[o] is nullptr for missing cells)
	virtual void restrictCell(int ind, Block* const* src, const int* sInd);
	// mass added to the cells at the interface of the levels to keep the total mass (computed blocks with finer neighbours)
	double* massCorr = nullptr;
	// how the cell exchanges mass with the neighbours: 0 - none, 1 - by the populations, 2 - free surface mass exchange
	virtual int massExchangeKind(int ind) const { return boundary[ind] ? 0 : 1; }
	// clears the state of a ghost block that is not copied from the other level
	virtual void resetGhost() {}
	void local2globalShift(int indBlock, int ind, int& indBlockNeigh, int& indLocalNeigh, int ex, int ey, int ez);
	void computeShiftIndexes();
	inline int nbBlock(int ind, int dir) const {
		int b = nbIndex[ind];
		return b < 0 ? indBlock : nbBlockTab[nf * b + dir];
	}
	inline int nbLocal(int ind, int dir) const {
		int b = nbIndex[ind];
		return b < 0 ? ind + offs[dir] : nbLocalTab[nf * b + dir];
	}
	// type of the neighbour in direction dir, exterior of the domain is a wall (or the moving lid)
	int neighbourType(int ind, int dir, int& bl, int& loc) const;
	// population entering cell ind in direction opp[dir] from a non-fluid neighbour of type ctype
	double boundaryIncoming(int ind, int dir, int ctype) const;

	static int getParInt(const string& quantity, const ParamMap& par, bool mandatory, int defaultValue);
	static double getParDouble(const string& quantity, const ParamMap& par, bool mandatory, double defaultValue);
	static double* getParDoubleArray(const string& quantity, const ParamMap& par);
	static bool hasPar(const string& quantity, const ParamMap& par);
	// relaxation time of the Smagorinsky model for the non-equilibrium part fneq (Hou et al. 1996)
	double lesTau(double tau0, const double* fneq, double rho) const;
	void equilibrium(double* feq, double Rho, double dRho, double U, double V, double W);
	void equilibrium(double* feq, double Rho, double dRho, double U, double V, double W, double eps);
};

#endif /* BLOCK_H */
