/*
 * File:   main.cpp
 * Author: obublik
 *
 * Created on 31. ledna 2018, 14:07
 *
 * Usage:
 *   cLBMnet                               reads "geometry simulation" from args.txt
 *   cLBMnet <geometry> <simulation>       case in simulations/<geometry>/<simulation>/
 *   cLBMnet --paths <geomDir> <simulDir>  explicit directories
 *
 * The run stops after maxIter iterations or when the file STOP appears in the simulation directory.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <cmath>
#include <omp.h>
#include "BlockBasic.h"
#include "Mesh.h"
#include "SimulationIO.h"

using namespace std;

static void loadArgs(string& geomPath, string& simulPath) {
	ifstream myfile("args.txt");
	string geom, sim;
	if (!(myfile >> geom >> sim)) {
		throw runtime_error("args.txt not found or invalid (expected: geometryName simulationName)");
	}
	geomPath = "simulations/" + geom + "/";
	simulPath = "simulations/" + geom + "/" + sim + "/";
	cout << "Geometry name: " << geom << endl;
	cout << "Simulation: " << sim << endl;
}

static string withSlash(string p) {
	if (!p.empty() && p.back() != '/' && p.back() != '\\') {
		p += "/";
	}
	return p;
}

static int run(int argc, char** argv) {
	cout << "" << endl;
	cout << "L       BBBBB   MM    MM     44   U     U" << endl;
	cout << "L       B    B  M M  M M    4 4   U     U" << endl;
	cout << "L       BBBBB   M  MM  M   4  4   U     U" << endl;
	cout << "L       B    B  M      M  444444  U     U" << endl;
	cout << "LLLLLL  BBBBB   M      M      4    UUUUU " << endl;
	cout << "" << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << "Welcome to the Lattice Boltzmann computational software!" << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << "" << endl;

	string geomPath;
	string simulPath;
	if (argc == 4 && string(argv[1]) == "--paths") {
		geomPath = withSlash(argv[2]);
		simulPath = withSlash(argv[3]);
	}
	else if (argc == 3) {
		geomPath = "simulations/" + string(argv[1]) + "/";
		simulPath = geomPath + string(argv[2]) + "/";
	}
	else {
		loadArgs(geomPath, simulPath);
	}
	cout << "Geometry path: " << geomPath << endl;
	cout << "Simulation path: " << simulPath << endl;

	ParamMap par = loadParamFile(simulPath + "parameters.txt");
	if (par.empty()) {
		throw runtime_error("parameters file " + simulPath + "parameters.txt not found or empty");
	}
	int nThreads = Block::getParInt("threads", par, false, 1);
	omp_set_num_threads(nThreads);
	int maxIter = Block::getParInt("maxIter", par, true, 0);
	int saveRate = max(Block::getParInt("saveRate", par, false, 1000000000), 1);
	int problemID = Block::getParInt("problemID", par, true, 0);
	cout << "Problem ID: " << problemID << endl;
	cout << "Number of threads: " << nThreads << endl;
	cout << "Total number of iteration: " << maxIter << endl;
	cout << "Save rate: " << saveRate << endl;

	filesystem::create_directories(simulPath + "result");
	filesystem::create_directories(simulPath + "out");
	string stopFile = simulPath + "STOP";
	filesystem::remove(stopFile);

	Mesh* mesh = new Mesh(geomPath, simulPath, par);
	mesh->initBlocks();
	mesh->writeBoundaryToVTK();
	cout << "" << endl;
	cout << "Computation start ..." << endl;
	ofstream outfile(simulPath + "result/residuum.txt", ios::out);

	double Ct = mesh->getCt();
	long long totalCells = mesh->cellUpdatesPerIteration();
	double tStart = omp_get_wtime();
	double tLast = tStart;
	int iterLast = 0;
	int exitCode = 0;
	int iter = 0;
	bool velocityWarning = false;
	for (iter = 0; iter < maxIter; iter++) {
		mesh->timeStep(iter);
		if (iter % saveRate == 0) {
			double now = omp_get_wtime();
			double elapsed = max(now - tLast, 1e-9);
			int iters = max(iter - iterLast, 1);
			double uMax = mesh->getUmax();
			double massTot = 0;
			long long activeCells = totalCells;
			if (mesh->isFreeSurface()) {
				mesh->evalBlocks(iter, 105);   // mass total
				massTot = mesh->massTot;
				activeCells = mesh->activeCellsNumTot;
			}
			double mlups = (iter == 0) ? 0.0 : (double)activeCells * iters / elapsed / 1e6;
			cout << "Iteration: " << iter << ", t = " << Ct * iter << ", uMax = " << uMax;
			if (mesh->isFreeSurface()) {
				cout << ", mass = " << massTot;
			}
			if (problemID == 4) {
				cout << ", nP = " << mesh->nParticles << ", maxP = " << mesh->maxParticlesInCell;
			}
			cout << ", MLUPS = " << mlups;
			if (mesh->isFreeSurface() && iter > 0) {
				// MLUPS counts the liquid cells only, the work is done for all cells of the non empty blocks
				cout << " (active blocks " << (double)mesh->activeBlockCellUpdates() * iters / elapsed / 1e6 << " Mcells/s)";
			}
			cout << endl;
			if (problemID == 2) {
				mesh->particles->getParticlesInfo();
			}
			outfile << iter << " " << uMax;
			if (mesh->isFreeSurface()) {
				outfile << " " << massTot;
			}
			outfile << "\n";
			outfile.flush();
			if (mesh->isFreeSurface()) {
				mesh->writeFreeSurfaceToVTK(iter);
			}
			mesh->writeResultsToVTK(iter);
			tLast = omp_get_wtime();
			iterLast = iter;
			double uLattice = uMax / mesh->blocks[mesh->blockId[0]]->Cu;
			if (uLattice > 0.3 && !velocityWarning) {
				cout << "warning: lattice velocity " << uLattice << " is high (Mach number), increase the characteristic velocity uMax or refine the mesh" << endl;
				velocityWarning = true;
			}
			if (std::isnan(uMax)) {
				cout << "Error: the solution diverged (NaN velocity)" << endl;
				exitCode = 3;
				break;
			}
			if (filesystem::exists(stopFile)) {
				cout << "Stop requested" << endl;
				filesystem::remove(stopFile);
				iter++;
				break;
			}
		}
	}
	if (exitCode == 0) {
		mesh->writeResultsToVTK(iter);
		mesh->writeResultsToOutFile();
		if (mesh->isFreeSurface()) {
			mesh->writeFreeSurfaceToVTK(iter);
		}
	}
	cout << "Total time: " << omp_get_wtime() - tStart << " s" << endl;
	delete mesh;

	if (exitCode == 0) {
		cout << "Computation was sucesfull!" << endl;
	}
	return exitCode;
}

int main(int argc, char** argv) {
	try {
		return run(argc, argv);
	}
	catch (const std::exception& e) {
		cout << "Error: " << e.what() << endl;
		return 1;
	}
}
