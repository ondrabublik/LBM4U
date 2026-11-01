/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: obublik
 *
 * Created on 31. ledna 2018, 14:07
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "BlockBasic.h"
#include "Mesh.h"
#include <omp.h>
#include <unordered_map>
#include <windows.h>
#include <ctime>

using namespace std;

string trim(const string& str) {
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first) {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

void loadArgs(string &geomPath, string &simulPath) {

    string line, delimiter = " ";
    ifstream *myfile = new ifstream("args.txt");

    if (myfile->is_open()) {
        getline(*myfile, line);
        size_t pos = 0;
        pos = line.find(delimiter);
        stringstream gp, sp;
        gp << "simulations/" << line.substr(0, pos) << "/";
        geomPath = gp.str();
        sp << "simulations/" << line.substr(0, pos) << "/" << line.substr(pos + 1, line.length()) << "/";
        simulPath = sp.str();
        cout << "Geometry name: " << line.substr(0, pos) << endl;
        cout << "Simulation: " << line.substr(pos + 1, line.length()) << endl;
    }
}

unordered_map<string, string> loadParamFile(string simulPath) {
    unordered_map<string, string> param;

    string line, name, value, delimiter = "=";
    stringstream ss;
    ss << simulPath << "parameters.txt";
    ifstream *myfile = new ifstream(ss.str());

    if (myfile->is_open()) {
        while (getline(*myfile, line)) {
            size_t pos = 0;
            while ((pos = line.find(delimiter)) != string::npos) {
                name = trim(line.substr(0, pos));
                line.erase(0, pos + delimiter.length());
                value = trim(line);
            }

            param[name] = value;
        }
    }

    return param;
}

int main(int argc, char** argv) {
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

    // loading arguments
    string geomPath;
    string simulPath;
    loadArgs(geomPath, simulPath);

    // loading parameters
	Block b;
    unordered_map<string, string> par = loadParamFile(simulPath);
	int nThreads = b.getParInt("threads", par, false, 1);
    omp_set_num_threads(nThreads);
    int maxIter = b.getParInt("maxIter",par, true, 0);
    int saveRate = b.getParInt("saveRate",par, false, 1e10);
    int problemID = b.getParInt("problemID",par, true, 0);
    cout << "Problem ID: " << problemID << endl;
    cout << "Number of threads: " << nThreads << endl;
    cout << "Total number of iteration: " << maxIter << endl;
    cout << "Save rate: " << saveRate << endl;
    
	stringstream ss;
	ss << simulPath << "result";
	if (CreateDirectory((ss.str()).c_str(), NULL)) {
		cout << "Results directory was create." << endl;
	}

	stringstream ss2;
	ss2 << simulPath << "out";
	if (CreateDirectory((ss2.str()).c_str(), NULL)) {
		cout << "Out directory was create." << endl;
	}

    Mesh *mesh = new Mesh(geomPath, simulPath, par);
    mesh->initBlocks();
    mesh->writeBoundaryToVTK();
    cout << "" << endl;
    cout << "Computation start ..." << endl;
    ofstream outfile;
    outfile.open(simulPath + "result/residuum.txt", ios::out);
    
	double massTot;
	double Ct = mesh->getCt();
    for (int iter = 0; iter < maxIter; iter++) {
        clock_t begin = clock();
        switch (problemID) {
            case 0: // classic LBM
                mesh->evalBlocks(iter, 0);
                mesh->evalBlocks(iter, 1);
                break;
            case 1: // free surface LBM
                mesh->evalBlocks(iter, 2); // mass exchange, find filled/emptied cells
                mesh->evalBlocks(iter, 3);
                mesh->evalBlocks(iter, 4);
                mesh->evalBlocks(iter, 5); // smooth
                mesh->evalBlocks(iter, 6); // normals and wet angle
                mesh->evalBlocks(iter, 7); // curvature
                mesh->evalBlocks(iter, 0); // propagation
                mesh->evalBlocks(iter, 8); // macro, collision, new fs
                mesh->evalBlocks(iter, 9); // destroy artefacts
                mesh->evalBlocks(iter, 10); // mass redistribution
                mesh->evalBlocks(iter, 11); // check empty block flag
                break;
            case 2: // free surface LBM with particles
				mesh->evalBlocks(iter, 202);
                mesh->evalBlocks(iter, 2);
                mesh->evalBlocks(iter, 3);
                mesh->evalBlocks(iter, 4);
                mesh->evalBlocks(iter, 5);
                mesh->evalBlocks(iter, 6); // normals and wet angle
                mesh->evalBlocks(iter, 7); // curvature
                mesh->evalBlocks(iter, 0); // propagation
                mesh->evalBlocks(iter, 20); // particles
                mesh->evalBlocks(iter, 8); // macro, collision, new fs
                mesh->evalBlocks(iter, 9); // destroy artefacts
                mesh->evalBlocks(iter, 10);
                break;
			case 21: // free surface LBM with fibers
				mesh->evalBlocks(iter, 2);
				mesh->evalBlocks(iter, 3);
				mesh->evalBlocks(iter, 4);
				mesh->evalBlocks(iter, 5);
				mesh->evalBlocks(iter, 6); // normals and wet angle
				mesh->evalBlocks(iter, 7); // curvature
				mesh->evalBlocks(iter, 0); // propagation			
				mesh->evalBlocks(iter, 21); // fibers dynamic
				mesh->evalBlocks(iter, 8); // macro, collision, new fs
				mesh->evalBlocks(iter, 9); // destroy artefacts
				mesh->evalBlocks(iter, 10);
				break;
			case 3: // free surface LBM with euler particles - solve transport eq. using FVM
				mesh->evalBlocks(iter, 2);  // mass exchange, find filled/emptied cells
				mesh->evalBlocks(iter, 3);
				mesh->evalBlocks(iter, 4);
				mesh->evalBlocks(iter, 5);  // smooth
				mesh->evalBlocks(iter, 6);  // normals and wet angle
				mesh->evalBlocks(iter, 7);  // curvature
				mesh->evalBlocks(iter, 0);  // propagation
				mesh->evalBlocks(iter, 30); // compute GradP + interpolate fluid velocity
				mesh->evalBlocks(iter, 31); // compute new particle velocity
				mesh->evalBlocks(iter, 32); // smooth force + upwind
				mesh->evalBlocks(iter, 8);  // macro, collision, new fs
				mesh->evalBlocks(iter, 9);  // destroy artefacts
				mesh->evalBlocks(iter, 10); // mass redistribution
				mesh->evalBlocks(iter, 11); // check empty block flag
				break;
			case 4: // free surface LBM with euler particles - solve transport eq. using LBM
				mesh->evalBlocks(iter, 402);
				mesh->evalBlocks(iter, 2);  // mass exchange, find filled/emptied cells
				mesh->evalBlocks(iter, 3);
				mesh->evalBlocks(iter, 4);
				mesh->evalBlocks(iter, 5);  // smooth
				mesh->evalBlocks(iter, 6);  // normals and wet angle
				mesh->evalBlocks(iter, 7);  // curvature
				mesh->evalBlocks(iter, 0);  // propagation
				mesh->evalBlocks(iter, 40); // compute GradP + interpolate fluid velocity  
				mesh->evalBlocks(iter, 41); // compute new particle velocity            
				mesh->evalBlocks(iter, 42); // smooth force + propagation of scalar     
				mesh->evalBlocks(iter, 43); // compute concentration                    				
				mesh->evalBlocks(iter, 45); // collision of scalar			
				mesh->evalBlocks(iter, 8);  // macro, collision, new fs
				mesh->evalBlocks(iter, 9);  // destroy artefacts
				mesh->evalBlocks(iter, 10); // mass redistribution
				mesh->evalBlocks(iter, 11); // check empty block flag			
				break;
        }
        clock_t end = clock();
        if (iter % saveRate == 0) {
			double uMax = mesh->getUmax();
			double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            double massTot, activeCells;
			switch (problemID) {
				case 0:
					cout << "Itteration: " << iter << ", t = " << Ct * iter << ", uMax = " << uMax << ", Mlups = " << mesh->nBlock*mesh->ni* mesh->nj* mesh->nk / elapsed_secs / 1e6 <<endl;
					outfile << iter << " " << uMax << "\n";
					break;
				case 2:
					mesh->particles->getParticlesInfo();
				case 21:
				case 1:
				case 3:
                    mesh->evalBlocks(iter, 105);   // mass total
                    massTot = mesh->massTot;
                    activeCells = mesh->activeCellsNumTot;
                    cout << "Itteration: " << iter << ", t = " << Ct * iter << ", uMax = " << uMax << ", mass = " << massTot << ", Mlups = " << activeCells / elapsed_secs / 1e6 << endl;
                    outfile << iter << " " << uMax << " " << massTot << "\n";
                    mesh->writeFreeSurfaceToVTK(iter);
                    break;
				case 4:
					mesh->evalBlocks(iter, 105);   // mass total
					massTot = mesh->massTot;
                    activeCells = mesh->activeCellsNumTot;
                    int nParticles = mesh->nParticles;
			        int maxParticlesInCell = mesh->maxParticlesInCell;
					cout << "Itteration: " << iter << ", t = " << Ct * iter << ", uMax = " << uMax << ", mass = " << massTot << ", nP = " << nParticles << ", maxP = " << maxParticlesInCell << ", Mlups = " << activeCells / elapsed_secs / 1e6 << endl;
					outfile << iter << " " << uMax << " " << massTot << "\n";
					mesh->writeFreeSurfaceToVTK(iter);
					break;
			}
			mesh->writeResultsToVTK(iter);
        }
    }
    mesh->writeResultsToVTK(maxIter);
	mesh->writeResultsToOutFile();
	if (problemID == 1 || problemID == 2 || problemID == 21 || problemID == 3 || problemID == 4) {
		mesh->writeFreeSurfaceToVTK(maxIter);
	}
    delete mesh;
    outfile.close();

    cout << "Computation was sucesfull!" << endl;
    return 0;
}
