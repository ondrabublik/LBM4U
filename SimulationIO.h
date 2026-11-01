/*
 * File:   SimulationIO.h
 *
 * Reading of the simulation parameter file (name = value per line, # or % starts a comment).
 */

#ifndef SIMULATIONIO_H
#define SIMULATIONIO_H

#include <string>
#include "Block.h"

std::string trim(const std::string& str);
ParamMap parseParams(const std::string& text);
ParamMap loadParamFile(const std::string& fileName);

#endif /* SIMULATIONIO_H */
