/*
 * File:   SimulationIO.cpp
 */

#include "SimulationIO.h"
#include <fstream>
#include <sstream>

using namespace std;

string trim(const string& str) {
	const char* ws = " \t\r\n";
	size_t first = str.find_first_not_of(ws);
	if (first == string::npos) {
		return "";
	}
	size_t last = str.find_last_not_of(ws);
	return str.substr(first, last - first + 1);
}

ParamMap parseParams(const string& text) {
	ParamMap param;
	istringstream in(text);
	string line;
	while (getline(in, line)) {
		size_t c = line.find_first_of("#%");
		if (c != string::npos) {
			line = line.substr(0, c);
		}
		size_t pos = line.find('=');
		if (pos == string::npos) {
			continue;
		}
		string name = trim(line.substr(0, pos));
		if (!name.empty()) {
			param[name] = trim(line.substr(pos + 1));
		}
	}
	return param;
}

ParamMap loadParamFile(const string& fileName) {
	ifstream myfile(fileName);
	if (!myfile.is_open()) {
		return ParamMap();
	}
	stringstream ss;
	ss << myfile.rdbuf();
	return parseParams(ss.str());
}
