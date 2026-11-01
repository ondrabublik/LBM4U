#include "TestFramework.h"
#include <iostream>
#include <chrono>
#include <cstring>
#include <omp.h>

std::vector<TestEntry>& testRegistry() {
	static std::vector<TestEntry> r;
	return r;
}

int main(int argc, char** argv) {
	// the test cases are small, many threads only add overhead
	omp_set_num_threads(4);
	int failed = 0, run = 0;
	for (auto& t : testRegistry()) {
		if (argc > 1 && std::strstr(t.name, argv[1]) == nullptr) {
			continue;
		}
		run++;
		auto t0 = std::chrono::steady_clock::now();
		try {
			t.fn();
			double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
			std::cout << "[PASS] " << t.name << " (" << s << " s)" << std::endl;
		}
		catch (const std::exception& e) {
			failed++;
			std::cout << "[FAIL] " << t.name << ": " << e.what() << std::endl;
		}
	}
	std::cout << run - failed << "/" << run << " tests passed" << std::endl;
	return failed == 0 ? 0 : 1;
}
