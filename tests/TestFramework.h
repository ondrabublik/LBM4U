/*
 * Minimal unit test framework (no external dependencies).
 */

#ifndef TESTFRAMEWORK_H
#define TESTFRAMEWORK_H

#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cmath>

struct TestFailure : public std::runtime_error {
	explicit TestFailure(const std::string& m) : std::runtime_error(m) {}
};

struct TestEntry {
	const char* name;
	void (*fn)();
};

std::vector<TestEntry>& testRegistry();

struct TestRegistrar {
	TestRegistrar(const char* name, void (*fn)()) {
		testRegistry().push_back({ name, fn });
	}
};

#define TEST(name) \
	static void name(); \
	static TestRegistrar registrar_##name(#name, name); \
	static void name()

#define CHECK(cond) \
	do { \
		if (!(cond)) { \
			std::ostringstream os_; \
			os_ << __FILE__ << ":" << __LINE__ << ": CHECK(" #cond ") failed"; \
			throw TestFailure(os_.str()); \
		} \
	} while (0)

#define CHECK_NEAR(a, b, tol) \
	do { \
		double va_ = (a), vb_ = (b), vt_ = (tol); \
		if (!(std::fabs(va_ - vb_) <= vt_)) { \
			std::ostringstream os_; \
			os_.precision(12); \
			os_ << __FILE__ << ":" << __LINE__ << ": CHECK_NEAR(" #a ", " #b ") failed: " << va_ << " vs " << vb_ << " (tol " << vt_ << ")"; \
			throw TestFailure(os_.str()); \
		} \
	} while (0)

#define CHECK_LT(a, b) 	do { 		double va_ = (a), vb_ = (b); 		if (!(va_ < vb_)) { 			std::ostringstream os_; 			os_.precision(10); 			os_ << __FILE__ << ":" << __LINE__ << ": CHECK_LT(" #a ", " #b ") failed: " << va_ << " >= " << vb_; 			throw TestFailure(os_.str()); 		} 	} while (0)

#define CHECK_THROWS(expr) \
	do { \
		bool thrown_ = false; \
		try { expr; } catch (...) { thrown_ = true; } \
		if (!thrown_) { \
			std::ostringstream os_; \
			os_ << __FILE__ << ":" << __LINE__ << ": expected exception from " #expr; \
			throw TestFailure(os_.str()); \
		} \
	} while (0)

#endif
