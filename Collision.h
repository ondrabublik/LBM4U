/*
 * File:   Collision.h
 *
 * Collision operators in velocity space: BGK, TRT and MRT.
 * relax() returns delta = Omega * fneq, the post-collision state is f - delta.
 */

#ifndef COLLISION_H
#define COLLISION_H

#include <vector>

enum CollisionModel {
	COLLISION_BGK = 0,
	COLLISION_MRT = 1,
	COLLISION_TRT = 2
};

class CollisionOperator {
public:
	int model;
	int nf;
	const int* opp;
	double magic;       // TRT magic parameter Lambda (3/16 places bounce-back walls exactly halfway)
	bool mrtMagicSq;    // MRT: relate the rates of the odd moments to s_nu (TRT magic 3/16) so the wall position does not depend on tau
	std::vector<double> M, iM, S0; // MRT transformation, its inverse and the default diagonal
	std::vector<int> stressMoments; // MRT moments relaxed with s_nu
	std::vector<int> fluxMoments;   // odd MRT moments relaxed with s_q

	CollisionOperator();
	void init(int model, int nf, const int* ex, const int* ey, const int* ez, const int* opp, double magic, bool mrtMagicSq, double tau0);
	// delta[p] = (Omega * fneq)[p] for relaxation time tau
	void relax(const double* fneq, double tau, double* delta) const;
	double tauMinus(double tau) const;

private:
	double cachedTau;
	std::vector<double> cachedT; // M^-1 S M for cachedTau
	void buildMatrices(const int* ex, const int* ey, const int* ez);
	void buildT(double tau, std::vector<double>& T) const;
	void diagonal(double tau, double* S) const;
};

// Gauss-Jordan inversion of an n x n row-major matrix, returns false for singular input
bool invertMatrix(const std::vector<double>& A, std::vector<double>& inv, int n);

#endif /* COLLISION_H */
