#include "../graph.h"
#include "helper.hpp"

#include <iostream>
#include <cmath>
#include <tuple>
#include <iomanip>

using namespace std;


// any set of orthogonal polynomials satisfies a three term recurrence relationship 
// p_n+1(x) = (a_n * x + b_n) * p_n(x) - c_n * p_n-1(x) with a_n, c_n > 0
void recurrenceRelation(Vec& Pnext, const Vec& Pnow, const Vec& Plast, Func a, Func b, Func c, int n = 0) {
	Pnext[0] = b(n) * Pnow[0] - c(n) * Plast[0];
	for (int i = 0; i < n - 1; ++i) {
		Pnext[i + 1] = a(n) * Pnow[i] + b(n) * Pnow[i + 1] - c(n) * Plast[i + 1];
	}
	Pnext[n] = a(n) * Pnow[n - 1] + b(n) * Pnow[n];
	Pnext[n + 1] = a(n) * Pnow[n];
}

Vec generateOrthogonalPolynomials(Func a, Func b, Func c, int n, const Vec& p0 = { 1.0 }) {
	Vec Plast(n + 1, 0.0), Pnow(n + 1, 0.0), Pnext(n + 1, 0.0);
	Pnow[0] = p0[0];

	for (int i = 0; i < n; ++i) {
		recurrenceRelation(Pnext, Pnow, Plast, a, b, c, i);
		std::swap(Plast, Pnow);
		std::swap(Pnow, Pnext);
	}

	return Pnow;
}




// mu0 is the integral of the weight function in the integration interval
std::tuple<Vec, Vec> GaussQuadRule(Vec a, Vec b, double mu0) {
	int n = a.size();

	// simulate 1-based indexing
	a.insert(a.begin(), 0.0);
	b.insert(b.begin(), 0.0);

	int i, j, k, m, ml;
	double norm, eps, ct, st, f, q, aa, aj, a2, eigmax,
		lambda, lambda1, lambda2, rho, r, det, bi, bj, b2, wj, cj;

	Vec w(n + 1, 0.0);

	// find the maximum row sum norm and initialize w
	b[0] = 0.0;
	norm = 0.0;

	for (i = 1; i <= n - 1; ++i) {
		norm = std::max(norm, std::abs(a[i]) + std::abs(b[i - 1]) + std::abs(b[i]));
	}
	norm = std::max(norm, std::abs(a[n]) + std::abs(b[n - 1]));
	w[1] = 1.0;

	m = n;
	eps = norm * 8e-13; // relative zero tolerance
	lambda = lambda1 = lambda2 = rho = norm;


	while (m > 0) {
		i = k = ml = m - 1;

		if (std::abs(b[ml]) <= eps) {
			w[m] = mu0 * w[m] * w[m];
			rho = std::min(lambda1, lambda2);
			m = ml;
			continue;
		}

		// small off-diagonal element means matrix can be split
		while (i > 0 && std::abs(b[i]) > eps) {
			k = i--;
		}

		b2 = b[ml] * b[ml];
		det = std::sqrt((a[ml] - a[m]) * (a[ml] - a[m]) + 4.0 * b2);
		aa = a[ml] + a[m];
		lambda2 = 0.5 * (aa >= 0 ? aa + det : aa - det);
		lambda1 = (a[ml] * a[m] - b2) / lambda2;
		eigmax = std::max(lambda1, lambda2);

		if (std::abs(eigmax - rho) <= 0.125 * std::abs(eigmax)) {
			lambda = rho = eigmax;
		} else {
			rho = eigmax;
		}

		// transform block from k to m
		cj = b[k];
		b[k - 1] = a[k] - lambda;
		for (j = k; j <= ml; ++j) {
			r = std::sqrt(cj * cj + b[j - 1] * b[j - 1]);
			st = cj / r;
			ct = b[j - 1] / r;
			aj = a[j];
			b[j - 1] = r;
			cj = b[j + 1] * st;
			b[j + 1] = -b[j + 1] * ct;
			f = aj * ct + b[j] * st;
			q = b[j] * ct + a[j + 1] * st;
			a[j] = f * ct + q * st;
			b[j] = f * st - q * ct;
			wj = w[j];
			a[j + 1] = aj + a[j + 1] - a[j];
			w[j] = wj * ct + w[j + 1] * st;
			w[j + 1] = wj * st - w[j + 1] * ct;
		}

		b[k - 1] = 0;
	}

	// diagonal entries will converge to the eigenvalues
	return { a, w };
}





int main() {

	Graph graph;

	cout << std::fixed << std::setprecision(15);

	// test with Legendre polynomials
	Func a = [](double n) { return (2.0 * n + 1.0) / (n + 1.0); };
	Func b = constFunc(0.0);
	Func c = [](double n) { return n / (n + 1.0); };

	int n = 5;

	// will just be used for plotting
	Vec coeffs = generateOrthogonalPolynomials(a, b, c, n);

	Vec diag(n, 1.0);
	Vec offdiag(n - 1);
	for (int i = 0; i < n - 1; ++i) {
		diag[i] = -b(i) / a(i);
		offdiag[i] = std::sqrt(c(i + 1) / (a(i) * a(i + 1)));
	}
	diag[n - 1] = -b(n - 1) / a(n - 1);


	auto [roots, weights] = GaussQuadRuleNoGoto(diag, offdiag, 2.0);
	cout << roots << "\n" << weights << "\n";

	plotPolynomial(coeffs, &graph, -1.0, 1.0, 100);

	for (size_t i = 1; i <= n; ++i) {
		graph.addPoint(Point(roots[i], 0.0, 2, olc::RED));
		cout << "Poly at root " << i << ": " << evaluatePolynomial(coeffs, roots[i]) << "\n";
	}

	graph.waitFinish();

	return 0;
}
