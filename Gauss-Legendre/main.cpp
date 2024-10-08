#include "../helper.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <iomanip>

using namespace std;


// Rodrigues formula (works untill about n = 40, then rapidly looses precision (big numbers))
Vec LegendreCoefficientsRodrigues(int n) {

	Vec coeffs(2 * n + 1);

	// binomial expansion for (x^2 - 1)^n
	double sign = n & 1 ? -1.0 : 1.0;
	for (int k = 0; k <= n; ++k) {
		coeffs[2 * k] += sign * binomialCoefficient(n, k);
		sign = -sign;
	}

	// derive n times
	for (int i = 0; i < n; ++i) {
		 // take advantage of the fact that we start with only even powers having nonzero coefficients
		for (int j = 1 + (~i & 1); j < coeffs.size() - i; j += 2) {
			coeffs[j - 1] = coeffs[j] * j;
			coeffs[j] = 0.0;
		}
		coeffs[i & 1] = 0.0;
	}

	coeffs.resize(n + 1);

	double term1 = 1.0 / static_cast<double>(factorial(std::min(20, n)));
	double term2 = 1.0 / std::pow(2.0, n);

	// this is because factorial(21) overflows uint64_t
	// calculating the factorial with doubles (or via std::tgamma(n + 1)) seem to reduce precision, but it's a little more efficient
	for (double j = 21.0; j <= n; ++j) {
		term1 /= j;
	}

	for (int i = n & 1; i < coeffs.size(); i += 2) {
		coeffs[i] *= term1;
		coeffs[i] *= term2;
	}

	return coeffs;
}


// uses the recurrent relation (n+1) * P_n+1 = (2n + 1)x * P_n - n * P_n-1, with P_0 = 1 and P_1 = x
// using this over Rodrigues formula gives much greater numerical stability, being able to go beyond n = 100
// this is also more easily adaptable. For example, generating Chebyshev polynomials would be almost the same,
// but with term1 = -1.0 and term2 = 2.0 (not tested but should be true)
void recurrenceRelation(Vec& Pnext, const Vec& Pnow, const Vec& Plast, int n = 0) {
	double term1 = -n / (n + 1.0);
	double term2 = (2.0 * n + 1.0) / (n + 1.0);

	Pnext[0] = term1 * Plast[0];
	for (int i = 0; i < n - 1; ++i) {
		Pnext[i + 1] = term2 * Pnow[i] + term1 * Plast[i + 1];
	}
	Pnext[n] = term2 * Pnow[n - 1];
	Pnext[n + 1] = term2 * Pnow[n];
}

Vec LegendreCoefficients(int n) {

	Vec Plast(n + 1, 0.0), Pnow(n + 1, 0.0), Pnext(n + 1, 0.0);

	Pnow[0] = 1.0;

	for (int i = 0; i < n; ++i) {
		recurrenceRelation(Pnext, Pnow, Plast, i);
		std::swap(Plast, Pnow);
		std::swap(Pnow, Pnext);
	}

	return Pnow;
}



Vec estimateLegendreRoots(int n) {
	double term1 = PI / (n + 0.5);
	double term2 = 1.0 / (8.0 * static_cast<double>(n * n * n)) - 1.0 / (8.0 * static_cast<double>(n * n)) + 1.0;

	Vec roots(n, 0.0);
	for (int i = 0; i < n >> 1; ++i) {
		double theta = (static_cast<double>(n) - i - 0.25) * term1;
		roots[i] = term2 * std::cos(theta);
		roots[n - i - 1] = -roots[i]; // roots are symmetrical around 0
	}

	return roots;
}


tuple<Vec, Vec> getLegendreCoefficientsAndRoots(int n) {
	Vec coeffs = LegendreCoefficients(n);
	Vec roots = estimateLegendreRoots(n);

	for (int i = 0; i < n >> 1; ++i) {
		roots[i] = NewtonPolynomial(coeffs, roots[i]);
		roots[n - i - 1] = -roots[i];
	}

	return { coeffs, roots };
}



tuple<Vec, Vec> GaussLegendrePointsAndWeights(int n) {
	auto [coeffs, roots] = getLegendreCoefficientsAndRoots(n);

	Vec weights(n);

	for (int i = 0; i < n >> 1; ++i) {
		double derivative = evaluatePolynomialDerivative(coeffs, roots[i]);

		weights[i] = 2.0 / ((1.0 - roots[i] * roots[i]) * derivative * derivative);
		weights[n - i - 1] = weights[i];
	}

	// handles odd n (the root in the middle is 0)
	if (n & 1) {
		weights[n >> 1] = 2.0 / (coeffs[1] * coeffs[1]);
	}

	return { roots, weights };
}


double GaussLegendreQuadrature(const Func& f, double a, double b, int n) {
	auto [Xs, Ws] = GaussLegendrePointsAndWeights(n);

	double I = 0.0;
	for (int i = 0; i < n; ++i) {
		I += Ws[i] * f(0.5 * (a + b + Xs[i] * (a - b)));
	}

	return I * (b - a) * 0.5;
}







int main() {

	cout << std::fixed << std::setprecision(15);

	cout << GaussLegendreQuadrature([](double x) { return std::exp(-x * x); }, -10.0, 10.0, 40) << "\n";

	return 0;
}
