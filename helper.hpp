#ifndef HELPER_HPP
#define HELPER_HPP

#include <vector>
#include <iostream>
#include <functional>
#include <cstdint>

constexpr const double PI = 3.14159265359;

using Vec = std::vector<double>;
using Mat = std::vector<std::vector<double>>;
using Func = std::function<double(double)>;

Mat makeMat(size_t n, size_t m, double p = 0.0) {
	return std::vector<std::vector<double>>(n, std::vector<double>(m, p));
}

Func constFunc(double n = 0.0) {
	return [n](double x) { return n; };
}

std::ostream& operator << (std::ostream& os, const Vec& v) {

	os << "(";

	for (size_t i = 0; i < v.size(); ++i) {
		os << v[i];

		if (i + 1 < v.size()) {
			os << ", ";
		}
	}

	os << ")";

	return os;
}


Vec multPolynomials(const Vec& coeffs1, const Vec& coeffs2) {
	Vec coeffs(coeffs1.size() + coeffs2.size() - 1, 0.0);

	for (size_t i = 0; i < coeffs1.size(); ++i) {
		for (size_t j = 0; j < coeffs2.size(); ++j) {
			coeffs[i + j] += coeffs1[i] * coeffs2[j];
		}
	}

	return coeffs;
}

double evaluatePolynomial(const Vec& coeffs, double x) {
	double v = coeffs.back();

	for (size_t i = coeffs.size() - 1; i > 0; --i) {
		v = v * x + coeffs[i - 1];
	}

	return v;
}

double evaluatePolynomialDerivative(const Vec& coeffs, double x) {
	double v = coeffs.back() * (coeffs.size() - 1.0);

	for (size_t i = coeffs.size() - 2; i > 0; --i) {
		v = v * x + coeffs[i] * i;
	}

	return v;
}

void plotPolynomial(const Vec& coeffs, Graph* graph, double a = -1.0, double b = 1.0, int n = 20) {

	double step = (b - a) / (n - 1);

	Line l;
	for (int i = 0; i < n; ++i) {
		double x = a + step * i;
		l.addPoint(Point(x, evaluatePolynomial(coeffs, x)));
	}

	graph->addLine(l);
}



uint64_t factorial(int n) {
	uint64_t fac = 1;

	for (int i = 2; i <= n; ++i) {
		fac *= i;
	}

	return fac;
}


// special function that performs Newton's method on polynomials
double NewtonPolynomial(const Vec& coeffs, double x, int maxIter = 10, double eps = 1e-10) {

	size_t n = coeffs.size();

	for (int iter = 1; iter <= maxIter; ++iter) {
		double b = coeffs.back();
		double c = b;

		for (int i = n - 2; i > 0; --i) {
			b = b * x + coeffs[i];
			c = c * x + b;
		}

		b = b * x + coeffs[0];

		x -= b / c;

		if (b <= eps) break;
	}

	return x;
}


uint64_t binomialCoefficient(int a, int b) {

//	if (a < b) return 0; // this should never happen

	b = b > a - b ? a - b : b;

	if (b == 0) return 1;

	uint64_t res = a - b + 1;
	for (int i = 1; i < b; ++i) {
		res = res * (a - b + 1 + i) / (i + 1);
	}

	return res;
}


#endif
