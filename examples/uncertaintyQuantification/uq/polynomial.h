// polynomial.h
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "utils.h"
#include "quadrature_base.h"

namespace Polynomials {

// Abstract base class for all polynomial basis types
class PolynomialBasis {
public:
    virtual ~PolynomialBasis() = default;

    // Compute polynomial coefficients of order n
    virtual std::vector<double> computeCoefficients(size_t n) const = 0;

    // Construct the Jacobi matrix of size n
    virtual std::vector<std::vector<double>> constructJacobiMatrix(size_t n) const = 0;

    // Evaluate the polynomial of order n at point x
    double evaluatePolynomial(size_t n, double x) const;

    // Compute the derivative of the polynomial of order n at point x
    double derivativePolynomial(size_t n, double x) const;

    // Dynamically create and return a quadrature object
    virtual std::shared_ptr<Quadrature::QuadratureBase> getQuadrature(size_t nq, Quadrature::QuadratureMethod method) const = 0;

};

// Implement evaluatePolynomial using Horner's method for efficiency
inline double PolynomialBasis::evaluatePolynomial(size_t n, double x) const {
    std::vector<double> coeffs = computeCoefficients(n);

    // Evaluate polynomial using Horner's method
    double result = coeffs.back();
    for (size_t i = coeffs.size() - 1; i-- > 0;) {
        result = result * x + coeffs[i];
    }
    return result;
}

// Implement derivativePolynomial using Horner's method
inline double PolynomialBasis::derivativePolynomial(size_t n, double x) const {
    std::vector<double> coeffs = computeCoefficients(n);
    double result = 0.0;
    for (size_t i = coeffs.size() - 1; i > 0; --i) {
        result = result * x + i * coeffs[i];
    }
    return result;
}

// LegendreBasis class that inherits from PolynomialBasis
class LegendreBasis : public PolynomialBasis {
public:
    // Implement all pure virtual functions
    std::vector<double> computeCoefficients(size_t n) const override;
    std::vector<std::vector<double>> constructJacobiMatrix(size_t n) const override;
    std::shared_ptr<Quadrature::QuadratureBase> getQuadrature(size_t nq, Quadrature::QuadratureMethod method) const override;
};

// HermiteBasis class that inherits from PolynomialBasis
class HermiteBasis : public PolynomialBasis {
public:
    // Implement all pure virtual functions
    std::vector<double> computeCoefficients(size_t n) const override;
    std::vector<std::vector<double>> constructJacobiMatrix(size_t n) const override;
    std::shared_ptr<Quadrature::QuadratureBase> getQuadrature(size_t nq, Quadrature::QuadratureMethod method) const override;
};

} // namespace Polynomials

#endif // POLYNOMIAL_H
