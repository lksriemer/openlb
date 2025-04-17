// quadrature.h
#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "utils.h"
#include "quadrature_base.h"
#include "polynomial.h"
#include "matrix_operation.h"

// Uncomment the following line if you have GSL installed and want to use it
// #define USE_GSL

#ifdef USE_GSL
#include <gsl/gsl_integration.h>
#endif

namespace Quadrature {

// Add this function to matrix_operation.h
double Pythag(double a, double b) {
    double absa = std::fabs(a);
    double absb = std::fabs(b);
    return (absa > absb ? absa * std::sqrt(1.0 + (absb / absa) * (absb / absa)) :
            (absb == 0.0 ? 0.0 : absb * std::sqrt(1.0 + (absa / absb) * (absa / absb))));
}

template <typename PolynomialBasis>
class Quadrature : public QuadratureBase {
public:
    Quadrature(size_t nq, QuadratureMethod method = QuadratureMethod::WilkinsonShiftQR)
        : nq(nq), basis(std::make_shared<PolynomialBasis>()), method(method) {
        computeQuadrature();
    }

    const std::vector<double>& getPoints() const override {
        return points;
    }

    const std::vector<double>& getWeights() const override {
        return weights;
    }

private:
    size_t nq;
    std::shared_ptr<PolynomialBasis> basis;
    QuadratureMethod method;
    std::vector<double> points;
    std::vector<double> weights;

    void computeQuadrature() {
#ifdef USE_GSL
        if (method == QuadratureMethod::GSL) {
            // std::cout << "Using GSL quadrature method." << std::endl;
            computeQuadratureGSL();
        } else {
            performQRDecomposition();  // Use Householder or Wilkinson QR based on the method
        }
#else
        if (method == QuadratureMethod::GSL) {
            std::cerr << "Warning: GSL is not enabled. Falling back to QR decomposition." << std::endl;
        }
        performQRDecomposition();  // Use Householder or Wilkinson QR based on the method
#endif
    }

    void performQRDecomposition() {
        // Step 1: Construct the Jacobi matrix using the basis
        auto J = basis->constructJacobiMatrix(nq);

        // Step 2: Perform the appropriate QR decomposition
        if (method == QuadratureMethod::HouseholderQR) {
            // std::cout << "Using HouseholderQR method." << std::endl;
            points = HouseholderQRDecomposition(J); // Perform Householder QR and get eigenvalues
        } else if (method == QuadratureMethod::WilkinsonShiftQR) {
            std::cerr << "Using Wilkinson's Shift QR method." << std::endl;
            points = WilkinsonShiftQRDecomposition(J); // Perform Wilkinson Shift QR
        } else {
            std::cerr << "Warning: Unsupported method. Defaulting to Wilkinson's Shift QR." << std::endl;
            points = WilkinsonShiftQRDecomposition(J); // Default to Wilkinson Shift QR
        }

        // Step 3: Compute quadrature weights
        computeWeights();
    }

#ifdef USE_GSL
    void computeQuadratureGSL() {
        if constexpr (std::is_same_v<PolynomialBasis, Polynomials::LegendreBasis>) {
            computeQuadraturePointsWeightsLegendreGSL(nq, points, weights);

            // Normalize weights if necessary
            double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
            for (auto& w : weights) {
                w /= sum;
            }
        } else if constexpr (std::is_same_v<PolynomialBasis, Polynomials::HermiteBasis>) {
            computeQuadraturePointsWeightsHermiteGSL(nq, points, weights);
            // No need to normalize weights for Gauss-Hermite quadrature adapted to normal distribution
        } else {
            throw std::runtime_error("GSL quadrature method is only supported for Legendre polynomials.");
        }
    }

    void computeQuadraturePointsWeightsLegendreGSL(size_t n, std::vector<double>& points, std::vector<double>& weights) {
        if (n <= 0) {
            throw std::invalid_argument("Number of quadrature points must be positive.");
        }

        points.resize(n);
        weights.resize(n);

        gsl_integration_glfixed_table* table = gsl_integration_glfixed_table_alloc(n);

        double xi, wi;
        for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
            gsl_integration_glfixed_point(-1.0, 1.0, i, &xi, &wi, table);
            points[i] = xi;
            weights[i] = wi * 0.5; // Adjust weights if necessary
        }

        gsl_integration_glfixed_table_free(table);
    }

    void computeQuadraturePointsWeightsHermiteGSL(size_t n, std::vector<double>& points, std::vector<double>& weights) {
        // TODO: Implement Gauss-Hermite quadrature using GSL
    }

#endif

    // Wilkinson Shift QR decomposition to compute the eigenvalues of a tridiagonal matrix
    std::vector<double> WilkinsonShiftQRDecomposition(const std::vector<std::vector<double>>& J) {
        size_t n = J.size();
        std::vector<double> d(n, 0.0); // Diagonal elements (eigenvalues)
        std::vector<double> e(n, 0.0); // Off-diagonal elements

        // Initialize d and e from the input matrix J
        for (size_t i = 0; i < n; ++i) {
            d[i] = J[i][i];
            if (i > 0) e[i] = J[i][i - 1];
        }

        const double eps = std::numeric_limits<double>::epsilon();
        double g, r, p, s, c, f, b;
        int iter, l, m, i;

        // Reduce subdiagonal elements
        for (i = 1; i < static_cast<int>(n); ++i) e[i - 1] = e[i];
        e[n - 1] = 0.0;

        for (l = 0; l < static_cast<int>(n); ++l) {
            iter = 0;
            do {
                // Find small subdiagonal element to isolate a diagonal block
                for (m = l; m < static_cast<int>(n) - 1; ++m) {
                    double dd = std::fabs(d[m]) + std::fabs(d[m + 1]);
                    if (std::fabs(e[m]) <= eps * dd) break;
                }

                // If no convergence yet, apply QR with Wilkinson shift
                if (m != l) {
                    if (++iter == 30) {
                        std::cerr << "Too many iterations in QR decomposition!\n";
                        return d;  // Early return if no convergence
                    }

                    // Wilkinson's shift calculation
                    g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                    r = Pythag(g, 1.0);
                    g = d[m] - d[l] + e[l] / (g + std::copysign(r, g));

                    s = c = 1.0;
                    p = 0.0;
                    for (i = m - 1; i >= l; --i) {
                        f = s * e[i];
                        b = c * e[i];
                        e[i + 1] = (r = Pythag(f, g));
                        if (r == 0.0) {
                            d[i + 1] -= p;
                            e[m] = 0.0;
                            break;
                        }
                        s = f / r;
                        c = g / r;
                        g = d[i + 1] - p;
                        r = (d[i] - g) * s + 2.0 * c * b;
                        d[i + 1] = g + (p = s * r);
                        g = c * r - b;
                    }

                    if (r == 0.0 && i >= l) continue;
                    d[l] -= p;
                    e[l] = g;
                    e[m] = 0.0;
                }
            } while (m != l);
        }

        // Return the diagonal elements as the eigenvalues
        std::sort(d.begin(), d.end());
        return d;
    }

        // Function to perform Householder QR decomposition
    void householderQR(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& Q, std::vector<std::vector<double>>& R) {
        size_t m = A.size(); // Number of rows
        size_t n = (A.empty() ? 0 : A[0].size()); // Number of columns (assuming all rows have the same number of columns)
        R = A;
        Q = ::MatrixOperations::generateIdentityMatrix(m);
        for (size_t k = 0; k < n; ++k) {
            std::vector<double> x(m - k, 0.0);
            for (size_t i = k; i < m; ++i) {
                x[i - k] = R[i][k];
            }

            double sign = 1.0;
            if(x[0] < 0) {
                sign = -1.0;
            }
            double norm = std::sqrt(::MatrixOperations::dotProduct(x, x));

            std::vector<double> v = ::MatrixOperations::vectorScalarProduct(x, 1.0 / (x[0] + sign * norm));
            v[0] = 1;
            double tau = 2.0 / ::MatrixOperations::dotProduct(v, v);

            std::vector<std::vector<double>> H = ::MatrixOperations::generateIdentityMatrix(m);

            for (size_t i = k; i < m; ++i) {
                for (size_t j = k; j < m; ++j) {
                    H[i][j] -= tau * v[i - k] * v[j - k];
                }
            }

            // Update R with Householder transformation
            R = ::MatrixOperations::matrixMultiplication(H, R);

            // Update Q with Householder transformation
            Q = ::MatrixOperations::matrixMultiplication(Q, ::MatrixOperations::transposeMatrix(H));
        }

    }

    // Householder QR decomposition to compute the eigenvalues of a matrix
    std::vector<double> HouseholderQRDecomposition(const std::vector<std::vector<double>>& J) {
        size_t n = J.size();
        std::vector<std::vector<double>> Q(n, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> J_new = J;
        std::vector<std::vector<double>> J_old = J;

        size_t iter = 1;
        size_t max_iter = 1000;
        double tol = 1e-12;

        while (iter < max_iter) {
            householderQR(J_new, Q, R);
            J_new = ::MatrixOperations::matrixMultiplication(R, Q);

            if (iter % 100 == 0) {
                if (isConverged(J_old, J_new, tol)) {
                    // std::cout << "Householder QR converged in " << iter << " iterations. Tol is " << tol << std::endl;
                    break;
                }
            }

            J_old = J_new;
            iter++;
        }

        std::vector<double> eigenvalues(n);
        for (size_t i = 0; i < n; ++i) {
            eigenvalues[i] = J_old[i][i]; // Extract diagonal as eigenvalues
        }
        std::sort(eigenvalues.begin(), eigenvalues.end());
        return eigenvalues;
    }

    bool isConverged(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B, double tol) {
        size_t n = A.size();
        double norm_diff = 0.0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                norm_diff += std::abs(A[i][j] - B[i][j]);
            }
        }
        // std::cout << "Norm diff: " << norm_diff << std::endl;
        return norm_diff < tol;
    }

    void computeWeights() {
        weights.resize(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
            double x = points[i];
            weights[i] = computeQuadratureWeight(*basis, x);
        }

        // Normalize weights if necessary
        double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        for (double& w : weights) {
            w /= sum;
        }
    }

    // Templated computeQuadratureWeight function
    template <typename Basis>
    double computeQuadratureWeight(const Basis& basis, double x) {
        static_assert(!std::is_same_v<Basis, Basis>, "Unsupported polynomial basis for quadrature weights.");
        return 0.0;
    }

    // Specialization for LegendreBasis
    double computeQuadratureWeight(const Polynomials::LegendreBasis& basis, double x) {
        double Pn_prime = basis.derivativePolynomial(nq, x);
        return 2.0 / ((1.0 - x * x) * Pn_prime * Pn_prime);
    }

    // Specialization for HermiteBasis
    // potential bug here
    double computeQuadratureWeight(const Polynomials::HermiteBasis& basis, double x) {
        double Hn_minus1 = basis.evaluatePolynomial(nq - 1, x);
        return std::pow(2, nq-1) * std::tgamma(nq+1) * std::sqrt(M_PI) / pow(nq * Hn_minus1, 2);
    }
};

} // namespace Quadrature

#endif // QUADRATURE_H
