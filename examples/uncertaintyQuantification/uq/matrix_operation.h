// matrix_operation.h
#ifndef MATRIX_OPERATION_H
#define MATRIX_OPERATION_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace MatrixOperations {

    // Function declarations
    double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must be the same size for dot product.");
        }
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            result += a[i] * b[i];
        }
        return result;
    }

    // Function to subtract one vector from another
    std::vector<double> vectorSubtraction(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must be the same size for subtraction.");
        }
        std::vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    // Function to multiply a vector by a scalar
    std::vector<double> vectorScalarProduct(const std::vector<double>& v, double scalar) {
        std::vector<double> result(v.size(), 0.0);
        for (size_t i = 0; i < v.size(); ++i) {
            result[i] = v[i] * scalar;
        }
        return result;
    }

    // Function to normalize a vector
    void normalize(std::vector<double>& v) {
        double norm = std::sqrt(dotProduct(v, v));
        if (norm == 0.0) {
            throw std::runtime_error("Cannot normalize a zero vector.");
        }
        for (size_t i = 0; i < v.size(); ++i) {
            v[i] /= norm;
        }
    }

    // Function to get the column of a matrix
    std::vector<double> getColumn(const std::vector<std::vector<double>>& matrix, size_t j) {
        if (matrix.empty() || j >= matrix[0].size()) {
            throw std::out_of_range("Invalid column index.");
        }
        std::vector<double> column(matrix.size());
        for (size_t i = 0; i < matrix.size(); ++i) {
            column[i] = matrix[i][j];
        }
        return column;
    }

    std::vector<std::vector<double>> generateIdentityMatrix(int n) {
        if (n <= 0) {
            throw std::invalid_argument("Size of identity matrix must be positive.");
        }
        std::vector<std::vector<double>> identityMatrix(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            identityMatrix[i][i] = 1.0;
        }
        return identityMatrix;
    }

    // Function to add one vector to another
    std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must be the same size for addition.");
        }
        std::vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    // Function to multiply two matrices
    std::vector<std::vector<double>> matrixMultiplication(const std::vector<std::vector<double>>& A,
                                                        const std::vector<std::vector<double>>& B) {
        if (A.empty() || B.empty()) {
            throw std::invalid_argument("Input matrices cannot be empty.");
        }
        size_t m = A.size();
        size_t n = A[0].size();
        size_t p = B[0].size();

        // Check dimensions
        if (n != B.size()) {
            throw std::invalid_argument("Number of columns in A must equal number of rows in B.");
        }

        // Initialize result matrix
        std::vector<std::vector<double>> C(m, std::vector<double>(p, 0.0));

        // Perform multiplication
        for (size_t i = 0; i < m; ++i) {
            if (A[i].size() != n) {
                throw std::invalid_argument("All rows in A must have the same number of columns.");
            }
            for (size_t j = 0; j < p; ++j) {
                for (size_t k = 0; k < n; ++k) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }

    // Function to transpose a matrix
    std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>>& A) {
        if (A.empty()) {
            return {};
        }
        size_t m = A.size();
        size_t n = A[0].size();
        std::vector<std::vector<double>> result(n, std::vector<double>(m));
        for (size_t i = 0; i < m; ++i) {
            if (A[i].size() != n) {
                throw std::invalid_argument("All rows in A must have the same number of columns.");
            }
            for (size_t j = 0; j < n; ++j) {
                result[j][i] = A[i][j];
            }
        }
        return result;
    }

    // Function to compute the Frobenius norm of a matrix
    double frobeniusNorm(const std::vector<std::vector<double>>& matrix) {
        double sum = 0.0;
        for (const auto& row : matrix) {
            for (double val : row) {
                sum += val * val;
            }
        }
        return std::sqrt(sum);
    }

    // Function to check if matrix A is converging to matrix B
    bool isConverging(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B,
                    double& former_norm, double tolerance) {
        // Ensure A and B are the same size
        if (A.size() != B.size() || A.empty() || B.empty() || A[0].size() != B[0].size()) {
            throw std::invalid_argument("Matrices must be the same size and non-empty.");
        }
        std::vector<std::vector<double>> diff(A.size(), std::vector<double>(A[0].size()));
        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = 0; j < A[0].size(); ++j) {
                diff[i][j] = A[i][j] - B[i][j];
            }
        }
        double norm = frobeniusNorm(diff);

        // Handle the first iteration
        if (former_norm == 0.0) {
            former_norm = norm;
            return false;
        }

        double relative_change = std::fabs(norm - former_norm) / former_norm;

        bool converged = relative_change <= tolerance;
        former_norm = norm;
        return converged;
    }

} // namespace MatrixOperations
#endif // MATRIX_OPERATION_H
