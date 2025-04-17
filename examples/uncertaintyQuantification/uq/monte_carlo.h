#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <vector>
#include <memory>
#include <random>
#include <stdexcept>
#include <cmath>
#include <functional>
#include "distribution.h" // Include Distribution and DistributionType

class MonteCarlo {
public:
    // Constructor for univariate distributions per dimension
    MonteCarlo(size_t numSamples, size_t randomNumberDimension, const std::vector<Distribution>& distributions, unsigned int seed)
        : numSamples(numSamples), randomNumberDimension(randomNumberDimension), distributions(distributions), rng(seed) {
        if (distributions.size() != randomNumberDimension && distributions.size() != 1) {
            throw std::invalid_argument("Number of distributions must match the random number dimension or be 1 for multivariate distribution.");
        }
    }

    // Constructor for single distribution applied to all dimensions
    MonteCarlo(size_t numSamples, size_t randomNumberDimension, Distribution distribution, unsigned int seed)
        : numSamples(numSamples), randomNumberDimension(randomNumberDimension), rng(seed) {
        distributions = std::vector<Distribution>(randomNumberDimension, distribution);
    }

    // Generate samples
    void generateSamples(std::vector<std::vector<double>>& samples) {
        samples.resize(numSamples, std::vector<double>(randomNumberDimension));

        // if (distributions.size() == 1 && distributions[0].type == DistributionType::MultivariateNormal) {
        //     // Multivariate Normal Distribution
        //     const Distribution& dist = distributions[0];

        //     // Validate the mean vector and covariance matrix dimensions
        //     if (dist.mean.size() != randomNumberDimension || dist.covariance.size() != randomNumberDimension) {
        //         throw std::runtime_error("Mean vector and covariance matrix dimensions must match the random number dimension.");
        //     }
        //     for (const auto& row : dist.covariance) {
        //         if (row.size() != randomNumberDimension) {
        //             throw std::runtime_error("Covariance matrix must be square and match the random number dimension.");
        //         }
        //     }

        //     // Perform Cholesky decomposition on the covariance matrix
        //     std::vector<std::vector<double>> L;
        //     if (!choleskyDecomposition(dist.covariance, L)) {
        //         throw std::runtime_error("Covariance matrix is not positive definite.");
        //     }

        //     // Generate samples
        //     std::normal_distribution<double> standard_normal(0.0, 1.0);
        //     for (size_t i = 0; i < numSamples; ++i) {
        //         // Generate standard normal random vector z
        //         std::vector<double> z(randomNumberDimension);
        //         for (size_t j = 0; j < randomNumberDimension; ++j) {
        //             z[j] = standard_normal(rng);
        //         }
        //         // Compute L * z
        //         std::vector<double> sample = multiplyMatrixVector(L, z);
        //         // Add the mean vector
        //         for (size_t j = 0; j < randomNumberDimension; ++j) {
        //             samples[i][j] = sample[j] + dist.mean[j];
        //         }
        //     }
        // } else {
            // Univariate distributions per dimension
            if (distributions.size() != randomNumberDimension) {
                throw std::runtime_error("Number of distributions must match the random number dimension.");
            }

            // Prepare distributions per dimension
            std::vector<std::function<double()>> randomGenerators(randomNumberDimension);

            for (size_t j = 0; j < randomNumberDimension; ++j) {
                const Distribution& dist = distributions[j];

                if (dist.type == DistributionType::Uniform) {
                    std::uniform_real_distribution<double> uniform_dist(dist.param1, dist.param2);
                    randomGenerators[j] = [this, uniform_dist]() mutable {
                        return uniform_dist(rng);
                    };
                } else if (dist.type == DistributionType::Normal) {
                    std::normal_distribution<double> normal_dist(dist.param1, dist.param2);
                    randomGenerators[j] = [this, normal_dist]() mutable {
                        return normal_dist(rng);
                    };
                } else {
                    throw std::runtime_error("Unsupported distribution type for univariate distribution.");
                }
            }

            // Generate samples
            for (size_t i = 0; i < numSamples; ++i) {
                for (size_t j = 0; j < randomNumberDimension; ++j) {
                    samples[i][j] = randomGenerators[j]();
                }
            }
        // }
    }

    size_t getSamplesNumber() const {
        return numSamples;
    }

private:
    size_t numSamples;
    size_t randomNumberDimension;
    std::vector<Distribution> distributions;
    std::mt19937 rng;  // Random number generator

    // Helper functions
    // bool choleskyDecomposition(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L) {
    //     size_t n = A.size();
    //     L.resize(n, std::vector<double>(n, 0.0));

    //     for (size_t i = 0; i < n; ++i) {
    //         for (size_t j = 0; j <= i; ++j) {
    //             double sum = 0.0;

    //             if (j == i) { // Diagonal elements
    //                 for (size_t k = 0; k < j; ++k)
    //                     sum += L[j][k] * L[j][k];
    //                 double val = A[j][j] - sum;
    //                 if (val <= 0.0)
    //                     return false; // Not positive definite
    //                 L[j][j] = std::sqrt(val);
    //             } else {
    //                 for (size_t k = 0; k < j; ++k)
    //                     sum += L[i][k] * L[j][k];
    //                 L[i][j] = (A[i][j] - sum) / L[j][j];
    //             }
    //         }
    //     }
    //     return true;
    // }

    // Helper function: Matrix-vector multiplication
    // std::vector<double> multiplyMatrixVector(const std::vector<std::vector<double>>& M, const std::vector<double>& v) {
    //     size_t n = M.size();
    //     std::vector<double> result(n, 0.0);

    //     for (size_t i = 0; i < n; ++i) {
    //         double sum = 0.0;
    //         for (size_t j = 0; j < v.size(); ++j) {
    //             sum += M[i][j] * v[j];
    //         }
    //         result[i] = sum;
    //     }
    //     return result;
    // }
};

#endif // MONTE_CARLO_H
