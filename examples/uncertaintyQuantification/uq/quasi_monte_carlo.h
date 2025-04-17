#ifndef QUASI_MONTE_CARLO_H
#define QUASI_MONTE_CARLO_H

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include "distribution.h" // For Distribution and DistributionType

class QuasiMonteCarlo {
public:
    // Constructor for multiple uniform distributions per dimension
    QuasiMonteCarlo(size_t numSamples, size_t randomNumberDimension,
                    const std::vector<Distribution>& distributions,
                    const std::string& dir_file = "new-joe-kuo-6.21201")
        : numSamples(numSamples),
          randomNumberDimension(randomNumberDimension),
          distributions(distributions),
          sobol_dir_file(dir_file) {
        if (distributions.size() != randomNumberDimension) {
            throw std::invalid_argument("Number of distributions must match the random number dimension.");
        }
        for (const auto& dist : distributions) {
            if (dist.type != DistributionType::Uniform) {
                throw std::invalid_argument("QuasiMonteCarlo now only supports Uniform distributions.");
            }
        }
    }

    // Constructor for a single uniform distribution applied to all dimensions
    QuasiMonteCarlo(size_t numSamples, size_t randomNumberDimension,
                    Distribution distribution,
                    const std::string& dir_file = "new-joe-kuo-6.21201")
        : numSamples(numSamples),
          randomNumberDimension(randomNumberDimension),
          sobol_dir_file(dir_file) {
        if (distribution.type != DistributionType::Uniform) {
            throw std::invalid_argument("QuasiMonteCarlo now only supports Uniform distributions.");
        }
        distributions = std::vector<Distribution>(randomNumberDimension, distribution);
    }

    // Generate samples
    void generateSamples(std::vector<std::vector<double>>& samples) {
        // Generate Sobol sequence points
        auto sobolPoints = sobol_points(static_cast<unsigned>(numSamples),
                                        static_cast<unsigned>(randomNumberDimension),
                                        sobol_dir_file);

        // Map Sobol points to desired uniform distributions
        samples.resize(numSamples, std::vector<double>(randomNumberDimension));

        for (size_t i = 0; i < numSamples; ++i) {
            for (size_t j = 0; j < randomNumberDimension; ++j) {
                double u = sobolPoints[i][j]; // Sobol point in [0,1]
                const Distribution& dist = distributions[j];
                // Map u to Uniform(param1, param2)
                samples[i][j] = dist.param1 + (dist.param2 - dist.param1) * u;
            }
        }
    }

    // Get the number of samples
    size_t getSamplesNumber() const {
        return numSamples;
    }

private:
    size_t numSamples;
    size_t randomNumberDimension;
    std::vector<Distribution> distributions;
    std::string sobol_dir_file;

    // Helper function to generate Sobol sequence points
    std::vector<std::vector<double>> sobol_points(unsigned N, unsigned D, const std::string& dir_file) {
        std::ifstream infile(dir_file);
        if (!infile) {
            throw std::runtime_error("Input file containing direction numbers cannot be found!");
        }

        // Read header line
        std::string buffer;
        std::getline(infile, buffer);

        // L = max number of bits needed
        unsigned L = static_cast<unsigned>(std::ceil(std::log(static_cast<double>(N)) / std::log(2.0)));

        // C[i] = index from the right of the first zero bit of i
        std::vector<unsigned> C(N);
        C[0] = 1;
        for (unsigned i = 1; i < N; i++) {
            C[i] = 1;
            unsigned value = i;
            while (value & 1) {
                value >>= 1;
                C[i]++;
            }
        }

        // POINTS[i][j] = the jth component of the ith point
        std::vector<std::vector<double>> POINTS(N, std::vector<double>(D));

        // Initialize direction numbers V
        std::vector<std::vector<unsigned>> V(D, std::vector<unsigned>(L + 1));

        // Read direction numbers for each dimension
        for (unsigned j = 0; j < D; ++j) {
            if (j == 0) {
                // First dimension
                for (unsigned i = 1; i <= L; ++i) {
                    V[j][i] = 1U << (32 - i); // All m's = 1
                }
            } else {
                // Read in parameters from file
                unsigned d, s, a;
                infile >> d >> s >> a;
                std::vector<unsigned> m(s + 1);
                for (unsigned i = 1; i <= s; ++i) {
                    infile >> m[i];
                }

                // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
                if (L <= s) {
                    for (unsigned i = 1; i <= L; ++i) {
                        V[j][i] = m[i] << (32 - i);
                    }
                } else {
                    for (unsigned i = 1; i <= s; ++i) {
                        V[j][i] = m[i] << (32 - i);
                    }
                    for (unsigned i = s + 1; i <= L; ++i) {
                        V[j][i] = V[j][i - s] ^ (V[j][i - s] >> s);
                        for (unsigned k = 1; k <= s - 1; ++k) {
                            V[j][i] ^= (((a >> (s - 1 - k)) & 1U) * V[j][i - k]);
                        }
                    }
                }
            }
        }

        // Evaluate X[0] to X[N-1], scaled by pow(2,32)
        std::vector<unsigned> X(D, 0);
        for (unsigned i = 0; i < N; ++i) {
            if (i > 0) {
                for (unsigned j = 0; j < D; ++j) {
                    X[j] ^= V[j][C[i - 1]];
                }
            }
            for (unsigned j = 0; j < D; ++j) {
                POINTS[i][j] = static_cast<double>(X[j]) / std::pow(2.0, 32);
            }
        }

        return POINTS;
    }
};

#endif // QUASI_MONTE_CARLO_H
