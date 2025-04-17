#ifndef UNCERTAINTY_QUANTIFICATION_HH
#define UNCERTAINTY_QUANTIFICATION_HH

#include "uncertainty_quantification.h"
#include <numeric>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

// Constructor
UncertaintyQuantification::UncertaintyQuantification(UQMethod uqMethod)
    : uqMethod(uqMethod),
      rng(std::random_device{}()) {
    // Initialize variables to default values
    order = 0;
    nq = 0;
    numSamples = 0;
    randomNumberDimension = 0;
    quadratureMethod = Quadrature::QuadratureMethod::HouseholderQR;
    ops = nullptr; // Initialize the unique_ptr to nullptr
    monteCarlo = nullptr;
    quasiMonteCarlo = nullptr;
    lhs = nullptr;
}

// Initialization function for GPC
void UncertaintyQuantification::initializeGPC(size_t order, size_t nq,
                                              const std::vector<Distribution>& distributions,
                                              Quadrature::QuadratureMethod quadratureMethod) {
    this->order = order;
    this->nq = nq;
    this->quadratureMethod = quadratureMethod;
    this->randomNumberDimension = distributions.size();
    this->distributions = distributions;


    // Initialize the GeneralizedPolynomialChaos object
    ops = std::make_unique<GeneralizedPolynomialChaos>(order, nq, distributions, quadratureMethod);
    // Get quadrature points and weights from ops
    ops->getPointsAndWeights(points, weights);

    No = ops->getPolynomialsOrder();

    // Compute the combined weights
    numSamples = ops->getQuadraturePointsNumber();
    weightsMultiplied = ops->getWeightsMultiplied();
    // Get multi-indices from ops
    multiIndices = ops->getMultiIndices();
}


// Initialization function for Monte Carlo with vector of distributions
void UncertaintyQuantification::initializeMonteCarlo(size_t numSamples, size_t randomNumberDimension, const std::vector<Distribution>& distributions, unsigned int seed) {
    this->numSamples = numSamples;
    this->randomNumberDimension = randomNumberDimension;
    this->distributions = distributions;

    // Create MonteCarlo instance
    monteCarlo = std::make_unique<MonteCarlo>(numSamples, randomNumberDimension, distributions, seed);
}

// Overload for a single distribution applied to all dimensions
void UncertaintyQuantification::initializeMonteCarlo(size_t numSamples, size_t randomNumberDimension, Distribution distribution, unsigned int seed) {
    this->numSamples = numSamples;
    this->randomNumberDimension = randomNumberDimension;
    this->distributions = distributions;

    // Create MonteCarlo instance
    monteCarlo = std::make_unique<MonteCarlo>(numSamples, randomNumberDimension, distribution, seed);
}

// Initialization function for Quasi-Monte Carlo
void UncertaintyQuantification::initializeQuasiMonteCarlo(size_t numSamples, size_t randomNumberDimension,
                                                          const std::vector<Distribution>& distributions,
                                                          const std::string& dir_file) {
    this->numSamples = numSamples;
    this->randomNumberDimension = randomNumberDimension;
    this->distributions = distributions;

    // Create QuasiMonteCarlo instance
    quasiMonteCarlo = std::make_unique<QuasiMonteCarlo>(numSamples, randomNumberDimension, distributions, dir_file);
}

// Overload for a single distribution applied to all dimensions
void UncertaintyQuantification::initializeQuasiMonteCarlo(size_t numSamples, size_t randomNumberDimension,
                                                          Distribution distribution,
                                                          const std::string& dir_file) {
    this->numSamples = numSamples;
    this->randomNumberDimension = randomNumberDimension;

    // Create QuasiMonteCarlo instance
    quasiMonteCarlo = std::make_unique<QuasiMonteCarlo>(numSamples, randomNumberDimension, distribution, dir_file);
}

// Initialization function for Latin Hypercube Sampling
void UncertaintyQuantification::initializeLatinHypercubeSampling(size_t numSamples, size_t randomNumberDimension) {
    this->numSamples = numSamples;
    this->randomNumberDimension = randomNumberDimension;

    // Create LatinHypercubeSampling instance
    lhs = std::make_unique<LatinHypercubeSampling>(static_cast<int>(numSamples), static_cast<int>(randomNumberDimension));
}

// Function to get sampling points
void UncertaintyQuantification::getSamplingPoints(std::vector<std::vector<double>>& out_points) {
    switch (uqMethod) {
        case UQMethod::GPC:
            if (!ops) {
                throw std::runtime_error("GPC has not been initialized. Call initializeGPC() first.");
            }
            ops->getStochasticCollocationSample(out_points);
            break;
        case UQMethod::MonteCarlo:
            if (!monteCarlo) {
                throw std::runtime_error("Monte Carlo has not been initialized. Call initializeMonteCarlo() first.");
            }
            monteCarlo->generateSamples(out_points);
            break;
        case UQMethod::QuasiMonteCarlo:
            if (!quasiMonteCarlo) {
                throw std::runtime_error("Quasi-Monte Carlo has not been initialized. Call initializeQuasiMonteCarlo() first.");
            }
            quasiMonteCarlo->generateSamples(out_points);
            break;
        case UQMethod::LatinHypercubeSampling:
            if (!lhs) {
                throw std::runtime_error("Latin Hypercube Sampling has not been initialized. Call initializeLatinHypercubeSampling() first.");
            }
            out_points = lhs->generateSamples();
            break;
        default:
            throw std::runtime_error("Invalid UQ method.");
    }
}

size_t UncertaintyQuantification::getSamplesNumber() {
    switch (uqMethod) {
        case UQMethod::GPC:
            if (!ops) {
                throw std::runtime_error("GPC has not been initialized. Call initializeGPC() first.");
            }
            return ops->getQuadraturePointsNumber();
            break;
        case UQMethod::MonteCarlo:
            if (!monteCarlo) {
                throw std::runtime_error("Monte Carlo has not been initialized. Call initializeMonteCarlo() first.");
            }
            return monteCarlo->getSamplesNumber();
            break;
        case UQMethod::QuasiMonteCarlo:
            if (!quasiMonteCarlo) {
                throw std::runtime_error("Quasi-Monte Carlo has not been initialized. Call initializeQuasiMonteCarlo() first.");
            }
            return quasiMonteCarlo->getSamplesNumber();
            break;
        case UQMethod::LatinHypercubeSampling:
            if (!lhs) {
                throw std::runtime_error("Latin Hypercube Sampling has not been initialized. Call initializeLatinHypercubeSampling() first.");
            }
            return lhs->getSamplesNumber();
            break;
        default:
            throw std::runtime_error("Invalid UQ method.");
    }
}

// Example function evaluation (to be customized)
double UncertaintyQuantification::mean(const std::vector<double>& input) {
    double mean = 0.0;
    switch (uqMethod) {
        case UQMethod::GPC: {
            if (!ops) {
                throw std::runtime_error("GPC has not been initialized. Call initializeGPC() first.");
            }
            std::vector<double> chaos(No, 0.0);
            ops->randomToChaos(input, chaos);
            mean = ops->mean(chaos);
            break;
        }
        case UQMethod::MonteCarlo: {
            if (!monteCarlo) {
                throw std::runtime_error("Monte Carlo has not been initialized. Call initializeMonteCarlo() first.");
            }
            double sum = std::accumulate(input.begin(), input.end(), 0.0);
            mean = sum / input.size();
            break;
        }
        case UQMethod::QuasiMonteCarlo: {
            if (!quasiMonteCarlo) {
                throw std::runtime_error("Quasi-Monte Carlo has not been initialized. Call initializeQuasiMonteCarlo() first.");
            }
            double sum = std::accumulate(input.begin(), input.end(), 0.0);
            mean = sum / input.size();
            break;
        }
        default:
            throw std::runtime_error("Invalid UQ method.");
    }
    return mean;
}

double UncertaintyQuantification::std(const std::vector<double>& input) {
    double std = 0.0;
    switch (uqMethod) {
        case UQMethod::GPC: {
            if (!ops) {
                throw std::runtime_error("GPC has not been initialized. Call initializeGPC() first.");
            }
            std::vector<double> chaos(No, 0.0);
            ops->randomToChaos(input, chaos);
            std = ops->std(chaos);
            break;
        }
        case UQMethod::MonteCarlo: {
            if (!monteCarlo) {
                throw std::runtime_error("Monte Carlo has not been initialized. Call initializeMonteCarlo() first.");
            }
            double avg = mean(input);
            double sumSq = std::accumulate(input.begin(), input.end(), 0.0, [avg](double acc, double val) {
                return acc + (val - avg) * (val - avg);
            });
            std = std::sqrt(sumSq / (input.size() - 1));
            break;
        }
        case UQMethod::QuasiMonteCarlo: {
            if (!quasiMonteCarlo) {
                throw std::runtime_error("Quasi-Monte Carlo has not been initialized. Call initializeQuasiMonteCarlo() first.");
            }
            double avg = mean(input);
            double sumSq = std::accumulate(input.begin(), input.end(), 0.0, [avg](double acc, double val) {
                return acc + (val - avg) * (val - avg);
            });
            std = std::sqrt(sumSq / (input.size() - 1));
            break;
        }
        default:
            throw std::runtime_error("Invalid UQ method.");
    }
    return std;
}

#endif // UNCERTAINTY_QUANTIFICATION_HH
