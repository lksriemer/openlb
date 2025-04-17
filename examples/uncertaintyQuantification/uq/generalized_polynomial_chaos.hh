// generalized_polynomial_chaos.hh
#ifndef GENERALIZED_POLYNOMIAL_CHAOS_HH
#define GENERALIZED_POLYNOMIAL_CHAOS_HH

#include "generalized_polynomial_chaos.h"

#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <functional>

// Include the polynomial basis and quadrature headers
// #include "quadrature.h"

// Namespace aliases
using LegendreBasis = Polynomials::LegendreBasis;
using HermiteBasis = Polynomials::HermiteBasis;
// using Quadrature = Quadrature::Quadrature;


// Constructor
GeneralizedPolynomialChaos::GeneralizedPolynomialChaos(size_t _order,
                                                       size_t _nq,
                                                       const std::vector<Distribution>& _distributions,
                                                       Quadrature::QuadratureMethod _quadratureMethod)
    : pointsWeightsMethod(0),
      No(0),
      nq(_nq),
      totalNq(0),
      order(_order),
      randomNumberDimension(0),
      quadratureMethod(_quadratureMethod),
      distributions(_distributions) // Initialize distributions
{
    // Map distributions to polynomial bases
    polynomialBases.clear();
    for (const auto& dist : distributions) {
        polynomialBases.push_back(createPolynomialBasis(dist));
    }

    // Set randomNumberDimension to the size of distributions
    randomNumberDimension = distributions.size();

    // Calculate multi-indices
    calculateMultiIndices(randomNumberDimension, order, inds);
    No = inds.size();

    // Initialize polynomial coefficients, quadratures, and matrices
    // initializePolynomialCoefficients();
    initializeQuadratures();
    initializeMatrices();

    // Evaluate polynomials at quadrature points
    evaluatePhiRan();

    // Compute tensors
    computeTensors();
}


// Initialize quadratures
void GeneralizedPolynomialChaos::initializeQuadratures() {
    points.resize(randomNumberDimension);
    weights.resize(randomNumberDimension);

    totalNq = std::pow(nq, randomNumberDimension);

    // std::cout << "Initializing quadratures..." << std::endl;

    for (size_t i = 0; i < randomNumberDimension; ++i) {
        auto quadrature = polynomialBases[i]->getQuadrature(nq, quadratureMethod);
        points[i] = quadrature->getPoints();
        weights[i] = quadrature->getWeights();

        // Print points and weights
        // std::cout << "Dimension " << i << std::endl;
        // for (size_t j = 0; j < nq; ++j) {
        //     std::cout << points[i][j] << " " << weights[i][j] << std::endl;
        // }
    }
}



// Initialize matrices
void GeneralizedPolynomialChaos::initializeMatrices() {
    // std::cout << "Initializing matrices..." << std::endl;
    // Resize vectors
    phiRan.resize(totalNq * No, 0.0);
    phiRan_T.resize(totalNq * No, 0.0);
    t2Product.resize(No, 0.0);
    t2Product_inv.resize(No, 0.0);
    t3Product.resize(No * No * No, 0.0);


    // auto a = findIndex(1, randomNumberDimension, nq);
    // std::cout << "findIndex i = 1, nq = 5: " << a[0] << std::endl;

    // Generate pointsWeightsIndexList
    pointsWeightsIndexList.resize(totalNq, std::vector<size_t>(randomNumberDimension));
    for (size_t i = 0; i < totalNq; ++i) {
        pointsWeightsIndexList[i] = findIndex(i, randomNumberDimension, nq);
    }


    // Compute weightsMultiplied and pointsTensor
    weightsMultiplied.resize(totalNq, 1.0);
    pointsTensor.resize(totalNq, std::vector<double>(randomNumberDimension));
    for (size_t k = 0; k < totalNq; ++k) {
        pointsTensor[k].resize(randomNumberDimension);
        for (size_t dim = 0; dim < randomNumberDimension; ++dim) {
            size_t idx = pointsWeightsIndexList[k][dim];
            weightsMultiplied[k] *= weights[dim][idx];
            pointsTensor[k][dim] = points[dim][idx];
        }
    }
}

// Initialize polynomial coefficients
void GeneralizedPolynomialChaos::initializePolynomialCoefficients() {
    // std::cout << "Initializing polynomial coefficients..." << std::endl;
    coefficients.resize(randomNumberDimension);
    for (size_t phi_i = 0; phi_i < randomNumberDimension; ++phi_i) {
        auto basis = std::static_pointer_cast<LegendreBasis>(polynomialBases[phi_i]);

        coefficients[phi_i].resize(No);
        // std::cout << "Dimension " << phi_i << std::endl;
        for (size_t i = 0; i < No; ++i) {
            coefficients[phi_i][i] = basis->computeCoefficients(i);
            for (size_t j = 0; j <= order; ++j) {
                // std::cout << coefficients[phi_i][i][j] << " ";
            }
            // std::cout << std::endl;
        }
    }

}

// Evaluate n_order polynomial at point k
double GeneralizedPolynomialChaos::evaluate(size_t n_order, size_t k) {
    double result = 1.0;
    for (size_t i = 0; i < randomNumberDimension; ++i) {
        result *= evaluate(inds[n_order][i], k, i);
    }
    return result;
}

// Evaluate n_order polynomial at point k and dimension phi_i
double GeneralizedPolynomialChaos::evaluate(size_t n_order, size_t k, size_t phi_i) {
    double x = points[phi_i][pointsWeightsIndexList[k][phi_i]];
    return evaluate(n_order, x, phi_i);
}

// Evaluate n_order polynomial at given multi-index
double GeneralizedPolynomialChaos::evaluate(size_t n_order, const std::vector<size_t>& idx) {
    double result = 1.0;
    for (size_t i = 0; i < randomNumberDimension; ++i) {
        result *= evaluate(inds[n_order][i], points[i][idx[i]], i);
    }
    return result;
}

// Evaluate polynomial basis at given order, point x, and dimension phi_i
double GeneralizedPolynomialChaos::evaluate(size_t n_order, double x, size_t phi_i) {
    if (phi_i < 0 || phi_i >= polynomialBases.size()) {
        throw std::out_of_range("Invalid dimension index phi_i.");
    }
    return polynomialBases[phi_i]->evaluatePolynomial(n_order, x);
}

// Evaluate the polynomial at kth point up to order_max
double GeneralizedPolynomialChaos::evaluate_polynomial(size_t order_max, size_t k) {
    double sum = 0.0;
    for (size_t i = 0; i <= order_max; ++i) {
        sum += evaluate(i, k);
    }
    return sum;
}

// Evaluate phiRan matrix
void GeneralizedPolynomialChaos::evaluatePhiRan() {
    // std::cout << "Evaluating phiRan matrix..." << std::endl;
    for (size_t k = 0; k < totalNq; ++k) {
        for (size_t i = 0; i < No; ++i) {
            phiRan[k * No + i] = evaluate(i, pointsWeightsIndexList[k]);
            // std::cout << phiRan[k * No + i] << " ";
            phiRan_T[i * totalNq + k] = phiRan[k * No + i];
        }
        // std::cout << std::endl;
    }
}

// Helper functions
void GeneralizedPolynomialChaos::calculateMultiIndices(size_t d, size_t n, std::vector<std::vector<size_t>>& indices) {

    std::vector<size_t> index(d, 0);

    std::function<void(size_t, size_t, size_t)> recursiveFunction = [&](size_t pos, size_t sum, size_t maxOrder) {
        if (pos == d - 1) {
            index[pos] = maxOrder - sum;
            indices.push_back(index);
            return;
        }

        for (size_t i = 0; i <= maxOrder - sum; ++i) {
            index[pos] = i;
            recursiveFunction(pos + 1, sum + i, maxOrder);
        }
    };

    for (size_t order = 0; order <= n; ++order) {
        recursiveFunction(0, 0, order);
    }
}

std::vector<size_t> GeneralizedPolynomialChaos::findIndex(size_t idx, size_t dimension, size_t nq) {
    if (dimension == 1) {
        return {idx};
    }

    // General case for dimension > 1
    std::vector<size_t> index(dimension);
    for (size_t i = dimension; i-- > 0;) {  // Loop from dimension-1 to 0
        index[i] = idx % nq;
        idx /= nq;
    }
    return index;
}

// Compute tensors (t2Product and t3Product)
void GeneralizedPolynomialChaos::computeTensors() {
    // std::cout << "Computing tensors..." << std::endl;
    // File paths for saved matrices
    const std::string directoryT2Product = "./t2Product/";
    if (!directoryExists(directoryT2Product)) {
        createDirectory(directoryT2Product);
    }
    const std::string directoryT3Product = "./t3Product/";
    if (!directoryExists(directoryT3Product)) {
        createDirectory(directoryT3Product);
    }
    const std::string t2ProductFile = directoryT2Product + "dims_" + std::to_string(randomNumberDimension) + "_order_" + std::to_string(order) + "_nq_" + std::to_string(nq) + ".bin";
    const std::string t3ProductFile = directoryT3Product + "dims_" + std::to_string(randomNumberDimension) + "_order_" + std::to_string(order) + "_nq_" + std::to_string(nq) + ".bin";

    // Compute t2Product
    // std::cout << "Computing t2Product..." << std::endl;
    // if (fileExists(t2ProductFile)) {
        // std::cout << "Loading t2Product from file." << std::endl;
        // readVector1D(t2ProductFile, t2Product);
    // } else {
        // std::cout << "Calculating t2Product." << std::endl;

        // std::vector<std::vector<double>> tensor2d(No, std::vector<double>(No, 0.0));
        // for (size_t i = 0; i < No; ++i) {
        //     for (size_t j = 0; j < No; ++j) {
        //         for (size_t m = 0; m < totalNq; ++m) {
        //             tensor2d[i][j] +=  phiRan[m * No + i] * phiRan[m * No + j] * weightsMultiplied[m];
        //         }
        //     }
        // }
        // // std::cout << "t2Product" << std::endl;
        // for (size_t i = 0; i < No; ++i) {
        //     t2Product[i] = tensor2d[i][i];
        //     // std::cout << t2Product[i] << " ";
        // }
        // // std::cout << std::endl;

        for (size_t i = 0; i < No; ++i) {
            for (size_t m = 0; m < totalNq; ++m) {
                t2Product[i] += phiRan[m * No + i] * phiRan[m * No + i] * weightsMultiplied[m];
            }
        }

        // saveVector1D(t2ProductFile, t2Product);
    // }

    for (size_t i = 0; i < No; ++i) {
        t2Product_inv[i] = 1.0 / t2Product[i];
    }

    // Compute t3Product
    // std::cout << "Computing t3Product..." << std::endl;
    // if (fileExists(t3ProductFile)) {
        // std::cout << "Loading t3Product from file." << std::endl;
        // readVector1D(t3ProductFile, t3Product);
    // } else {
        // std::cout << "Calculating t3Product." << std::endl;
        for (size_t i = 0; i < No; ++i) {
            for (size_t j = 0; j < No; ++j) {
                for (size_t k = 0; k < No; ++k) {
                    double sum = 0.0;
                    for (size_t m = 0; m < totalNq; ++m) {
                        sum += phiRan[m * No + i] * phiRan[m * No + j] * phiRan[m * No + k] * weightsMultiplied[m];
                    }
                    t3Product[i * No * No + j * No + k] = sum;
                }
            }
        }
        // saveVector1D(t3ProductFile, t3Product);
    // }
    // std::cout << "Tensors computed." << std::endl;
}

// Transformation functions
void GeneralizedPolynomialChaos::chaosToRandom(const std::vector<double>& chaosCoefficients, std::vector<double>& randomVariables) {
    randomVariables.resize(totalNq, 0.0);

    for (size_t k = 0; k < totalNq; ++k) {
        auto startIt = phiRan.begin() + k * No;
        randomVariables[k] = std::inner_product(chaosCoefficients.begin(), chaosCoefficients.end(), startIt, 0.0);
    }
}

void GeneralizedPolynomialChaos::randomToChaos(const std::vector<double>& randomVariables, std::vector<double>& chaosCoefficients) {
    chaosCoefficients.resize(No, 0.0);
    std::vector<double> weightedRandomVariables(totalNq);

    // Compute weighted random variables
    for (size_t k = 0; k < totalNq; ++k) {
        weightedRandomVariables[k] = weightsMultiplied[k] * randomVariables[k];
    }

    // Compute chaos coefficients
    for (size_t i = 0; i < No; ++i) {
        auto startIt = phiRan_T.begin() + i * totalNq;
        chaosCoefficients[i] = std::inner_product(weightedRandomVariables.begin(), weightedRandomVariables.end(), startIt, 0.0);
        chaosCoefficients[i] *= t2Product_inv[i];
    }
}

// Chaos operations
void GeneralizedPolynomialChaos::chaosProduct(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>& product) {
    product.resize(No, 0.0);
    std::vector<double> precomputedProductsFlat(No * No);

    for (size_t j = 0; j < No; ++j) {
        for (size_t k = 0; k < No; ++k) {
            precomputedProductsFlat[j * No + k] = chaos1[j] * chaos2[k];
        }
    }

    for (size_t i = 0; i < No; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < No; ++j) {
            for (size_t k = 0; k < No; ++k) {
                size_t flatIndex = i * No * No + j * No + k;
                sum += precomputedProductsFlat[j * No + k] * t3Product[flatIndex];
            }
        }
        product[i] = sum * t2Product_inv[i];
    }
}

void GeneralizedPolynomialChaos::chaosSum(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>& sum) {
    sum.resize(No);
    for (size_t i = 0; i < No; ++i) {
        sum[i] = chaos1[i] + chaos2[i];
    }
}

// Statistical moments
double GeneralizedPolynomialChaos::mean(const std::vector<double>& chaosCoefficients) {
    return chaosCoefficients[0];
}

double GeneralizedPolynomialChaos::std(const std::vector<double>& chaosCoefficients) {
    double variance = 0.0;
    for (size_t i = 1; i < No; ++i) {
        variance += t2Product[i] * chaosCoefficients[i] * chaosCoefficients[i];
    }
    return std::sqrt(variance);
}

void GeneralizedPolynomialChaos::convert2affinePCE(const Distribution& distribution, std::vector<double>&chaos)
{
    switch (distribution.type) {
        case DistributionType::Uniform: {
            double a1 = 0.5 * (distribution.param1 + distribution.param2);
            double a2 = 0.5 * (distribution.param2 - distribution.param1);
            chaos[0] = a1;
            chaos[1] = a2;
        }
        case DistributionType::Normal: {
            chaos[0] = distribution.param1;
            chaos[1] = distribution.param2;
        }
        // Add cases for other distributions
        default:
            throw std::runtime_error("Unsupported distribution type for GPC.");
    }
}

// Getters
size_t GeneralizedPolynomialChaos::getPolynomialsOrder() const {
    return No;
}

size_t GeneralizedPolynomialChaos::getQuadraturePointsNumber() const {
    return totalNq;
}


void GeneralizedPolynomialChaos::getPointsAndWeights(std::vector<std::vector<double>>& points, std::vector<std::vector<double>>& weights) {
    points = this->points;
    weights = this->weights;
}

void GeneralizedPolynomialChaos::getStochasticCollocationSample(std::vector<std::vector<double>>& samples) {
    // Check that 'samples' has the correct size
    if (samples.size() != totalNq || samples[0].size() != randomNumberDimension) {
        samples.resize(totalNq, std::vector<double>(randomNumberDimension));
    }


    for (size_t j = 0; j < randomNumberDimension; ++j) {
        for (size_t i = 0; i < totalNq; ++i) {
            samples[i][j] = affine(points[j][pointsWeightsIndexList[i][j]], distributions[j]);
        }

    }

}

std::vector<double> GeneralizedPolynomialChaos::getWeightsMultiplied() const {
    return weightsMultiplied;
}

void GeneralizedPolynomialChaos::getTensors(std::vector<double>& t2Product, std::vector<double>& t2Product_inv, std::vector<double>& t3Product) {
    t2Product = this->t2Product;
    t2Product_inv = this->t2Product_inv;
    t3Product = this->t3Product;
}

template <typename PolynomialBasis>
std::shared_ptr<PolynomialBasis> GeneralizedPolynomialChaos::getPolynomialBasis(size_t dimension) const {
    if (dimension < 0 || dimension >= randomNumberDimension) {
        throw std::out_of_range("Dimension is out of bounds");
    }

    // Cast the void pointer back to the correct polynomial basis type
    return std::static_pointer_cast<PolynomialBasis>(polynomialBases[dimension]);
}

std::vector<std::vector<size_t>> GeneralizedPolynomialChaos::getMultiIndices() const {
    return inds;
}

void GeneralizedPolynomialChaos::getPhiRan(std::vector<double>& phiRan) {
    phiRan = this->phiRan;
}

void GeneralizedPolynomialChaos::getCoefficients(std::vector<std::vector<std::vector<double>>>& polynomialCoeffs) {
    polynomialCoeffs = this->coefficients;
}

#endif // GENERALIZED_POLYNOMIAL_CHAOS_HH