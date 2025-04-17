// generalized_polynomial_chaos.h
#ifndef GENERALIZED_POLYNOMIAL_CHAOS_H
#define GENERALIZED_POLYNOMIAL_CHAOS_H

#include <vector>
#include <memory>
#include <cmath>
#include <string>

#include "utils.h"

#include "distribution.h"

// Include the polynomial basis and quadrature headers
#include "polynomial.h"

// #include "quadrature.h"

class GeneralizedPolynomialChaos {
public:
    // Constructor
    GeneralizedPolynomialChaos(size_t order,
                               size_t nq,
                               const std::vector<Distribution>& distributions,
                               Quadrature::QuadratureMethod quadratureMethod);

    // Evaluation functions
    double evaluate(size_t n_order, size_t k);
    double evaluate(size_t n_order, size_t k, size_t phi_i);
    double evaluate(size_t n_order, const std::vector<size_t>& idx);
    double evaluate(size_t n_order, double x, size_t phi_i);
    double evaluate_polynomial(size_t order_max, size_t k);

    // Compute phiRan matrix
    void evaluatePhiRan();

    // Compute tensors
    void computeTensors();

    // Transformation functions
    void chaosToRandom(const std::vector<double>& chaosCoefficients, std::vector<double>& randomVariables);
    void randomToChaos(const std::vector<double>& randomVariables, std::vector<double>& chaosCoefficients);

    // Chaos operations
    void chaosProduct(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>& product);
    void chaosSum(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>& sum);

    // Statistical moments
    double mean(const std::vector<double>& chaosCoefficients);
    double std(const std::vector<double>& chaosCoefficients);

    void convert2affinePCE(const Distribution& distribution, std::vector<double>&chaos);

    // Getters
    size_t getPolynomialsOrder() const;
    size_t getQuadraturePointsNumber() const;
    void getPointsAndWeights(std::vector<std::vector<double>>& points, std::vector<std::vector<double>>& weights);
    void getStochasticCollocationSample(std::vector<std::vector<double>>& points);
    void getTensors(std::vector<double>& t2Product, std::vector<double>& t2Product_inv, std::vector<double>& t3Product);
    std::vector<double> getWeightsMultiplied() const;
    // Template function to get the polynomial basis at a specific dimension (i)
    template <typename PolynomialBasis>
    std::shared_ptr<PolynomialBasis> getPolynomialBasis(size_t i) const;
    std::vector<std::vector<size_t>> getMultiIndices() const;

    void getPhiRan(std::vector<double>& phiRan);
    void getCoefficients(std::vector<std::vector<std::vector<double>>>& polynomialCoeffs);


    // void get_polynomial_coefficients(std::vector<std::vector<double>>& polynomialCoeffs) {
    //     polynomialCoeffs = this->polynomialCoeffs;
    // }


private:
    size_t pointsWeightsMethod;
    size_t No; // Number of polynomials
    size_t nq; // Number of quadrature points per dimension
    size_t totalNq; // Total number of quadrature points
    size_t order;
    size_t randomNumberDimension;
    std::vector<std::vector<size_t>> inds; // Multi-indices
    std::vector<std::vector<double>> points;  // Points for each dimension
    std::vector<std::vector<double>> weights; // Weights for each dimension
    std::vector<std::vector<double>> pointsTensor; // Tensor product of points
    std::vector<double> weightsMultiplied;    // Combined weights
    std::vector<std::vector<size_t>> pointsWeightsIndexList;
    std::vector<std::vector<std::vector<double>>> coefficients; // Coefficients of polynomials

    Quadrature::QuadratureMethod quadratureMethod;

    std::vector<double> phiRan;   // Evaluated polynomials at quadrature points
    std::vector<double> phiRan_T; // Transpose of phiRan
    std::vector<double> t2Product;
    std::vector<double> t2Product_inv;
    std::vector<double> t3Product;

    // Distributions for each dimension
    std::vector<Distribution> distributions;

    // Polynomial bases for each dimension
    std::vector<std::shared_ptr<Polynomials::PolynomialBasis>> polynomialBases;

    // Initialization functions
    void initializeQuadratures();
    void initializeMatrices();

    void initializePolynomialCoefficients();

    // Helper functions
    std::vector<size_t> findIndex(size_t idx, size_t dimension, size_t nq);
    void calculateMultiIndices(size_t d, size_t n, std::vector<std::vector<size_t>>& indices);

    std::shared_ptr<Polynomials::PolynomialBasis> createPolynomialBasis(const Distribution& dist) {
        switch (dist.type) {
            case DistributionType::Uniform:
                return std::make_shared<Polynomials::LegendreBasis>();
            case DistributionType::Normal:
                return std::make_shared<Polynomials::HermiteBasis>();
            // Add cases for other distributions
            default:
                throw std::runtime_error("Unsupported distribution type for GPC.");
        }
    }

    // File I/O for tensor storage
    // bool fileExists(const std::string& filename);
    // void createDirectory(const std::string& directory);
    // void readVector1D(const std::string& filename, std::vector<double>& data);
    // void saveVector1D(const std::string& filename, const std::vector<double>& data);
};

#endif // GENERALIZED_POLYNOMIAL_CHAOS_H
