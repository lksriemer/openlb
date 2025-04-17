#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>

// Enumeration to specify the distribution type
enum class DistributionType {
    Uniform,
    Normal,
    // MultivariateNormal
};

// Struct to hold distribution information for each dimension
struct Distribution {
    DistributionType type;

    // Parameters for univariate distributions
    double param1; // For Uniform: lower bound, for Normal: mean
    double param2; // For Uniform: upper bound, for Normal: standard deviation

    // Parameters for multivariate distributions
    std::vector<double> mean;                       // Mean vector
    std::vector<std::vector<double>> covariance;    // Covariance matrix

    // Constructor for univariate distributions
    Distribution(DistributionType type, double param1 = 0.0, double param2 = 1.0)
        : type(type), param1(param1), param2(param2) {}

    // Constructor for multivariate normal distribution
    // Distribution(const std::vector<double>& mean, const std::vector<std::vector<double>>& covariance)
    //     : type(DistributionType::MultivariateNormal), param1(0.0), param2(0.0), mean(mean), covariance(covariance) {}
};

double affine(double x, const Distribution& dist) {
    switch (dist.type) {
        case DistributionType::Uniform: {
            double a = dist.param1; // lower bound
            double b = dist.param2; // upper bound
            // Affine transformation from [-1, 1] to [a, b]
            return 0.5 * (b - a) * x + 0.5 * (a + b);
        }
        case DistributionType::Normal: {
            double mean = dist.param1;
            double stddev = dist.param2;
            // Affine transformation from standard normal to N(mean, stddev^2)
            // Assuming x is from standard normal distribution
            return mean + stddev * x;
        }
        // case DistributionType::MultivariateNormal:
        //     // Implement transformation for multivariate normal distribution
        //     break;
        default:
            throw std::runtime_error("Unsupported distribution type for affine transformation.");
    }
}

#endif // DISTRIBUTION_H
