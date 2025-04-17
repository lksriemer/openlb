#ifndef LATIN_HYPERCUBE_SAMPLING_H
#define LATIN_HYPERCUBE_SAMPLING_H

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>

class LatinHypercubeSampling {
public:
    // Constructor: numSamples is the number of samples, numDimensions is the dimensionality of the space
    LatinHypercubeSampling(int numSamples, int numDimensions)
        : numSamples(numSamples), numDimensions(numDimensions) {
        // Initialize the random generator
        rng.seed(std::chrono::system_clock::now().time_since_epoch().count());
    }

    // Generate Latin Hypercube samples
    std::vector<std::vector<double>> generateSamples() {
        // Create a matrix to store samples
        std::vector<std::vector<double>> samples(numSamples, std::vector<double>(numDimensions));

        // Step 1: Divide the range [0, 1] into numSamples intervals for each dimension
        for (int i = 0; i < numDimensions; ++i) {
            // Generate evenly spaced intervals
            std::vector<double> intervals(numSamples);
            for (int j = 0; j < numSamples; ++j) {
                intervals[j] = (j + randomUniform(0.0, 1.0)) / numSamples;
            }

            // Step 2: Randomly permute the intervals to ensure a random distribution in each dimension
            std::shuffle(intervals.begin(), intervals.end(), rng);

            // Step 3: Assign the permuted intervals to the samples
            for (int j = 0; j < numSamples; ++j) {
                samples[j][i] = intervals[j];
            }
        }

        return samples;
    }

    // Get the number of samples
    int getSamplesNumber() const {
        return numSamples;
    }

private:
    int numSamples;
    int numDimensions;
    std::mt19937 rng;  // Random number generator

    // Generate a random number in the range [low, high)
    double randomUniform(double low, double high) {
        std::uniform_real_distribution<double> dist(low, high);
        return dist(rng);
    }
};

#endif // LATIN_HYPERCUBE_SAMPLING_H
