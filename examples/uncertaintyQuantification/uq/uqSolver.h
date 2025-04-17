#ifndef UQ_SOLVER_H
#define UQ_SOLVER_H

#include "lbSolver.h"
#include "unitConverter.h" // Assuming this contains UnitConverterFromResolutionAndRelaxationTime

namespace olb {

template <typename T, typename PARAMETERS, typename LATTICES>
class uqSolver : public LbSolver<T, PARAMETERS, LATTICES> {

private:
    mutable OstreamManager clout {std::cout, "uqSolver"};
    UnitConverterFromResolutionAndRelaxationTime<T, LATTICES> customConverter;

public:
    uqSolver(utilities::TypeIndexedSharedPtrTuple<PARAMETERS> params)
        : LbSolver<T, PARAMETERS, LATTICES>(params) {
        // You can initialize additional members here if needed
    }

    /// Override initialize() to use custom converter
    void initializeConverter() override {
        // Set your custom parameters for the converter
        auto geomParams = this->parameters(names::Geometry());  // Assuming geomParams can be accessed here
        int N = 100;  // Set resolution value
        T Re = 1000;  // Set Reynolds number or other relevant parameters

        // Custom converter definition
        customConverter = UnitConverterFromResolutionAndRelaxationTime<T, LATTICES>(
            N,                                      // resolution: number of voxels per charPhysL
            static_cast<T>(0.56),                   // latticeRelaxationTime: relaxation time
            static_cast<T>(2.0 * geomParams.radiusCylinder),   // charPhysLength: reference length
            static_cast<T>(0.2),                    // charPhysVelocity: maximal expected velocity
            static_cast<T>(0.2 * 2.0 * geomParams.radiusCylinder / Re),  // physViscosity
            static_cast<T>(1.0)                     // physDensity: physical density
        );

        clout << "uqSolver initialized with custom converter." << std::endl;

        // Call the base class initialize to set up the remaining things like geometry and lattices
        LbSolver<T, PARAMETERS, LATTICES>::initialize();
    }

    /// New function: outputField
    void outputField(std::size_t iT) {
        clout << "Output field at time step " << iT << std::endl;
        // Logic to output the field at a given timestep
    }

    // You can access the custom converter like this
    auto& converter() {
        return customConverter;
    }
};

} // namespace olb

#endif  // UQSOLVER_H
