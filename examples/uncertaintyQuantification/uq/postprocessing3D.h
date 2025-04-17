#ifndef POSTPROCESSING_3D_H
#define POSTPROCESSING_3D_H

#include <olb.h>
// #include "../../src/io/vtiReader.h"
#include "uq.h"

using namespace olb;
using namespace olb::descriptors;

// Template over T and DESCRIPTOR
template <typename T, typename DESCRIPTOR>
void readData3D(int samples_number,
                const std::string& uqFolder,
                const std::string& _name,
                const std::string& dataName,
                const std::vector<std::size_t>& iT,
                std::size_t iC,
                std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>& data,
                int& nx, int& ny, int& nz, unsigned& size, Vector<T,3>& spacing) {

    OstreamManager clout( std::cout,"readData3D" );
    // Read the first sample's VTI file to get geometry information
    std::string subFoldPath = uqFolder + std::to_string(0) + "/tmp/";
    singleton::directories().setOutputDir( subFoldPath );

    BlockVTIreader3D<T, T> reader0(singleton::directories().getVtkOutDir() + "data/" + createFileName( _name, iT[0], iC) + ".vti", dataName);
    BlockData<3, T, T>& blockDataSample0 = reader0.getBlockData();

    // Read spacing from vti file
    XMLreader MRI_config(singleton::directories().getVtkOutDir() + "data/" + createFileName( _name, iT[0], iC) + ".vti");
    std::stringstream stream_val(MRI_config["ImageData"].getAttribute("Spacing"));

    for (int i=0; i<3; i++) {
        stream_val >> spacing[i];
    }

    nx = blockDataSample0.getNx();
    ny = blockDataSample0.getNy();
    nz = blockDataSample0.getNz();
    size = blockDataSample0.getSize();

    data.resize(nx);
    for (int ix = 0; ix < nx; ++ix) {
        data[ix].resize(ny);
        for (int iy = 0; iy < ny; ++iy) {
            data[ix][iy].resize(nz);
            for (int iz = 0; iz < nz; ++iz) {
                data[ix][iy][iz].resize(size);
                for (unsigned iSize = 0; iSize < size; ++iSize) {
                    data[ix][iy][iz][iSize].resize(samples_number);
                }
            }
        }
    }

    // For each sample
    for (int iSample = 0; iSample < samples_number; ++iSample) {
        // Open the sample .vti file
        subFoldPath = uqFolder + std::to_string(iSample) + "/tmp/";
        singleton::directories().setOutputDir( subFoldPath );
        BlockVTIreader3D<T, T> reader(singleton::directories().getVtkOutDir() + "data/" + createFileName( _name, iT[iSample], iC) + ".vti", dataName);

        clout << "Reading sample: " << singleton::directories().getVtkOutDir() + "data/" + createFileName( _name, iT[iSample], iC) + ".vti" << std::endl;

        BlockData<3, T, T>& blockDataSample = reader.getBlockData();

        // For each lattice point
        for (int iz = 0; iz < nz; ++iz) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int ix = 0; ix < nx; ++ix) {
                    for (unsigned iSize = 0; iSize < size; ++iSize) {
                        data[ix][iy][iz][iSize][iSample] = blockDataSample.get({ix, iy, iz}, iSize);
                    }
                }
            }
        }
    }
}

// Function to compute mean and std of the field
template <typename T, typename DESCRIPTOR>
void computeMeanAndStd3D(UncertaintyQuantification& uq,
                         const std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>& data,
                         int nx, int ny, int nz, unsigned size,
                         BlockData<3, T, T>& blockDataMean,
                         BlockData<3, T, T>& blockDataStd) {

    int samples_number = data[0][0][0][0].size(); // Assuming data is non-empty

    // For each lattice point
    for (int iz = 0; iz < nz; ++iz) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                for (unsigned iSize = 0; iSize < size; ++iSize) {

                    // Extract the vector of samples for this point
                    std::vector<T> samples(samples_number);
                    for (int iSample = 0; iSample < samples_number; ++iSample) {
                        samples[iSample] = data[ix][iy][iz][iSize][iSample];
                    }

                    // Use UQ class to compute mean and standard deviation
                    blockDataMean.get({ix, iy, iz}, iSize) = uq.mean(samples);
                    blockDataStd.get({ix, iy, iz}, iSize) = uq.std(samples);
                }
            }
        }
    }
}

template <typename T, typename DESCRIPTOR>
void computeMeanAndStdAndWriteVTI3D(UncertaintyQuantification& uq,
                                    const std::string& uqFolder,
                                    const std::string& _name,
                                    const std::string& dataName,
                                    SuperGeometry<T,3>& sGeometry) {

    OstreamManager clout( std::cout,"computeMeanAndStdAndWriteVTI3D" );
    SuperLattice<T,DESCRIPTOR> sLattice(sGeometry);
    sLattice.defineDynamics<NoDynamics>(sGeometry, 0);
    sLattice.initialize();

    // Load iT lists for each sample
    std::vector<std::vector<size_t>> iTList(uq.getSamplesNumber());

    for (int n = 0; n < uq.getSamplesNumber(); ++n) {
        std::string filePath = uqFolder + std::to_string(n) + "/tmp/" + "iteration_log.txt";
        std::ifstream inFile(filePath);
        if (inFile.is_open()) {
            size_t i;
            while (inFile >> i) {
                iTList[n].push_back(i);
            }
            inFile.close();
        } else {
            clout << "Could not open file: " << filePath << std::endl;
        }
    }

    size_t numIterations = iTList[0].size();
    for (size_t n = 1; n < iTList.size(); ++n) {
        if (iTList[n].size() < numIterations) {
            numIterations = iTList[n].size();
        }
    }

    SuperVTMwriter3D<T> vtmWriter(_name);

    for (int iC=0; iC < singleton::mpi().getSize(); ++iC) {
        for (size_t i = 0; i < numIterations; ++i) {

            std::vector<size_t> iT;
            for (size_t n = 0; n < iTList.size(); ++n) {
                iT.push_back(iTList[n][i]);
            }

            clout << "Processing iteration " << iT[0] << std::endl;

            if ( iT[0] == 0 ) {
                // Writes the geometry, cuboid no. and rank no. as vti file for visualization
                clout << "Writing geometry, cuboid no. and rank no. as vti file for visualization" << std::endl;
                singleton::directories().setOutputDir( "./tmp/" );
                SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
                SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
                vtmWriter.write( cuboid );
                vtmWriter.write( rank );

                vtmWriter.createMasterFile();
            }

            // Read data
            std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>> data;
            int nx, ny, nz;
            unsigned size;
            Vector<T,3> spacing;

            readData3D<T, DESCRIPTOR>(uq.getSamplesNumber(), uqFolder, _name, dataName, iT, iC, data, nx, ny, nz, size, spacing);

            std::string subFoldPath = uqFolder + std::to_string(0) + "/tmp/";
            singleton::directories().setOutputDir( subFoldPath );

            BlockVTIreader3D<T, T> reader0(singleton::directories().getVtkOutDir() + "data/" + createFileName( _name, iT[0], iC) + ".vti", dataName);
            BlockData<3, T, T>& blockDataMean = reader0.getBlockData();
            BlockData<3, T, T>& blockDataStd = reader0.getBlockData();
            Cuboid3D<T>& cuboid3d = reader0.getCuboid();

            // Compute mean and std
            computeMeanAndStd3D<T, DESCRIPTOR>(uq, data, nx, ny, nz, size, blockDataMean, blockDataStd);

            // Create BlockDataF3D objects
            BlockDataF3D<T,T> blockDataFMean(blockDataMean);
            BlockDataF3D<T,T> blockDataFStd(blockDataStd);

            // Create SpecialAnalyticalFfromBlockF3D objects
            SpecialAnalyticalFfromBlockF3D<T,T> meanField(blockDataFMean, cuboid3d, spacing, 1.0);
            SpecialAnalyticalFfromBlockF3D<T,T> stdField(blockDataFStd, cuboid3d, spacing, 1.0);

            // Define fields on the lattice
            clout << "Defining mean field" << std::endl;
            sLattice.getBlock(iC).defineField<descriptors::VELOCITY>(
                sGeometry.getMaterialIndicator({1})->getBlockIndicatorF(iC),
                meanField);

            clout << "Defining std field" << std::endl;
            sLattice.getBlock(iC).defineField<descriptors::VELOCITY>(
                sGeometry.getMaterialIndicator({1})->getBlockIndicatorF(iC),
                stdField);

            // Write data
            singleton::directories().setOutputDir( "./tmp/" );

            SuperLatticePhysField3D<T,DESCRIPTOR,VELOCITY> meanPhysField(sLattice, 1.0, dataName + "Mean");
            vtmWriter.addFunctor(meanPhysField);

            SuperLatticePhysField3D<T,DESCRIPTOR,VELOCITY> stdPhysField(sLattice, 1.0, dataName + "Std");
            vtmWriter.addFunctor(stdPhysField);

            vtmWriter.write( iT[0] );
        }
    }
}

#endif // POSTPROCESSING_3D_H
