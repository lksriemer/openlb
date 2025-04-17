#ifndef POSTPROCESSING_2D_H
#define POSTPROCESSING_2D_H
using namespace olb;
using namespace olb::descriptors;

template <typename T, typename DESCRIPTOR>
void readData2D(
    int samples_number,
    const std::string& uqFolder,
    const std::string& name,
    const std::string& dataName,
    const std::vector<std::size_t>& iT,
    std::size_t rank,
    std::vector<std::vector<std::vector<std::vector<T>>>>& data,
    int& nx, int& ny, unsigned& size, Vector<T, 2>& spacing
) {
    // Read the first sample's VTI file to get geometry information
    std::string basePath = uqFolder + "0/tmp/";
    singleton::directories().setOutputDir(basePath);

    std::string fileName = singleton::directories().getVtkOutDir() + "data/" +
                           createFileName(name, iT[0], rank) + ".vti";

    BlockVTIreader2D<T, T> reader0(fileName, dataName);
    BlockData<2, T, T>& blockDataSample0 = reader0.getBlockData();

    // Read spacing from VTI file
    XMLreader xmlConfig(fileName);
    std::stringstream spacingStream(xmlConfig["ImageData"].getAttribute("Spacing"));
    for (int i = 0; i < 2; ++i) {
        spacingStream >> spacing[i];
    }

    nx = blockDataSample0.getNx();
    ny = blockDataSample0.getNy();
    size = blockDataSample0.getSize();

    // Initialize data vector
    data.assign(nx, std::vector<std::vector<std::vector<T>>>(
                        ny, std::vector<std::vector<T>>(
                                size, std::vector<T>(samples_number))));

    // Read data from each sample
    for (int iSample = 0; iSample < samples_number; ++iSample) {
        std::string samplePath = uqFolder + std::to_string(iSample) + "/tmp/";
        singleton::directories().setOutputDir(samplePath);

        std::string sampleFileName = singleton::directories().getVtkOutDir() + "data/" +
                                     createFileName(name, iT[iSample], rank) + ".vti";

        BlockVTIreader2D<T, T> reader(sampleFileName, dataName);
        BlockData<2, T, T>& blockDataSample = reader.getBlockData();

        // Extract data for each lattice point
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                for (unsigned iSize = 0; iSize < size; ++iSize) {
                    data[ix][iy][iSize][iSample] = blockDataSample.get({ix, iy}, iSize);
                }
            }
        }
    }
}

template <typename T>
void computeMeanAndStd2D(
    UncertaintyQuantification& uq,
    const std::vector<std::vector<std::vector<std::vector<T>>>>& data,
    int nx, int ny, unsigned size,
    BlockData<2, T, T>& blockDataMean,
    BlockData<2, T, T>& blockDataStd
) {
    // Compute mean and standard deviation for each lattice point
    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            for (unsigned iSize = 0; iSize < size; ++iSize) {
                const auto& samples = data[ix][iy][iSize];
                blockDataMean.get({ix, iy}, iSize) = uq.mean(samples);
                blockDataStd.get({ix, iy}, iSize) = uq.std(samples);
            }
        }
    }
}

template <typename T, typename DESCRIPTOR>
void computeMeanAndStdAndWriteVTI2D(
    UncertaintyQuantification& uq,
    const std::string& uqFolder,
    const std::string& name,
    const std::string& dataName,
    CuboidGeometry2D<T>& cuboidGeometry,
    SuperGeometry<T, 2>& sGeometry
) {
    OstreamManager clout(std::cout, "computeMeanAndStdAndWriteVTI");
    SuperLattice<T, DESCRIPTOR> sLattice(sGeometry);
    sLattice.template defineDynamics<NoDynamics>(sGeometry, 0);
    sLattice.initialize();

    // Load iteration logs for each sample
    std::vector<std::vector<size_t>> iTList(uq.getSamplesNumber());
    size_t numIterations = std::numeric_limits<size_t>::max();

    for (size_t n = 0; n < uq.getSamplesNumber(); ++n) {
        std::string filePath = uqFolder + std::to_string(n) + "/tmp/iteration_log.txt";
        std::ifstream inFile(filePath);
        if (inFile.is_open()) {
            size_t i;
            while (inFile >> i) {
                iTList[n].push_back(i);
            }
            inFile.close();
            numIterations = std::min(numIterations, iTList[n].size());
        } else {
            clout << "Could not open file: " << filePath << std::endl;
        }
    }

    SuperVTMwriter2D<T> vtmWriter(name);

#ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
#else
    const int noOfCuboids = 1;
#endif

    // Rename materials for visualization
    for (int rank = 0; rank < noOfCuboids; ++rank) {
        Cuboid2D<T>& cuboid2d = cuboidGeometry.get(rank);
        Vector<T, 2> origin = cuboid2d.getOrigin();
        origin[0] -= cuboid2d.getDeltaR();
        origin[1] -= cuboid2d.getDeltaR();
        Vector<T, 2> extent = cuboid2d.getExtent() * cuboid2d.getDeltaR();
        // extent[0] += 2 * cuboid2d.getDeltaR();
        // extent[1] += 2 * cuboid2d.getDeltaR();
        IndicatorCuboid2D<T> readCuboid(extent, origin);
        sGeometry.rename(0, 10 + rank, readCuboid);
    }

    // Process each iteration
    for (size_t iter = 0; iter < numIterations; ++iter) {
        std::vector<size_t> iT;
        for (const auto& list : iTList) {
            iT.push_back(list[iter]);
        }

        clout << "Processing iteration " << iT[0] << std::endl;

        if (iT[0] == 0) {
            // Write initial geometry for visualization
            singleton::directories().setOutputDir("./tmp/");
            SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
            SuperLatticeRank2D<T, DESCRIPTOR> rank(sLattice);
            vtmWriter.write(cuboid);
            vtmWriter.write(rank);
            vtmWriter.createMasterFile();
        }

        // Process data for each rank
        for (int rank = 0; rank < noOfCuboids; ++rank) {
            std::string basePath = uqFolder + "0/tmp/";
            singleton::directories().setOutputDir(basePath);

            std::string fileName = singleton::directories().getVtkOutDir() + "data/" +
                                   createFileName(name, iT[0], rank) + ".vti";

            BlockVTIreader2D<T, T> readerMean(fileName, dataName);
            BlockVTIreader2D<T, T> readerStd(fileName, dataName);
            BlockData<2, T, T>& blockDataMean = readerMean.getBlockData();
            BlockData<2, T, T>& blockDataStd = readerStd.getBlockData();

            // Read and process data
            std::vector<std::vector<std::vector<std::vector<T>>>> data;
            int nx, ny;
            unsigned size;
            Vector<T, 2> spacing;

            readData2D<T, DESCRIPTOR>(
                uq.getSamplesNumber(), uqFolder, name, dataName, iT, rank, data, nx, ny, size, spacing
            );

            computeMeanAndStd2D<T>(uq, data, nx, ny, size, blockDataMean, blockDataStd);

            // Create analytical fields
            BlockDataF2D<T, T> blockDataFMean(blockDataMean);
            BlockDataF2D<T, T> blockDataFStd(blockDataStd);

            SpecialAnalyticalFfromBlockF2D<T, T> meanField(blockDataFMean, cuboidGeometry.get(rank), spacing, 1.0);
            SpecialAnalyticalFfromBlockF2D<T, T> stdField(blockDataFStd, cuboidGeometry.get(rank), spacing, 1.0);

            // Define fields on the lattice
            auto materialIndicator = sGeometry.getMaterialIndicator({10 + rank});
            for (int iBalancer = 0; iBalancer < sGeometry.getLoadBalancer().size(); ++iBalancer) {
                sLattice.getBlock(iBalancer).template defineField<descriptors::VELOCITY>(
                    materialIndicator->getBlockIndicatorF(iBalancer), meanField
                );
                sLattice.getBlock(iBalancer).template defineField<descriptors::VELOCITY2>(
                    materialIndicator->getBlockIndicatorF(iBalancer), stdField
                );
            }
        }

        // Write data to VTM file
        SuperLatticePhysField2D<T, DESCRIPTOR, VELOCITY> meanPhysField(sLattice, 1.0, dataName + "Mean");
        vtmWriter.addFunctor(meanPhysField);

        SuperLatticePhysField2D<T, DESCRIPTOR, VELOCITY2> stdPhysField(sLattice, 1.0, dataName + "Std");
        vtmWriter.addFunctor(stdPhysField);

        singleton::directories().setOutputDir("./tmp/");
        vtmWriter.write(iT[0]);
    }
}

#endif // POSTPROCESSING_2D_H