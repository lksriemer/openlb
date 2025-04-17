#include "tgv2d.h"
#include <time.h>
#include "../uq/uq.h"
#include "../uq/postprocessing2D.h"

#define SC
// #define MC

using namespace std;

std::vector<std::vector<double>> loadMatrix(const std::string &filename)
{
    std::ifstream inFile(filename);
    if (!inFile) {
        throw std::runtime_error("Error: could not open file " + filename);
    }

    std::vector<std::vector<double>> data;
    std::string line;

    while (std::getline(inFile, line)) {
        // Optionally skip empty lines
        if (line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        double val;
        std::vector<double> row;

        // Read the values in this line
        while (iss >> val) {
            row.push_back(val);
        }

        if (!row.empty()) {
            data.push_back(row);
        }
    }

    return data;
}


void totalKineticEnergy(std::vector<double>&tke, int t, double u0, double dx, vector<vector<vector<T>>> u, vector<vector<vector<T>>> v, std::unique_ptr<GeneralizedPolynomialChaos>& op)
{
  OstreamManager clout( std::cout,"output" );
  size_t No = op->getPolynomialsOrder();
  std::vector<T> u2(No, 0.0);
  std::vector<T> v2(No, 0.0);
  int nx = u.size();
  int ny = u.size();

  for (int i = 0; i < nx; ++i){
    for (int j = 0; j < ny; ++j){
      op->chaosProduct(u[i][j], u[i][j], u2);
      op->chaosProduct(v[i][j], v[i][j], v2);

      for (int alpha = 0; alpha < No; ++alpha) {
        tke[alpha] += ((u2[alpha] + v2[alpha]) *  2 / (nx*ny*u0*u0));
      }
    }
  }
}

std::vector<T> MonteCarloPostprocessing(string foldPath, int nq, int resolution, double physVelocity, double dx, UncertaintyQuantification& uq) {
    OstreamManager clout( std::cout,"main" );
    clout << "Starting Monte Carlo postprocessing" << endl;
    std::vector<T> tke(nq, 0.0);
    std::vector<std::vector<T>> u_mean(resolution, std::vector<T>(resolution, 0.0));
    std::vector<std::vector<T>> v_mean(resolution, std::vector<T>(resolution, 0.0));
    std::vector<std::vector<T>> u_var(resolution, std::vector<T>(resolution, 0.0));
    std::vector<std::vector<T>> v_var(resolution, std::vector<T>(resolution, 0.0));
    std::vector<std::vector<T>> u_mean_old(resolution, std::vector<T>(resolution, 0.0));
    std::vector<std::vector<T>> v_mean_old(resolution, std::vector<T>(resolution, 0.0));
    std::vector<std::vector<T>> u_var_old(resolution, std::vector<T>(resolution, 0.0));
    std::vector<std::vector<T>> v_var_old(resolution, std::vector<T>(resolution, 0.0));

    for (int i = 0; i < nq; ++i) {
      string file_u = foldPath + "dataFiles/5/u_" + std::to_string(i) + ".dat";
      string file_v = foldPath + "dataFiles/5/v_" + std::to_string(i) + ".dat";

      std::vector<std::vector<T>> u = loadMatrix(file_u);
      std::vector<std::vector<T>> v = loadMatrix(file_v);


      for (int ii = 0; ii < resolution; ++ii) {
        for (int jj = 0; jj < resolution; ++jj) {
          tke[i] += ((u[ii][jj] * u[ii][jj] + v[ii][jj] * v[ii][jj]) * 0.5 / ( resolution * resolution * physVelocity * physVelocity ));

          u_mean_old[ii][jj] = u_mean[ii][jj];
          u_mean[ii][jj] = u_mean_old[ii][jj] + (u[ii][jj] - u_mean_old[ii][jj])/(i+1);
          u_var_old[ii][jj] = u_var[ii][jj];
          u_var[ii][jj] = u_var_old[ii][jj] + (u[ii][jj] - u_mean_old[ii][jj])*(u[ii][jj] - u_mean[ii][jj]);

          v_mean_old[ii][jj] = v_mean[ii][jj];
          v_mean[ii][jj] = v_mean_old[ii][jj] + (v[ii][jj] - v_mean_old[ii][jj])/(i+1);
          v_var_old[ii][jj] = v_var[ii][jj];
          v_var[ii][jj] = v_var_old[ii][jj] + (v[ii][jj] - v_mean_old[ii][jj]) * (v[ii][jj] - v_mean[ii][jj]);
        }
      }
    }

    clout << std::fixed << std::setprecision(20);
    clout << "Total kinetic energy: mean = " << uq.mean(tke) << ", std = " << uq.std(tke) << endl;

    std::string filenameUAll = "uq/dataFiles/5/final/tgv2dMc.dat";
    std::ofstream outputFileUAll(filenameUAll);

    if (!outputFileUAll) {
      std::cerr << "Error opening the file: " << filenameUAll << std::endl;
    }

    outputFileUAll << "TITLE = \"TGV2D\"\n";
    outputFileUAll << "VARIABLES = \"X\", \"Y\", \"umean\", \"ustd\", \"vmean\", \"vstd\" \n";
    outputFileUAll << "ZONE T=\"Zone 1\", I=" << resolution << ", J=" << resolution << ", DATAPACKING=POINT\n";

    for (int i = 0; i < resolution; ++i) {
      for (int j = 0; j < resolution; ++j) {
          outputFileUAll.precision(20);
          outputFileUAll << i * dx << "\t" << j * dx << "\t" << u_mean[i][j] << "\t" << sqrt(u_var[i][j]) << "\t" << v_mean[i][j] << "\t" << sqrt(v_var[i][j]) << "\n";
      }
    }
    outputFileUAll.close();

    return tke;

}

std::vector<T> stochasticCollocationPostprocessing(string foldPath, int order, int nq, int resolution, double physVelocity, double dx, UncertaintyQuantification& uq) {
    OstreamManager clout( std::cout,"main" );

    std::vector<std::vector<std::vector<T>>> u(nq, std::vector<std::vector<T>>(resolution, std::vector<T>(resolution, 0.0)));
    std::vector<std::vector<std::vector<T>>> v(nq, std::vector<std::vector<T>>(resolution, std::vector<T>(resolution, 0.0)));
    std::vector<T> tke(order+1, 0.0);

    for (int i = 0; i < nq; ++i) {
      string file_u = foldPath + "dataFiles/5/u_" + std::to_string(i) + ".dat";
      string file_v = foldPath + "dataFiles/5/v_" + std::to_string(i) + ".dat";
      u[i] = loadMatrix(file_u);
      v[i] = loadMatrix(file_v);
    }
    // clout << "Loaded velocity fields" << endl;

    std::unique_ptr<GeneralizedPolynomialChaos> op = uq.getOps();

    // clout << "No: " << op->getPolynomialsOrder() << ", nq: " << op->getQuadraturePointsNumber() << endl;

    std::vector<std::vector<std::vector<T>>> uChaos(resolution, std::vector<std::vector<T>>(resolution, std::vector<T>(order+1, 0.0)));
    std::vector<std::vector<std::vector<T>>> vChaos(resolution, std::vector<std::vector<T>>(resolution, std::vector<T>(order+1, 0.0)));

    for (int i = 0; i < resolution; ++i) {
      for (int j = 0; j < resolution; ++j) {
        std::vector<T> utmp(nq, 0.0);
        std::vector<T> vtmp(nq, 0.0);
        for (int k = 0; k < nq; ++k) {
          utmp[k] = u[k][i][j];
          vtmp[k] = v[k][i][j];
        }
        op->randomToChaos(utmp, uChaos[i][j]);
        op->randomToChaos(vtmp, vChaos[i][j]);
      }
    }
    // clout << "Converted velocity fields to chaos" << endl;

    std::vector<T> u2(order+1, 0.0);
    std::vector<T> v2(order+1, 0.0);
    for (int i = 0; i < resolution; ++i){
      for (int j = 0; j < resolution; ++j){
        op->chaosProduct(uChaos[i][j], uChaos[i][j], u2);
        op->chaosProduct(vChaos[i][j], vChaos[i][j], v2);

        for (int alpha = 0; alpha < order+1; ++alpha) {
          tke[alpha] += ((u2[alpha] + v2[alpha]) * 0.5 / ( resolution * resolution * physVelocity * physVelocity ));
        }
      }
    }

    clout << std::fixed << std::setprecision(20);
    clout << "Total kinetic energy: mean = " << op->mean(tke) << ", std = " << op->std(tke) << endl;

    std::string filenameUAll = "uq/dataFiles/5/final/tgv2dSc.dat";
    std::ofstream outputFileUAll(filenameUAll);

    if (!outputFileUAll) {
      std::cerr << "Error opening the file: " << filenameUAll << std::endl;
    }

    outputFileUAll << "TITLE = \"TGV2D\"\n";
    outputFileUAll << "VARIABLES = \"X\", \"Y\", \"umean\", \"ustd\", \"vmean\", \"vstd\" \n";
    outputFileUAll << "ZONE T=\"Zone 1\", I=" << resolution << ", J=" << resolution << ", DATAPACKING=POINT\n";

    for (int i = 0; i < resolution; ++i) {
      for (int j = 0; j < resolution; ++j) {
        // std::vector<T> utmp(nq, 0.0);
        // std::vector<T> vtmp(nq, 0.0);
        // for (int n = 0; n < nq; ++n) {
        //   utmp[n] = u[n][i][j];
        //   vtmp[n] = v[n][i][j];
        // }
        // clout << endl;
          outputFileUAll.precision(20);
          // outputFileUAll << i * dx << "\t" << j * dx << "\t" << uq.mean(utmp) << "\t" << uq.std(utmp) << "\t" << uq.mean(vtmp) << "\t" << uq.std(vtmp) << "\n";
          outputFileUAll << i * dx << "\t" << j * dx << "\t" << op->mean(uChaos[i][j]) << "\t" << op->std(uChaos[i][j]) << "\t" << op->mean(vChaos[i][j]) << "\t" << op->std(vChaos[i][j]) << "\n";
      }
    }
    outputFileUAll.close();

  return tke;
}

bool createDirectories(const std::string& foldPath) {
    OstreamManager clout( std::cout,"main" );
    std::string command;
    int mkdirResult;

    // Create the first directory
    command = "mkdir -p " + foldPath;
    mkdirResult = std::system(command.c_str());
    if (mkdirResult != 0) {
        clout << "Failed to create directory: " << foldPath << std::endl;
        return false; // Early return on failure
    }

    for ( int i = 0; i <= 5; ++i ) {
      std::string subFoldPath = foldPath + "/" + std::to_string(i);
      command = "mkdir -p " + subFoldPath;
      mkdirResult = std::system(command.c_str());
      if (mkdirResult != 0) {
          clout << "Failed to create directory: " << subFoldPath << std::endl;
          return false; // Early return on failure
      }

      // Create the second directory
      std::string foldPathFinal = subFoldPath + "/final";
      clout << "foldPath: " << foldPath << std::endl;
      command = "mkdir -p " + foldPathFinal;
      mkdirResult = std::system(command.c_str());
      if (mkdirResult != 0) {
          clout << "Failed to create directory: " << foldPathFinal << std::endl;
          return false; // Early return on failure
      }
    }

    // If both directories are created successfully
    clout << "Directories created successfully:" << std::endl;
    clout << foldPath << std::endl;

    return true;
}

void save_viscosity_list(std::string dir, std::vector<double> list) {
  std::string filename = dir + "/viscosity_list.dat";
  std::ofstream outputFile(filename);
  if (!outputFile) {
    std::cerr << "Error opening the file: " << filename << std::endl;
    return;
  }
  for (int i = 0; i < list.size(); ++i) {
    outputFile << list[i] << "\n";
  }
  outputFile.close();
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  OstreamManager clout( std::cout,"main" );

  string foldPath = "./uq/";
  string command;
  command = "mkdir -p " + foldPath;
  int result = std::system(command.c_str());
  if (result != 0) {
      clout << "Command failed with return code: " << result << std::endl;
  }
  else {
      clout << "finish mkdir" << std::endl;
  }


    // Parameters params;
    // readParameters("./parameters.dat", params);
    double num_vortice = 2;
    double L = 1.0;
    double resolution = 33;
    double physVelocity = 0.01;
    double charL = num_vortice * M_PI * L;
    double Ma = 0.0125;
    double Re = 15;

    int order = 3;
    int nq = 7;

    // Parse command-line arguments
    if (argc > 1) {
      order = std::stoi(argv[1]);
    }
    if (argc > 2) {
      nq = std::stoi(argv[2]);
    }
    if (argc > 2) {
      resolution = std::stoi(argv[3]);
    }

    double dx = charL / resolution;
    double dt = dx / (physVelocity / (Ma / std::sqrt(3)));
    double conversionViscosity = dx * dx / dt;
    double physViscosity = physVelocity * charL / Re;
    double tau = physViscosity * 3 + 0.5;

  #ifdef SC
    std::vector<Distribution> distributions = { Distribution(DistributionType::Uniform, 0.8 * Re, 1.2 * Re) };
    // // std::vector<Distribution> distributions = { Distribution(DistributionType::Normal, Re, 0.1 * Re) };
    UncertaintyQuantification uq(UQMethod::GPC);
    uq.initializeGPC(order, nq, distributions, Quadrature::QuadratureMethod::WilkinsonShiftQR );
    foldPath += "sc/res" + std::to_string(static_cast<int>(resolution)) + "/";
  #elif defined(MC)
    UncertaintyQuantification uq(UQMethod::MonteCarlo);
    unsigned int seed = 123456;
    nq = 5000;
    uq.initializeMonteCarlo(nq, 1, Distribution(DistributionType::Uniform, 0.8 * Re, 1.2 * Re), seed);
    foldPath += "mc/res" + std::to_string(static_cast<int>(resolution)) + "/";
  #endif

  std::vector<std::vector<double>> samples;
  uq.getSamplingPoints(samples);

  if (createDirectories(foldPath + "dataFiles/")) {
    clout << "Both directories were created successfully." << std::endl;
  } else {
    clout << "Error creating directories." << std::endl;
  }

    singleton::directories().setOutputDir( foldPath + "/tmp/" );


    clout << "tau: " << tau << "; conversionViscosity: " << conversionViscosity <<std::endl;
    // double tau = params.tau;


    clock_t start,end;
    start = clock();
    clout << "Start simulation" << endl;

  for(size_t n = 0; n < samples.size(); ++n) {
    simulateTGV( samples[n][0], charL, tau, resolution, foldPath, n);
  }

    end = clock();
    clout << "total CPI time used: " << (double)(end-start)/CLOCKS_PER_SEC << "s" << endl;

    clout << "Finished simulation" << endl;

  // for (int order = 1; order <= 8; ++order) {
  //   #if defined(SC)
  //     uq.initializeGPC(order, nq, distributions, Quadrature::QuadratureMethod::GSL );
  //     std::vector<T> tke = stochasticCollocationPostprocessing(foldPath, order, nq, resolution, physVelocity, dx, uq);
  //   #elif defined(MC)
  //     std::vector<T> tke = MonteCarloPostprocessing(foldPath, nq, resolution, physVelocity, dx, uq);
  //   #endif
  // }
  #ifdef SC
    std::vector<T> tke = stochasticCollocationPostprocessing(foldPath, order, nq, resolution, physVelocity, dx, uq);
  #elif defined(MC)
    std::vector<int> nqList = {10, 100, 1000, 5000};
    for (int nq : nqList) {
        std::vector<T> tke = MonteCarloPostprocessing(foldPath, nq, resolution, physVelocity, dx, uq);
    }
  #endif

  return 0;
}