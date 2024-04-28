#include <boost/program_options.hpp>
#include <fstream>
#include <sys/time.h>

#include "eqChecker.h"
#include "memMeasure.h"

extern void qasmParser(std::ifstream &inFile, std::vector<GateType> &gates,
                       std::vector<std::vector<int>> &qubits, int &n);

void partialEquivalencCheckingBenchmarks(const int minN, const int maxN,
                                         const size_t reps,
                                         const bool addAncilla, EqType eqType,
                                         const std::string &filename,
                                         const std::string &directoryname) {

  std::string rootDirectoryBenchmarks = "../../benchmarks/";

  std::fstream log2(rootDirectoryBenchmarks + "sliqec_benchmark_log.txt",
                    std::ios::out | std::ios::app);
  log2 << "starting benchmark.\nminN: " << minN << ", maxN: " << maxN
       << ", reps: " << reps << ", addAncilla: " << addAncilla
       << ", filename: " << filename << "\n";
  log2.close();

  for (int d = minN; d < maxN; d += 5) {
    double totalTime{0};
    std::size_t timeouts{0};
    int n{0};
    if (addAncilla) {
      n = static_cast<int>(1.5 * d);
    } else {
      n = d;
    }
    const auto m = static_cast<int>(0.5 * d);
    for (size_t k = 0; k < reps; k++) {

      // Parse QASM files
      std::vector<std::vector<GateType>> gates(2);
      std::vector<std::vector<std::vector<int>>> qubits(2);
      int nQ1, nQ2;

      std::string circuitsFilename = std::to_string(d) + "_" +
                                     std::to_string(m) + "_" +
                                     std::to_string(k) + ".qasm";

      std::ifstream inFile;

      inFile.open(rootDirectoryBenchmarks + directoryname + "a/" +
                  circuitsFilename);
      if (!inFile) {
        std::fstream log2(rootDirectoryBenchmarks + "sliqec_benchmark_log.txt",
                          std::ios::out | std::ios::app);
        log2 << "Circuit1 file doesn't exist. "
             << rootDirectoryBenchmarks + directoryname + "a/" +
                    circuitsFilename
             << std::endl;
        log2.close();
        return;
      }
      qasmParser(inFile, gates[0], qubits[0], nQ1);
      inFile.close();

      inFile.open(rootDirectoryBenchmarks + directoryname + "b/" +
                  circuitsFilename);
      if (!inFile) {
        std::fstream log2(rootDirectoryBenchmarks + "sliqec_benchmark_log.txt",
                          std::ios::out | std::ios::app);
        log2 << "Circuit2 file doesn't exist." << std::endl;
        log2.close();
        return;
      }
      qasmParser(inFile, gates[1], qubits[1], nQ2);
      inFile.close();

      struct timeval tStart, tFinish;
      double elapsedTime;
      double runtime;
      size_t memPeak;

      gettimeofday(&tStart, NULL);

      EquivalenceChecker checker(gates, qubits, n, d, m, false, eqType);

      checker.check();

      gettimeofday(&tFinish, NULL);
      elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
      elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

      runtime = elapsedTime / 1000.0;
      std::fstream log(rootDirectoryBenchmarks + "sliqec_benchmark_log.txt",
                       std::ios::out | std::ios::app);
      if (runtime <= 600) {
        totalTime += runtime;
      } else {
        log << "TIMEOUT; ";
        timeouts++;
      }
      log << "k: " << k << ", d: " << d << ", m: " << m << ", time: " << runtime
          << "\n";
      log.close();
    }
    std::fstream resultsFile(rootDirectoryBenchmarks + "sliqec_" + filename,
                             std::ios::out | std::ios::app);
    resultsFile << "" << n << "," << d << "," << m << "," << reps << ","
                << (totalTime / static_cast<double>(reps - timeouts)) << ","
                << timeouts << "," << reps - timeouts << "\n";

    if (timeouts >= reps - 3) {
      resultsFile << "stopping because of too many timeouts\n";
      resultsFile.close();
      return;
    }
    resultsFile.close();
  }
}

int main(int argc, char **argv) {

  std::string directory1 = "benchmarkCircuitsConstruction";
  std::string directory2 = "benchmarkCircuitsConstructionNoAncilla";
  std::string directory3 = "benchmarkCircuitsAlternating";

  EqType eqType;

  // algorithm 2
  int minN = 5;
  int maxN = 31;
  int reps = 20;
  partialEquivalencCheckingBenchmarks(minN, maxN, reps, true, EqType::Peq,
                                      "construction_benchmarks.txt",
                                      directory1);
  partialEquivalencCheckingBenchmarks(minN, maxN, reps, false, EqType::Peq,
                                      "construction_no_ancilla.txt",
                                      directory2);

  // algorithm 3
  minN = 5;
  maxN = 101;
  reps = 20;
  partialEquivalencCheckingBenchmarks(minN, maxN, reps, false, EqType::PeqS,
                                      "alternating_benchmark.txt", directory3);
  return 0;
}
