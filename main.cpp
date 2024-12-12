#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "tmcmc.h"

// 自定义对数似然函数
class LogLikelihood : public LogDensity {
public:
    LogLikelihood() : LogDensity() {};

    double eval(RealVector &x) override {
        int ndim = x.size();
        double loglik = 0.0;
        double gamma = 1.0;

        for (int j = 0; j < ndim; j++) {
            loglik += -0.5 * log(2.0 * Pi) - log(gamma) - 0.5 * pow(x[j] / gamma, 2.0);
        }
        return loglik;
    }
};

int main(int argc, char *argv[]) {
    int pid, ierr, np;
    ierr = MPI_Init(&argc, &argv);
    if (ierr != 0) {
        std::cerr << "Error initializing MPI environment.\n";
        exit(1);
    }

    ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    const int ndim = 4;          // 4维问题
    const int nsamps = 500;      // 采样点数
    const int write_flag = 0;    // 不写入中间结果

    // 固定的 CoV 值
    double cov = 0.1;  // 选定单个 CoV 值

    // 固定 Proposal 配置
    std::string proposalParamsFile = "./proposal_files/proposal_1_params.dat";
    std::string proposalTypesFile = "./proposal_files/proposal_1_types.dat";

    // 固定 Prior 文件路径
    std::string priorParamsFile = "prior_params.dat";
    std::string priorTypesFile = "prior_types.dat";

    LogLikelihood logLikelihood;

    if (pid == 0) std::cout << "Running Proposal 1 with CoV=" << cov << "\n";

    // 调用 TMCMC
    double logevid = tmcmc(nsamps, 100 * pid + 4, ndim, cov,
                           write_flag, pid, np,
                           priorParamsFile, priorTypesFile, 
                           proposalParamsFile, proposalTypesFile, 
                           &logLikelihood);

    if (pid == 0) {
        std::cout << "Proposal 1 with CoV = " << cov
                  << " completed. Log Evidence = " << logevid << "\n";

        // 输出结果到文件
        std::ofstream outfile("result_summary.dat");
        outfile << "Proposal 1 with CoV = " << cov
                << " completed. Log Evidence = " << logevid << "\n";
        outfile.close();
    }

    MPI_Finalize();
    return 0;
}


