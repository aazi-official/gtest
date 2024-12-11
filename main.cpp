#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "tmcmc.h"

// �Զ��������Ȼ����
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

    const int ndim = 4;         
    const int nsamps = 5000;    
    const int write_flag = 1;   

    // �̶��� CoV ֵ
    std::vector<double> covValues = {pow(10, -1), pow(10, -0.75), pow(10, -0.5), pow(10, -0.25), 1.0};

    // �̶� Proposal ����
    std::vector<std::string> proposalConfigPaths = {
        "./proposal_files/proposal_1",
        "./proposal_files/proposal_2",
        "./proposal_files/proposal_3",
        "./proposal_files/proposal_4",
        "./proposal_files/proposal_5"
    };

    // �̶� Prior �ļ�·��
    std::string priorParamsFile = "prior_params.dat";
    std::string priorTypesFile = "prior_types.dat";

    LogLikelihood logLikelihood;

    // �洢�׶���
    std::vector<std::vector<int>> stages(proposalConfigPaths.size(), std::vector<int>(covValues.size(), 0));

    // ���� Proposal �� CoV ֵ
    for (size_t i = 0; i < proposalConfigPaths.size(); i++) {
        // ��ȡ Proposal �ļ�·��
        std::string proposalParamsFile = proposalConfigPaths[i] + "_params.dat";
        std::string proposalTypesFile = proposalConfigPaths[i] + "_types.dat";

        for (size_t j = 0; j < covValues.size(); j++) {
            double cov = covValues[j];
            if (pid == 0) std::cout << "Running Proposal " << i + 1 << " with CoV=" << cov << "\n";

            // ��ȡ��ǰ�׶���
            int initialStageCount = 0;
            if (pid == 0) {
                std::string command = "ls samples.dat.* 2>/dev/null | wc -l";
                FILE *pipe = popen(command.c_str(), "r");
                if (pipe) {
                    fscanf(pipe, "%d", &initialStageCount);
                    pclose(pipe);
                }
            }

            // ���� TMCMC
            double logevid = tmcmc(nsamps, 100 * pid + 4, ndim, cov,
                                   write_flag, pid, np,
                                   priorParamsFile, priorTypesFile, 
                                   proposalParamsFile, proposalTypesFile, 
                                   &logLikelihood);

            // ����׶���
            if (pid == 0) {
                int finalStageCount = 0;
                std::string command = "ls samples.dat.* 2>/dev/null | wc -l";
                FILE *pipe = popen(command.c_str(), "r");
                if (pipe) {
                    fscanf(pipe, "%d", &finalStageCount);
                    pclose(pipe);
                }
                stages[i][j] = finalStageCount - initialStageCount;
                std::cout << "Proposal " << i + 1 << " with CoV = " << cov
                          << " completed in " << stages[i][j] << " stages.\n";
            }
        }
    }

    // ����׶������ļ�
    if (pid == 0) {
        std::ofstream outfile("stages_vs_cov.dat");
        for (size_t i = 0; i < stages.size(); i++) {
            outfile << "Proposal " << i + 1 << ": ";
            for (size_t j = 0; j < stages[i].size(); j++) {
                outfile << stages[i][j] << " ";
            }
            outfile << "\n";
        }
        outfile.close();
    }

    MPI_Finalize();
    return 0;
}

