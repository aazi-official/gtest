#ifndef TMCMCHEADERSEEN
#define TMCMCHEADERSEEN
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>
#include <functional>
#include <string>
#include <random>
#include <mpi.h>
#include <assert.h>

extern double* boundary;
const double SMALL_LPRIOR = -1.0e10;
const double Pi = 3.14159265358979323846;
const double BETA_MAX = 1.0;
const int MAX_TMCMC_ITER = 1e3;

typedef std::vector<double>               RealVector;
typedef std::vector<std::vector<double> > RealMatrix;
typedef std::vector<int>                  IntVector;

RealVector flatten(const RealMatrix& v);

void readPriorParams(double *alphas, double *betas, char *types, int ndim, std::string priorParamsFile, std::string priorTypesFile);

void genPriorSamps(std::default_random_engine &generator,
  std::uniform_real_distribution<double> &u_distribution,
  std::normal_distribution<double> &n_distribution,
  RealMatrix &spls,
  int nsamps, int ndim,
  double *alphas, double *betas, char *types);

void outProcToFile(const RealVector spls, const int ndim, const int nspl, std::string fname);

void shuffle_spls(RealVector &spls, RealVector &llik, RealVector &lprior, RealVector &lproposal, std::default_random_engine &generator);
void cholesky(RealVector &A, int n);

/*************************************************
Will be redefined for likelihood PDF
*************************************************/
class LogDensity{
public:
  LogDensity(){};
  ~LogDensity(){};
  virtual double eval(RealVector&){return 0.0;};
};

/*************************************************
Prior PDF (Uniform or Gaussian)
*************************************************/
class LogPrior{
public:
  LogPrior(){};
  ~LogPrior(){};

  static double eval(RealVector& x, double* alphas, double* betas, char* types){

    int ndim = x.size();
    double log_prior = 0.0;

    for (int j = 0; j < ndim; j++) {
      switch (types[j]){
        case 'U':

          // std::cout << betas[j] << std::endl;
          if (x[j] < alphas[j] || x[j] > betas[j])
            return SMALL_LPRIOR;
          else
            log_prior += log(1.0/(betas[j]-alphas[j]));
          break;

        case 'G':
          log_prior += -0.5*log(2.0*Pi) - log(betas[j]) - 0.5*pow((x[j]-alphas[j])/betas[j],2.0);
          break;
        case 'J':
          log_prior += log(sqrt(2/pow(x[j],2)));
          break;
      }
    }
    return log_prior;
  }
};

double tmcmc(int nspl,
          int iseed, int ndim, double cv,
          int write_flag, int pid, int np,
          std::string priorParamsFile, std::string priorPTypesFile,
          std::string proposalParamsFile, std::string proposalPTypesFile,
          LogDensity* likelihood);

#endif