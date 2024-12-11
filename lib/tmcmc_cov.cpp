#include "tmcmc_cov.h"

double tmcmc(int nsamps,
          int iseed, int ndim, double cv,
          int write_flag, int pid, int np,
          std::string priorParamsFile, std::string priorPTypesFile,
          std::string proposalMeanFile,std::string proposalInvCovFile,
          std::string proposalCovSqrtFile,
          LogDensity* logLikelihood) {
  /* TMCMC Algorithm
      Input: nsamps - number of Samples
             iseed - random seed
             ndim - dimensionality
             cv - Coefficient of Variance threshold for adapting Beta
             write_flag - provide intermediate TMCMC artifacts
             pid - MPI process ID
             np - number of MPI processes
             priorParamsFile - file containing prior types
             priorPTypesFile - file containing prior hyperparameters
             logLikelihood - likelihood instance
      Output: evid - asymptotically unbiased model evidence estimator
  */

  int MFactor = 1;
  int nscal = ndim*nsamps;
  int CATSteps = 1;
  int nspl = nsamps*np;
  int nSteps;
  int nsplSt;

  RealMatrix spls;
  RealVector spls_flat;
  RealVector lprior;
  RealVector lproposal;
  RealVector llik;
  RealVector phi;
  llik.resize(nsamps);
  lprior.resize(nsamps);
  lproposal.resize(nsamps);
  spls_flat.resize(nsamps*ndim);
  phi.resize(ndim);

  RealVector all_spls;
  RealVector all_lprior;
  RealVector all_llik;
  RealVector all_lproposal;
  RealVector all_lprod;
  all_spls.resize(nsamps*np*ndim);
  all_lprior.resize(nspl);
  all_lproposal.resize(nspl);
  all_llik.resize(nspl);
  all_lprod.resize(nspl);

  RealVector splSt, llikSt, lpriorSt, lproposalSt;
  RealVector local_splSt;
  RealVector cvmat(ndim*ndim,0.0);
  RealVector splsComp;
  RealVector splSave, llikSave, lpriorSave, lproposalSave;
  IntVector splCard;

  /* initialize random number generator */
  std::default_random_engine generator(iseed);
  std::uniform_real_distribution<double> u_distribution(0.0,1.0);
  std::normal_distribution<double> n_distribution(0.0,1.0);

  // Roberts and Rosenthal 2011 - Initial gamma
  double gm = 2.38 / sqrt(ndim);
  double gm2 = gm * gm;

  


  // *********** Modified by Han, 09/17/2020 *************************
  // Proposal mean and covariace matrix
  // covMat^-1, det(covMat), Cholesky decomposition of covMat
  RealVector proposal_mean;
  proposal_mean.resize(ndim);
  double proposal_cov_det;
  RealMatrix proposal_cov;
  RealMatrix proposal_cov_inv; 
  
  RealMatrix proposal_C; // Cholesky decomposition of the covariance matrix
  proposal_cov.resize(ndim);
  proposal_cov_inv.resize(ndim);
  proposal_C.resize(ndim);
  for(int i = 0; i < ndim; i++){
  	proposal_cov[i].resize(ndim);
  	proposal_cov_inv[i].resize(ndim);
  	proposal_C[i].resize(ndim);
  }
  


  if (pid == 0){
    // Read inverse covariance matrix, det(cov), and cholesky decomposition of the covariance matrix
  	readProposalFiles(proposal_mean, proposal_cov_inv, proposal_C, proposal_cov_det,
  		ndim, proposalMeanFile, proposalInvCovFile, proposalCovSqrtFile);
  	//proposal_cov_det = Det(proposal_cov, ndim);
  	// matrixInversion(proposal_cov_inv, proposal_cov, ndim, proposal_cov_det);
  	// cholesky(proposal_cov, proposal_C, ndim);
  }
  MPI_Bcast(&proposal_mean[0], ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&proposal_cov_det, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for(int i = 0; i < ndim; i++){
  	MPI_Bcast(&proposal_cov[i][0], ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  	MPI_Bcast(&proposal_cov_inv[i][0], ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  	MPI_Bcast(&proposal_C[i][0], ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  // *********************************************************************************

  //===========================================================================
  // Read prior parameters
  double* alphas = new double[ndim];
  double* betas = new double[ndim];
  char* types = new char[ndim];

  if (pid == 0){
    readPriorParams(alphas, betas, types, ndim, priorParamsFile, priorPTypesFile);
  }
  MPI_Bcast(alphas, ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(betas, ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(types, ndim, MPI_CHAR, 0, MPI_COMM_WORLD);




  
  //********************* Modified by Han ***********************
  // Generate proposal samples, local to each process
  // Generate samples from correlated Gaussian 
  genPriorSamps(generator, n_distribution, spls, nsamps, ndim, 
  				proposal_mean, proposal_C);

  spls_flat = flatten(spls);

  // Compute proposal and prior PDF values for impoprtance weights
  for (int i = 0; i < nsamps; i++) {
    //lproposal[i] = LogPrior::eval(spls[i],proposal_alphas,proposal_betas,proposal_types);
    lproposal[i] = LogPrior::evalProposal(spls[i], proposal_mean, proposal_cov_inv, proposal_cov_det);
    lprior[i] = LogPrior::evalPrior(spls[i],alphas,betas,types);
    llik[i] = logLikelihood->eval(spls[i]);
  }

  MPI_Gather(&spls_flat[0], nsamps*ndim, MPI_DOUBLE,
    &all_spls[0], nsamps*ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&lprior[0], nsamps, MPI_DOUBLE,
    &all_lprior[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&lproposal[0], nsamps, MPI_DOUBLE,
    &all_lproposal[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&llik[0], nsamps, MPI_DOUBLE,
    &all_llik[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Output first set of samples to file*/
  if (pid == 0){
    outProcToFile(all_spls,ndim,nspl,std::string("samples.dat.0"));
    outProcToFile(all_llik,1,nspl,std::string("loglik.dat.0") );
    outProcToFile(all_lprior,1,nspl,std::string("logprior.dat.0"));
    outProcToFile(all_lproposal,1,nspl,std::string("logproposal.dat.0"));
  }

  RealVector Sm;
  double accRatio = 1.0;
  int iter = 0;

  double beta = 0.0, dBeta = 0.0, dBetaFinal, evid = 0.0;
  double max_llik;

  do { // Start algorithm
    iter++;

    if (pid == 0){
      /* shuffle samples */
      shuffle_spls(all_spls,all_llik,all_lprior,all_lproposal,generator);

      /* compute weights */
      RealVector w(nspl,0.0);
      double wsum, wmean, w2mean, wstd;

      dBeta = std::min(BETA_MAX,1.0-beta);

      // compute prior(theta)*likelihood(theta)/prposal(theta)
      for (int j=0; j < nspl; j++)
        all_lprod[j] = all_lprior[j] + all_llik[j] - all_lproposal[j];
      /* used to normalize log_prod values for weight computations */
      max_llik = *std::max_element(all_lprod.begin(), all_lprod.end());

      /* Adapt delta beta as needed */
      do {
        dBetaFinal = dBeta;

        for (int j=0; j < nspl; j++)
          w[j] = exp(dBeta*(all_lprod[j]-max_llik));

        wsum   = std::accumulate(w.begin(), w.end(), 0.0);
        wmean  = wsum / w.size();
        w2mean = std::inner_product(w.begin(), w.end(), w.begin(), 0.0)/ w.size();
        wstd   = sqrt(w2mean- pow(wmean, 2));

        if (wstd/wmean > (cv + 1.0) || wstd == 0) dBeta *= 0.9;
        else if (wstd/wmean > (cv + 0.5) || wstd == 0) dBeta *= 0.95;
        else if (wstd/wmean > (cv + 0.05) || wstd == 0) dBeta *= 0.99;
        else if (wstd/wmean > (cv + 0.005) || wstd == 0) dBeta *= 0.999;
        else if (wstd/wmean > (cv + 0.0005) || wstd == 0) dBeta *= 0.9999;
        else if (wstd/wmean > (cv + 0.00005) || wstd == 0) dBeta *= 0.99999;

        if (dBeta < 1.0e-10)
          break;

      } while (wstd/wmean > (cv + 0.00005) || wstd == 0);
      dBeta = dBetaFinal;
      if (write_flag == 1){
        std::cout<<"DBeta: " << dBeta<<" Wmean: "<<wmean;
        std::cout<<" Wstd: "<<wstd<<" Cv: "<<wstd/wmean<<std::endl;
      }

      beta += dBeta;
      evid += log(wsum) + dBeta*max_llik - log(w.size());
      if (write_flag == 1){
        std::cout<<"Iteration "<<iter<<" Beta= "<<beta;
        std::cout<<" wMean= "<< wmean << " Evid=   " << evid<<std::endl<<std::flush;
      }

      /* Save mean ll for Bayes factor */
      Sm.push_back(wmean);

      /* rescale w and do cumulative sum */
      RealVector wb(nspl);
      std::transform(w.begin(),w.end(),w.begin(),
        std::bind(std::multiplies<double>(), 1.0/wsum, std::placeholders::_1));


      std::partial_sum(&w[0],&w[0]+nspl,&wb[0]);
      wb.insert ( wb.begin() , 0.0 );

      /* Covariance matrix */
      RealVector theta0(ndim);
      for (int i=0; i<ndim; i++ ) {
        theta0[i] = 0.0;
        for (int j=0; j<nspl; j++ ) {
          theta0[i] += all_spls[j*ndim+i]*w[j];
        }
      }

      std::fill(cvmat.begin(), cvmat.end(), 0.0);
      for (int j=0; j < nspl; j++) {
        for (int i1=0; i1<ndim; i1++) {
          for (int i2=0; i2<ndim; i2++) {
            double outp=w[j] * (all_spls[j*ndim+i1]-theta0[i1]) *
                              (all_spls[j*ndim+i2]-theta0[i2]);
            cvmat[i1*ndim+i2] += outp;
            cvmat[i2*ndim+i1] += outp;
    	    }
        }
      }

      /* Control Parameter, Covariance rescaling */
      std::transform(cvmat.begin(), cvmat.end(), cvmat.begin(),
        std::bind(std::multiplies<double>(), gm2, std::placeholders::_1));

      /* Cholesky factorization of the proposal covariance, in-place */
      int chol_info=0;
      char lu='L';
      cholesky(cvmat, ndim);

      /* generate random samples into [0,1] */
      RealVector spl01(nspl);
      for (int j=0; j<nspl; j++)
        spl01[j] = u_distribution(generator);

      /* get bin IDs and count */
      IntVector pos(nspl,0);
      for (int j=0; j < nspl; j++) {
        pos[j] = std::lower_bound(wb.begin(),wb.end(),spl01[j])-wb.begin()-1;
      }

      /* Count number of times a sample is "picked" by PRNG */
      IntVector splPos, splCount;
      for (int j=0; j < nspl; j++) {
        int icount=0;
        for (int ispl=0; ispl<nspl; ispl++) {
          if (pos[ispl]==j) icount += MFactor;
        }
        if (icount>0) {
          splPos.push_back(j);
          splCount.push_back(icount);
        }
      }

      /* Initialize samples that were retained, cardinality, and
      likelihood values */
      splCard.clear();
      nsplSt = splPos.size();

      splSt.clear();
      llikSt.clear();
      lpriorSt.clear();
      lproposalSt.clear();
      /* Resampling Step */
      for (int ispl=0; ispl<nsplSt; ispl++) {
        int isplCount = splCount[ispl];
        for (size_t i = 0; i < isplCount; ++i) {
          for (int j = 0; j < ndim; ++j) {
            splSt.push_back(all_spls[splPos[ispl]*ndim + j]);
          }

          splCard.push_back(CATSteps);
          llikSt.push_back(all_llik[splPos[ispl]]);
          lpriorSt.push_back(all_lprior[splPos[ispl]]);
          lproposalSt.push_back(all_lproposal[splPos[ispl]]);
        }
      }

      nSteps = *std::max_element(splCard.begin(), splCard.end());
    }

    MPI_Bcast(&nSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&evid, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Run single steps of the Markov chains at a time, for the chains
        that need several jumps */
    for (int isbSteps=0; isbSteps < nSteps; isbSteps++ ) {
      if (pid == 0){
        nsplSt = nspl; // Post resampling, nspl # of chains

        RealVector splCand(nsplSt*ndim);
        for (int ispl=0; ispl<nsplSt; ispl++) {
          /* generate candidate */
          RealVector xi(ndim);
          for (size_t i=0; i < ndim; i++)  {
            xi[i] = n_distribution(generator);
          }

    	    for (size_t i=0; i < ndim; i++) {
            splCand[ispl*ndim+i] = splSt[ispl*ndim+i];
            double Lnrv=0.0;

            for (size_t j=0; j < (i+1); ++j) {
                Lnrv += cvmat[j*ndim+i] * xi[j];
            }

            splCand[ispl*ndim+i] += Lnrv;
          } /* done generating candidate */
        }


        /* Compute new likelihoods */
        splsComp.clear();
        // int compCount=0;
        for (int ispl=0; ispl<nsplSt; ispl++) {
          for (int i=0; i<ndim; i++ ) {
            splsComp.push_back(splCand[ispl*ndim+i]);
          }
          // compCount++;
        }
      }

      MPI_Scatter(&splsComp[0], nsamps*ndim, MPI_DOUBLE,
          &spls_flat[0], nsamps*ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      

      int offset = 0;
      for (int i = 0; i < nsamps; i++) {
        for (int j=0; j<ndim; j++ ){
          phi[j] = spls_flat[offset+j];
        }
        offset += ndim;
        llik[i] = logLikelihood->eval(phi);
        lprior[i] = LogPrior::evalPrior(phi,alphas,betas,types);
        lproposal[i] = LogPrior::evalProposal(phi, proposal_mean, proposal_cov_inv, proposal_cov_det);
        
      }

      MPI_Gather(&lprior[0], nsamps, MPI_DOUBLE,
        &all_lprior[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&llik[0], nsamps, MPI_DOUBLE,
        &all_llik[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&lproposal[0], nsamps, MPI_DOUBLE,
        &all_lproposal[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (pid == 0){

        /* decide who jumps */
        int icomp=0;
        int acceptCount = 0;
        RealVector splNew(nsplSt*ndim), llikNew(nsplSt), lpriorNew(nsplSt), lproposalNew(nsplSt);

        for (int ispl=0; ispl<nsplSt; ispl++) {
          double alpha = u_distribution(generator);
          double AcceptRatio = -1.0;

          AcceptRatio = beta * (all_llik[icomp] + all_lprior[icomp]
            - llikSt[icomp] - lpriorSt[icomp] )
            + (1.0 - beta) * (all_lproposal[icomp] - lproposalSt[icomp]);

          if (log(alpha) < AcceptRatio) { // Accept proposal
            for (int i=0; i < ndim; i++) {
              splNew[ispl*ndim+i] = splsComp[icomp*ndim+i];
            }
            lpriorNew[ispl] = all_lprior[icomp];
            llikNew[ispl] = all_llik[icomp];
            lproposalNew[ispl] = all_lproposal[icomp];
            acceptCount++;

          } else { // Reject Proposal
            for (int i=0; i<ndim; i++ ) {
              splNew[ispl*ndim+i] = splSt[ispl*ndim+i] ;
            }
            lpriorNew[ispl] = lpriorSt[ispl];
            llikNew[ispl] = llikSt[ispl];
            lproposalNew[ispl] = lproposalSt[ispl];
          }

          icomp++;
        }


        // /* Clear Proposals for next iteration */
        splSt = splNew;
        llikSt = llikNew;
        lpriorSt = lpriorNew;
        lproposalSt = lproposalNew;
        

        /* Reduce length of chains remaining in samples */
        for (int ispl=0; ispl<splCard.size();) {
          if (splCard[ispl]==0)
            splCard.erase(splCard.begin()+ispl) ;
          else
            ispl++;
        }

        accRatio = (double) acceptCount / (double) nsplSt;

        nsplSt = llikSt.size();
        assert(splCard.size()==llikSt.size());
        assert(splSt.size()==llikSt.size()*ndim);
      }
    }

    if (pid == 0){
      // Rescaling based on Catanach (Thesis 2017)
      // Assumptions made on optimal acceptance rate
      double G = 2.1;
      gm = gm * exp(G * (accRatio - 0.234));
      gm2 = pow(gm, 2.0);

      /* Set samples for next temperature iteration */
      all_spls = splSt;
      all_llik = llikSt;
      all_lprior = lpriorSt;
      all_lproposal = lproposalSt;
      assert(all_llik.size()==nspl);
      assert(all_lprior.size()==nspl);
      assert(all_lproposal.size()==nspl);
      assert(all_spls.size()==nspl*ndim);
        

      outProcToFile(splSt,ndim,nspl,
                      std::string("samples.dat.")+std::to_string(iter));
      outProcToFile(llikSt,1,   nspl,
                      std::string("loglik.dat." )+std::to_string(iter));
      outProcToFile(lpriorSt,1, nspl,
                      std::string("logprior.dat.") + std::to_string(iter));
      outProcToFile(lproposalSt,1, nspl,
                      std::string("logproposal.dat.") + std::to_string(iter));
    }

    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  } while ( ( iter < MAX_TMCMC_ITER ) && (beta<1-1.e-10) );

  if (pid == 0){
    if (write_flag == 1){
      std::cout << "TMCMC Algorithm Done" << std::endl;
    }
  }
  return (evid);
} /* done tmcmc */

void cholesky(RealVector &A, int n) {
  RealVector L;
  L.resize(n*n);

  for (int i = 0; i < n; i++)
      for (int j = 0; j < (i+1); j++) {
          double s = 0;
          for (int k = 0; k < j; k++)
              s += L[i * n + k] * L[j * n + k];
          L[i * n + j] = (i == j) ?
                         sqrt(A[i * n + i] - s) :
                         (1.0 / L[j * n + j] * (A[i * n + j] - s));
      }

  for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
          A[i * n + j] = L[i * n + j];
      }
}

// *********** Added by Han *******************
// Cholesky that keeps the original matrix
void cholesky(RealMatrix &A, RealMatrix &C, int n) {
	//RealMatrix L;
	C.resize(n);
	for (int i = 0; i < n; i++){
		C[i].resize(n);
	  for (int j = 0; j < (i+1); j++) {
	      double s = 0;
	      for (int k = 0; k < j; k++)
	          s += C[i][k] * C[j][k];
	      C[i][j] = (i == j) ?
	                     sqrt(A[i][i] - s) :
	                     (1.0 / C[j][j] * (A[i][j] - s));
	  }
    }
}

RealVector flatten(const RealMatrix& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    RealVector result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}



void readPriorParams(double *alphas, double *betas, char *types, int ndim, std::string priorParamsFile, std::string priorTypesFile) {
  std::string line;
  std::ifstream DAT;
  double temp;
  std::string::size_type sz;
  int i = 0;

  // Prior types
  DAT.open(priorTypesFile);
  while (std::getline(DAT, line)) {
      types[i] = line.at(0);
      i = i+1;
  }
  DAT.close();

  // Prior parameters
  i = 0;
  DAT.open(priorParamsFile);
  while (std::getline(DAT, line)) {
      try {
        temp = std::stod (line,&sz);
      } catch (const std::exception& ex) {
        break;
      }
      line = line.substr(sz);
      alphas[i] = temp;
      try {
        temp = std::stod (line,&sz);
      } catch (const std::exception& ex) {
        break;
      }
      line = line.substr(sz);
      betas[i] = temp;
      i = i+1;
  }
  DAT.close();
}

void genPriorSamps(std::default_random_engine &generator,
  std::uniform_real_distribution<double> &u_distribution,
  std::normal_distribution<double> &n_distribution,
  RealMatrix &spls,
  int nsamps, int ndim,
  double *alphas, double *betas, char *types) {

  spls.resize(nsamps);

  for (int i = 0; i < nsamps; i++) {
    spls[i].resize(ndim);

    for (int j = 0; j < ndim; j++) {
      switch (types[j]){
        case 'U':
          spls[i][j] = alphas[j] + (betas[j]-alphas[j])*u_distribution(generator);
          break;

        case 'G':
          spls[i][j] =  alphas[j] + betas[j]*n_distribution(generator);
          break;
      }
    }
  }
}

//************************ Added by Han *******************
// Overloaded function, generate correlated Gaussian samples
void genPriorSamps(std::default_random_engine &generator,
  std::normal_distribution<double> &n_distribution,
  RealMatrix &spls, int nsamps, int ndim,
  RealVector &mean, RealMatrix &C) {

  spls.resize(nsamps);
  RealVector temp_results;
  temp_results.resize(ndim);
  for (int i = 0; i < nsamps; i++) {
    spls[i].resize(ndim);
  	for (int j = 0; j < ndim; j++) 
  		spls[i][j] =  n_distribution(generator);	
  	// spls[i] = mu + C*spls[i], where CC^T = covarianceMatrix 
  	temp_results = matVec(C, spls[i]);
  	vectorAdd(mean, temp_results, spls[i]);
  }


}

void outProcToFile(const RealVector spls, const int ndim, const int
nspl, std::string fname) {

  FILE *myfile  = fopen(fname.c_str(),"w") ;
  for (int j = 0; j < nspl; j++) {
    for (int i = 0; i < ndim; i++)
      fprintf(myfile,"%24.18e ",spls[j*ndim+i]);
    fprintf(myfile,"\n");
  }
  fclose(myfile);

  return ;

}

void shuffle_spls(RealVector &spls, RealVector &llik, RealVector &lprior, RealVector &lproposal, std::default_random_engine &generator) {
  // Shuffle the samples randomly
  int nspl = llik.size();
  int ndim = spls.size()/nspl;

  IntVector idx(nspl);
  for (int j=0; j<nspl; j++) idx[j]=j;

  RealVector splsTmp(nspl*ndim), llikTmp(nspl), lpriorTmp(nspl), lproposalTmp(nspl);
  shuffle (idx.begin(), idx.end(), generator);
  for (int j = 0; j < nspl; j++) {
    llikTmp[j] = llik[idx[j]];
    lpriorTmp[j] = lprior[idx[j]];
    lproposalTmp[j] = lproposal[idx[j]];
    for (int i = 0; i < ndim; i++) splsTmp[j*ndim+i] = spls[idx[j]*ndim+i];
  }

  for (int j = 0; j < nspl; j++) {
    llik[j] = llikTmp[j];
    lprior[j] = lpriorTmp[j];
    lproposal[j] = lproposalTmp[j];
    for (int i = 0; i < ndim; i++) spls[j*ndim+i] = splsTmp[j*ndim+i];
  }

  return ;

}

// ******************** 09/17/2020 Added by Han **********************
// Read proposal mean, inverse matrix, determinant 
// and cholesky decomposition of the covariance matrix 
void readProposalFiles(RealVector& proposal_mean, RealMatrix& proposal_cov_inv, 
					RealMatrix& proposal_C, double& proposal_cov_det, int ndim, 
					std::string proposalMeanFile, std::string proposalInvCovFile, 
					std::string proposalCovSqrtFile){
	std::ifstream mean(proposalMeanFile);
	std::ifstream cov_inv(proposalInvCovFile);	
	std::ifstream cov_sqrt(proposalCovSqrtFile);
	proposal_mean.resize(ndim);
	proposal_cov_inv.resize(ndim);
	cov_inv >> proposal_cov_det;
	for(int i = 0; i < ndim; i++){
		proposal_cov_inv[i].resize(ndim);
		mean >> proposal_mean[i];
		for(int j = 0; j < ndim; j++){
			cov_inv >> proposal_cov_inv[i][j];
			cov_sqrt >> proposal_C[i][j];
		}
		
	}
	mean.close();
	cov_inv.close();
	cov_sqrt.close();
}