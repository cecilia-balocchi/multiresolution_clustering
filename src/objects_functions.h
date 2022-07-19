#ifndef OBJECTS_FUNCTIONS_H_
#define OBJECTS_FUNCTIONS_H_

#include <stdio.h>
#include <armadillo>
#include <random>
#include "partition.h"
#include "multires.h"

using namespace std;
using namespace arma;

arma::vec getloglike(arma::vec grid_eta, LPPartition part, double temp);
// arma::vec getlogpost_etaLR(arma::vec grid_eta, arma::vec prior_eta, int nchains, std::vector<LPMultires> chains, double Ts[]);
// arma::vec getlogpost_etaCT(arma::vec grid_eta, arma::vec prior_eta, int nchains, std::vector<LPMultires> chains, double Ts[]);
// arma::vec getlogpost_etaTD(arma::vec grid_eta, arma::vec prior_eta, int nchains, std::vector<LPMultires> chains, double Ts[]);
arma::vec getlogpost_etaLR(arma::vec grid_eta, arma::vec prior_eta, LPMultires chain, double Temp);
arma::vec getlogpost_etaLR(arma::vec grid_eta, arma::vec prior_eta, LPPartition chain, double Temp);
arma::vec getlogpost_etaCT(arma::vec grid_eta, arma::vec prior_eta, LPMultires chain, double Temp);
arma::vec getlogpost_etaTD(arma::vec grid_eta, arma::vec prior_eta, LPMultires chain, double Temp);
#endif /* OBJECTS_FUNCTIONS_H_ */