#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <random>
#include <algorithm>
#include <armadillo>
#include <unordered_set>
#include "various_functions.h"
#include "partition.h"
#include "multires.h"
using namespace std;
using namespace arma;

extern arma::mat Y;
extern arma::mat X;
// extern arma::mat A_block;
extern arma::imat Mapping;
extern arma::imat InvMapping;

arma::vec getloglike(arma::vec grid_eta, LPPartition part, double temp){
  arma::vec out = part->K * log(grid_eta) - (lgamma(grid_eta + part->nObs) - lgamma(grid_eta));
  out = temp * out;
  return out;
}
// arma::vec getlogpost_etaLR(arma::vec grid_eta, arma::vec prior_eta, int nchains, std::vector<LPMultires> chains, double Ts[]){
//   arma::vec out = log(prior_eta);
//   for(int c = 0; c < nchains; c++){
//     out += getloglike(grid_eta, chains[c]->lowresPart, Ts[c]);
//   }
//   return out;
// }
arma::vec getlogpost_etaLR(arma::vec grid_eta, arma::vec prior_eta, LPMultires chain, double Temp){
  arma::vec out = log(prior_eta);
  out += getloglike(grid_eta, chain->lowresPart, Temp);
  return out;
}
arma::vec getlogpost_etaLR(arma::vec grid_eta, arma::vec prior_eta, LPPartition chain, double Temp){
  arma::vec out = log(prior_eta);
  out += getloglike(grid_eta, chain, Temp);
  return out;
}
// arma::vec getlogpost_etaCT(arma::vec grid_eta, arma::vec prior_eta, int nchains, std::vector<LPMultires> chains, double Ts[]){
//   arma::vec out = log(prior_eta);
//   for(int c = 0; c < nchains; c++){
//     // each highresPart has many CT partitions (one for each rest)
//     for(int r = 0; r < chains[c]->highresPart->nRest; r++){
//       out += getloglike(grid_eta, chains[c]->highresPart->costumers_tables[r], Ts[c]);
//     }
//   }
//   return out;
// }
arma::vec getlogpost_etaCT(arma::vec grid_eta, arma::vec prior_eta, LPMultires chain, double Temp){
  arma::vec out = log(prior_eta);
  for(int r = 0; r < chain->highresPart->nRest; r++){
    out += getloglike(grid_eta, chain->highresPart->costumers_tables[r], Temp);
  }
  return out;
}

// arma::vec getlogpost_etaTD(arma::vec grid_eta, arma::vec prior_eta, int nchains, std::vector<LPMultires> chains, double Ts[]){
//   arma::vec out = log(prior_eta);
//   for(int c = 0; c < nchains; c++){
//     out += getloglike(grid_eta, chains[c]->highresPart->tables_dishes, Ts[c]);
//   }
//   return out;
// }
arma::vec getlogpost_etaTD(arma::vec grid_eta, arma::vec prior_eta, LPMultires chain, double Temp){
  arma::vec out = log(prior_eta);
  out += getloglike(grid_eta, chain->highresPart->tables_dishes, Temp);
  return out;
}