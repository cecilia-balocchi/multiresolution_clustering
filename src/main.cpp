#include <armadillo>
#include <iostream>
#include <vector>
#include <set>
#include "various_functions.h"
#include "objects_functions.h"
#include "partition.h"
#include "hdp.h"
#include "multires.h"
using namespace std;
using namespace arma;

arma::mat Y;
arma::mat YLR;
arma::imat Mapping;
arma::imat InvMapping;
// double eta_LR;
// double eta_CT;
// double eta_TD;
double mu_H = 0;
double mu_L = 0;
double k_H = 0.001;
double k_L = 0.001;
double alpha_H = 1;
double beta_H = 1;
double alpha_L = 1;
double beta_L = 1;
bool print_bool = false;
bool yLR_bool = false;
bool sampleLR_bool = true;
bool sampleHR_bool = true;
double alpha1 = 1.0;
double alpha2 = 0.0;


extern "C" {

void mainfunc(double Y_input[], double YLR_input[], 
              double Mapping_input[], double InvMapping_input[], int* n_ptr,
              int* iter_burnin_ptr, int* iter_adjust_ptr, int* iter_save_ptr, int* gibbs_iter_ptr,
              int* print_bool_ptr, int* sampleLR_bool_ptr, int* sampleHR_bool_ptr, 
              int* one_vs_many_LR_ptr, int* one_vs_many_HR_ptr,
              double* eta_LR_ptr, double* eta_CT_ptr, double* eta_TD_ptr, 
              int* ngrid_eta_ptr, double grid_eta_input[], double prior_eta_input[],
              int* sampleEtaLR_bool_ptr, int* sampleEtaCT_bool_ptr, int* sampleEtaTD_bool_ptr,
              double* k_H_ptr, double* alpha_H_ptr, double* beta_H_ptr, 
              double* k_L_ptr, double* alpha_L_ptr, double* beta_L_ptr, 
              double* alpha1_ptr, double* alpha2_ptr, 
              int seed_ptr[], int gammaLR_init_ptr[], int gammaHR_init_ptr[], int* KLR_init_ptr, int* KHR_init_ptr,
              double accept[], double MCMCchain_low[], double MCMCchain_high[], 
              double MCMCchain_highRest[], double MCMCchain_highTable[], //double MCMCchain_highDish[], 
              double MCMCchain_eta[], double switch_chain[],
              double postmean_low[], double postmean_high[],
              int* nchains_ptr, double Ts[]){ 

  // Initial setup
  int n = *n_ptr;
  Y.set_size(n,1);
  Mapping.set_size(n,2);
  InvMapping.set_size(n,2);
  for(int i = 0; i < n; i++){
    Y(i) = Y_input[i];
    Mapping(i,0) = Mapping_input[i];
    Mapping(i,1) = Mapping_input[n+i];
    InvMapping(i,0) = InvMapping_input[i];
    InvMapping(i,1) = InvMapping_input[n+i];
  }
  int ngrid_eta = *ngrid_eta_ptr;
  arma::vec grid_eta, prior_eta;
  grid_eta.set_size(ngrid_eta);
  prior_eta.set_size(ngrid_eta);
  for(int i = 0; i < ngrid_eta; i++){
    grid_eta(i) = grid_eta_input[i];
    prior_eta(i) = prior_eta_input[i];
  }
  bool sample_eta = true;
  if(ngrid_eta == 1){
    sample_eta = false;
  }

  // Setup of restaurants (or ct)
  std::set<int> rests;
  for(int i = 0; i < Mapping.n_rows; i++){
    int r = Mapping(i,1);
    rests.insert(r);
  }
  int nrest = rests.size();

  if(YLR_input[0] != 0){
    yLR_bool = true;
  }
  if(yLR_bool){
    YLR.set_size(nrest,1);
    for(int i = 0; i < nrest; i++){
      YLR(i) = YLR_input[i];
    }
    alpha1 = alpha1_ptr[0];
    alpha2 = 1 - alpha1;
  } else {
    alpha1 = 1.0;
    alpha2 = 0.0;
  }
  if(alpha2_ptr[0] != 0){
    alpha2 = alpha2_ptr[0];
  }
  
  // eta_LR = *eta_LR_ptr;
  // eta_CT = *eta_CT_ptr;
  // eta_TD = *eta_TD_ptr;
  k_H = *k_H_ptr;
  alpha_H = *alpha_H_ptr;
  beta_H = *beta_H_ptr;
  k_L = *k_L_ptr;
  alpha_L = *alpha_L_ptr;
  beta_L = *beta_L_ptr;
  int iter_burnin = *iter_burnin_ptr;
  int iter_adjust = *iter_adjust_ptr;
  int iter_save = *iter_save_ptr;
  int gibbs_iter = *gibbs_iter_ptr;
  int nchains = *nchains_ptr;

  std::vector<double> etas_LR(nchains);
  std::vector<double> etas_CT(nchains);
  std::vector<double> etas_TD(nchains);
  for(int c = 0; c < nchains; c++){
    etas_LR[c] = *eta_LR_ptr;
    etas_CT[c] = *eta_CT_ptr;
    etas_TD[c] = *eta_TD_ptr;
  }

  bool print_bool, sampleLR_bool, sampleHR_bool;
  if(*print_bool_ptr == 1){
    print_bool = true;
  } else {
    print_bool = false;
  }
  if(*sampleLR_bool_ptr == 1){
    sampleLR_bool = true;
  } else {
    sampleLR_bool = false;
  }
  if(*sampleHR_bool_ptr == 1){
    sampleHR_bool = true;
  } else {
    sampleHR_bool = false;
  }
  bool sample_etaLR, sample_etaCT, sample_etaTD;
  if(*sampleEtaLR_bool_ptr == 1){
    sample_etaLR = true;
  } else {
    sample_etaLR = false;
  }
  if(*sampleEtaCT_bool_ptr == 1){
    sample_etaCT = true;
  } else {
    sample_etaCT = false;
  }
  if(*sampleEtaTD_bool_ptr == 1){
    sample_etaTD = true;
  } else {
    sample_etaTD = false;
  }

  unsigned seed;
  if(seed_ptr[0] != 0){
    seed = seed_ptr[0];
  } else {
    seed = time(0);
    seed_ptr[0] = (int)seed;
  }
  std::default_random_engine eng(seed);
  cout << "Random seed: " << seed << endl;
  
  std::vector<int> rest_config(nrest); // tells me how many bg in each ct (how many costumers in each restaurant)
  std::vector< std::vector<int> > Rests(nrest); // tells me which bg is in each ct (which costumer in each restaurant)
  for(int i = 0; i < Mapping.n_rows; i++){
    rest_config[Mapping(i,1)]++;
    Rests[Mapping(i,1)].push_back(Mapping(i,0));
  }

  LPMultires point = new Multires();
  if(gammaHR_init_ptr[0] == 0){
    if(gammaLR_init_ptr[0] == 0){
      // If none are specified, let's use the bool
      if(*one_vs_many_LR_ptr == 1){ // 1 cluster in LR (one restaurant)
        if(*one_vs_many_HR_ptr == 1){ // 1 cluster in HR
          point->Initialize_Multires(nrest, n, Rests, etas_LR[0], etas_CT[0], etas_TD[0]); // Initialize_Multires sets everything to one cluster
          // SOFT START: Since one cluster everyone might get the MCMC stuck
          double lqratio, lpratio, llikratio;
          LPHdp point2 = new Hdp(point->highresPart);
          for(int t = 0; t < 9; t++){
            point->highresPart->ProposalSM_CostTable(point2, 0, gibbs_iter, lqratio, lpratio, llikratio, eng, etas_CT[0], etas_TD[0]); // point2 not changed
            delete point2;
            point2 = new Hdp(point->highresPart);
          }
          delete point2;
        } else {
          point->Initialize_Multires_oneLR_nclHR(nrest, n, Rests, etas_LR[0], etas_CT[0], etas_TD[0]);
        }
      } else { // then initialize n clusters in LR (many restaurant)
        if(*one_vs_many_HR_ptr == 1){ 
          point->Initialize_Multires_nclLR_oneHR(nrest, n, Rests, etas_LR[0], etas_CT[0], etas_TD[0]);
        } else {
          point->Initialize_Multires_nclLR_nclHR(nrest, n, Rests, etas_LR[0], etas_CT[0], etas_TD[0]);
        }
      }
    } else {
      // here we have gammaLR_init but we use bool for HR
      if(*one_vs_many_HR_ptr == 1){
        point->Initialize_Multires_partLR_oneHR(nrest, n, Rests, gammaLR_init_ptr, KLR_init_ptr, etas_LR[0], etas_CT[0], etas_TD[0]);
        // "soft-start": we need to split the costumer partition into various clusters to avoid staying stuck. 
        double lqratio, lpratio, llikratio;
        LPHdp point2 = new Hdp(point->highresPart);
        for(int i = 0; i < *KLR_init_ptr; i++){
          point->highresPart->ProposalSM_CostTable(point2, i, gibbs_iter, lqratio, lpratio, llikratio, eng, etas_CT[0], etas_TD[0]); // point2 not changed
          delete point2;
          point2 = new Hdp(point->highresPart);
        }
        delete point2;
      } else {
        point->Initialize_Multires_partLR_nclHR(nrest, n, Rests, gammaLR_init_ptr, KLR_init_ptr, etas_LR[0], etas_CT[0], etas_TD[0]);
      }
    }
  } else {
    
    if(gammaLR_init_ptr[0] == 0){
      if(*one_vs_many_LR_ptr == 1){
        point->Initialize_Multires_oneLR_partHR(nrest, n, Rests, gammaHR_init_ptr, KHR_init_ptr, etas_LR[0], etas_CT[0], etas_TD[0]);
        // double lqratio, lqratio2, lpratio, llikratio;
        // for(int t = 0; t < 9; t++){
        //   point->ProposalSM_LR(gibbs_iter, lqratio, lqratio2, lpratio, llikratio, eng); 
        // }
      } else {
        point->Initialize_Multires_nclLR_partHR(nrest, n, Rests, gammaHR_init_ptr, etas_LR[0], etas_CT[0], etas_TD[0]);
      }
    } else {
      // cout << "ERROR (partLR_partHR), you did not code this function yet!!" << endl;
      point->Initialize_Multires_partLR_partHR(nrest, n, Rests, gammaLR_init_ptr, KLR_init_ptr, gammaHR_init_ptr, etas_LR[0], etas_CT[0], etas_TD[0]);
    }
  }
  if(print_bool){
    cout << endl << "***** Starting multires object: " << endl;
    point->Print_Multires();
    cout << "*****" << endl << endl; 
  }
  

  // // "soft-start": we need to split the costumer partition into various clusters to avoid staying stuck. Can do this in a smarter way
  // double lqratio, lpratio, llikratio;
  // // With phillysim4 we don't need to split the restaurants
  // LPHdp point2 = new Hdp(point->highresPart);
  // for(int i = 0; i < nrest; i++){
  // // for(int i = 0; i < 1; i++){
  //   point->highresPart->ProposalSM_CostTable(point2, i, gibbs_iter, lqratio, lpratio, llikratio, eng); // point2 not changed
  //   delete point2;
  //   point2 = new Hdp(point->highresPart);
  // }
  // delete point2;

  
  // save the starting partition: when no index, it's 0
  point->lowresPart->format_partition(MCMCchain_low);
  if(yLR_bool){
    point->lowresPart->get_postmean(postmean_low); // lowres = true  
  }
  LPPartition CostDish;
  CostDish = new Partition();
  point->Get_Highres_CostDish(CostDish);
  CostDish->format_partition(MCMCchain_high);
  CostDish->get_postmean(postmean_high,0,false);
  if(print_bool){
    cout << endl << "***** Starting HR partition (costumers-dishes): " << endl;
    CostDish->Print_Partition();
    cout << "*****" << endl << endl;
  }
  delete CostDish; 
  point->highresPart->format_rest_table(MCMCchain_highRest, MCMCchain_highTable);

  MCMCchain_eta[0] = etas_LR[nchains-1];
  MCMCchain_eta[1] = etas_CT[nchains-1];
  MCMCchain_eta[2] = etas_TD[nchains-1];

  // Initialize the vector of chains
  std::vector<LPMultires> chains(nchains);
  for(int c = 0; c < nchains; c++){
    chains[c] = new Multires(point);
  }
  cout << "chains initialized! length: " << nchains << endl;
  if(print_bool){
    for(int c = 0; c < nchains; c++){
      cout << endl << "***** Printing the chains starting points: " << c << endl;
      chains[c]->Print_Multires();
      cout << "*****" << endl << endl; 
    }
  }
  // keep track of acceptance rate (three parameters: 0. lowres, 1. highres: cost_table , 2. highres: table_dishes )
  std::vector<int> accept_int(3*nchains,0);
  std::vector<int> accept0_int(nchains-1,0);
  std::vector<int> count0_int(nchains-1,0);

  std::vector<bool> v;
  
  for(int t = 0; t < iter_burnin; t++){  
    for(int c = 0; c < nchains; c++){
      v = chains[c]->SampleSM(gibbs_iter, Ts[c], eng, etas_LR[c], etas_CT[c], etas_TD[c], sampleLR_bool, sampleHR_bool);
      // v = chains[c]->SampleMix(gibbs_iter, Ts[c], eng, etas_LR[c], etas_CT[c], etas_TD[c], sampleLR_bool, sampleHR_bool);
      if(v[0]) {
        accept_int[c*3 + 0]++;
      }
      if(v[1]){ 
        accept_int[c*3 + 1]++;
      }
      if(v[2]) {
        accept_int[c*3 + 2]++;
      }
    }
  }

  int starting_index_high;
  int starting_index_low;
  int starting_index_eta;
  int ind_switch1, ind_switch2;
  double pr_i1, pr_i2, log_acc_switch, acc_switch, eta_switch;
  LPMultires point_switch;

  for(int t = 0; t < iter_save; t++){
    if(t % 100 == 0){
      cout << t << " ";
    }
    
    for(int c = 0; c < nchains; c++){
      v = chains[c]->SampleSM(gibbs_iter, Ts[c], eng, etas_LR[c], etas_CT[c], etas_TD[c], sampleLR_bool, sampleHR_bool);
      // v = chains[c]->SampleMix(gibbs_iter, Ts[c], eng, etas_LR[c], etas_CT[c], etas_TD[c], sampleLR_bool, sampleHR_bool);
      if(v[0]) {
        accept_int[c*3 + 0]++;
      }
      if(v[1]){ 
        accept_int[c*3 + 1]++;
      }
      if(v[2]) {
        accept_int[c*3 + 2]++;
      }
    }
    
    // here we switch states
    if(nchains > 1){
      if(nchains == 2){
        ind_switch1 = 0;
        ind_switch2 = 1;
      } else {
        std::uniform_int_distribution<int> distribution_int(0, nchains-2);
        ind_switch1 = distribution_int(eng);
        ind_switch2 = ind_switch1+1;
      }

      // cout << "ind_switch " << ind_switch;
      pr_i1 = chains[ind_switch1]->get_logposterior(etas_LR[ind_switch1], etas_CT[ind_switch1], etas_TD[ind_switch1]);
      pr_i2 = chains[ind_switch2]->get_logposterior(etas_LR[ind_switch2], etas_CT[ind_switch2], etas_TD[ind_switch2]);
      log_acc_switch = (Ts[ind_switch2] - Ts[ind_switch1]) *  (pr_i1 - pr_i2);
      acc_switch = exp(log_acc_switch);
      acc_switch = min(1.0,acc_switch);
      // cout << " acc_switch " << acc_switch << endl;
      uniform_real_distribution<double> distribution_unif(0.0,1.0);
      if(distribution_unif(eng) < acc_switch){
        accept0_int[ind_switch1]++;
        count0_int[ind_switch1]++;

        point_switch = new Multires(chains[ind_switch1]);
        delete chains[ind_switch1];
        chains[ind_switch1] = new Multires(chains[ind_switch2]);
        delete chains[ind_switch2];
        chains[ind_switch2] = new Multires(point_switch);
        delete point_switch;

        eta_switch = etas_LR[ind_switch1];
        etas_LR[ind_switch1] = etas_LR[ind_switch2];
        etas_LR[ind_switch2] = eta_switch;

        eta_switch = etas_CT[ind_switch1];
        etas_CT[ind_switch1] = etas_CT[ind_switch2];
        etas_CT[ind_switch2] = eta_switch;

        eta_switch = etas_TD[ind_switch1];
        etas_TD[ind_switch1] = etas_TD[ind_switch2];
        etas_TD[ind_switch2] = eta_switch;

        if(ind_switch1 == nchains-2){
          switch_chain[1 + t] = 1;
        }
      } else {
        count0_int[ind_switch1]++;
        // nothing for now
      }
    }
    
    // let's save the output partition even if sample bool is false
    starting_index_low = (1 + t) * nrest; // 1+ because we save the starting point
    chains[nchains-1]->lowresPart->format_partition(MCMCchain_low, starting_index_low);
    if(yLR_bool){
      chains[nchains-1]->lowresPart->get_postmean(postmean_low, starting_index_low);
    }
    starting_index_high = (1 + t) * n; // 1+ because we save the starting point
    // LPPartition CostDish;
    CostDish = new Partition();
    chains[nchains-1]->Get_Highres_CostDish(CostDish);
    CostDish->format_partition(MCMCchain_high, starting_index_high);
    CostDish->get_postmean(postmean_high, starting_index_high, false);
    delete CostDish;
    chains[nchains-1]->highresPart->format_rest_table(MCMCchain_highRest, MCMCchain_highTable, starting_index_high);
    
    if(sample_eta){
      std::vector<double> pr(ngrid_eta);
      arma::vec pr_vec;
      int i_eta;
      double norm;
      for(int c = 0; c < nchains; c++){
        // etaLR
        if(sample_etaLR){
          // pr_vec = getlogpost_etaLR(grid_eta, prior_eta, nchains, chains, Ts);
          pr_vec = getlogpost_etaLR(grid_eta, prior_eta, chains[c], Ts[c]);
          pr_vec = pr_vec - max(pr_vec);
          norm = sum(exp(pr_vec));
          for(int i = 0; i < ngrid_eta; i++){
            pr[i] = exp(pr_vec(i))/norm;
          }
          std::discrete_distribution<int> distribution1(pr.begin(), pr.end());
          i_eta = distribution1(eng);
          etas_LR[c] = grid_eta[i_eta];
        }
        

        // etaCT
        if(sample_etaCT){
          // pr_vec = getlogpost_etaCT(grid_eta, prior_eta, nchains, chains, Ts);
          pr_vec = getlogpost_etaCT(grid_eta, prior_eta, chains[c], Ts[c]);
          pr_vec = pr_vec - max(pr_vec);
          norm = sum(exp(pr_vec));
          for(int i = 0; i < ngrid_eta; i++){
            pr[i] = exp(pr_vec(i))/norm;
          }
          std::discrete_distribution<int> distribution2(pr.begin(), pr.end());
          i_eta = distribution2(eng);
          etas_CT[c] = grid_eta[i_eta];
        }

        // etaTD
        if(sample_etaTD){
          // pr_vec = getlogpost_etaTD(grid_eta, prior_eta, nchains, chains, Ts);
          pr_vec = getlogpost_etaTD(grid_eta, prior_eta, chains[c], Ts[c]);
          pr_vec = pr_vec - max(pr_vec);
          norm = sum(exp(pr_vec));
          for(int i = 0; i < ngrid_eta; i++){
            pr[i] = exp(pr_vec(i))/norm;
          }
          std::discrete_distribution<int> distribution3(pr.begin(), pr.end());
          i_eta = distribution3(eng);
          etas_TD[c] = grid_eta[i_eta];
        }
      }
      

      starting_index_eta = (1 + t) * 3;
      MCMCchain_eta[starting_index_eta + 0] = etas_LR[nchains-1];
      MCMCchain_eta[starting_index_eta + 1] = etas_CT[nchains-1];
      MCMCchain_eta[starting_index_eta + 2] = etas_TD[nchains-1];
    }
  }
  cout << endl;
  for(int c = 0; c < nchains; c++){
    accept[c*3 + 0] = (double) accept_int[c*3 + 0]/(iter_save);
    accept[c*3 + 1] = (double) accept_int[c*3 + 1]/((iter_save)*(1+iter_adjust));
    accept[c*3 + 2] = (double) accept_int[c*3 + 2]/((iter_save)*(1+iter_adjust));
  }
  for(int c = 0; c < nchains-1; c++){
    double d1 = (double) (accept0_int[c]+1);
    double d2 = (double) (count0_int[c]+1);
    accept[nchains*3 + c] = (double) d1/d2;
    // accept[nchains*3 + c] = (double) (accept0_int[c]+1)/(count0_int[c]+1);
  }

  if(print_bool){
    cout << endl << "******** LAST ITERATION ********" << endl;
    // cout << "--- LowRes ---" << endl;
    // point->lowresPart->Print_Partition();
    // cout << endl;
    // cout << "--- HighRes ---" << endl;
    // for(int i = 0; i < point->highresPart->nRest; i++){
    //   point->highresPart->costumers_tables[i]->Print_Partition();
    //   point->highresPart->costumers_tables[i]->Print_Y(false);
    //   cout << endl;
    // }
    // cout << "-- TableDish --" << endl;
    // point->highresPart->tables_dishes->Print_Partition();
    // cout << endl;
    // point->highresPart->Print_Y();

    chains[nchains-1]->Print_Multires();
    cout << "*****" << endl;  
  }

  // cout << point->nObs_lowres << " " << point->nObs_highres << endl;
  // why are we overwriting MCMCchain?? when index is not specified, it means 0!
  // point->lowresPart->format_partition(MCMCchain_low);
  // LPPartition CostDish;
  // CostDish = new Partition();
  // point->Get_Highres_CostDish(CostDish);
  // CostDish->format_partition(MCMCchain_high);
  for(int c = 0; c < nchains; c++){
    delete chains[c];
  }
  delete point;
  // return 0;
}

}
