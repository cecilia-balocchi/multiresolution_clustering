// I NEED TO CHECK THAT NOBS IS UPDATED BECAUSE NOW IT CAN CHANGE!! 
#include <armadillo>  
#include <math.h>
#include <iostream>
#include <string.h>
#include <armadillo>
#include <list>
#include <random>
#include "partition.h"
#include "various_functions.h"
#include "objects_functions.h"
using namespace std;
using namespace arma;

extern arma::mat Y;
extern arma::mat YLR;
extern arma::imat Mapping;
extern arma::imat InvMapping;
// extern arma::mat sX;
// Set the hyper-parameters
// extern double rho;
// extern double a1;
// extern double a2;
// extern double b1;
// extern double b2;
// extern double a_cohes;
// extern double b_cohes;
// extern double nu_sigma;
// extern double alpha_sigma;
// extern double eta;
// extern double eta_LR;
// extern double eta_CT;
// extern double eta_TD;
// extern int method;
// extern int* sigma;
extern double mu_H;
extern double mu_L;
extern double k_H;
extern double k_L;
extern double alpha_H;
extern double beta_H;
extern double alpha_L;
extern double beta_L;

// stupid comment

//class Partition::Begins
Partition::Partition(){
  nObs = 0;
  K = 0;
  cluster_config = NULL;
  clusters = NULL;
  cluster_assignment = NULL;
  pairwise_assignment = NULL;
  log_prior = NULL;
  indexing = std::vector<int>();
  // log_det_Omegay = NULL;
  // quad_forms = NULL;
  return;
}
// initializes a partition
Partition::Partition(LPPartition initial_partition){
  nObs = initial_partition->nObs;
  K = initial_partition->K;
  cluster_config = new int[K];
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  pairwise_assignment = new int*[nObs];
  log_prior = new double[K];
  indexing = std::vector<int>(initial_partition->indexing);
  // log_det_Omegay = new double[K];
  // quad_forms = new double[K];

  for(int k = 0; k < K; k++){
    cluster_config[k] = initial_partition->cluster_config[k];
    clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      clusters[k][i] = initial_partition->clusters[k][i];
    }
    log_prior[k] = initial_partition->log_prior[k];
    // log_det_Omegay[k] = initial_partition->log_det_Omegay[k];
    // quad_forms[k] = initial_partition->quad_forms[k];
  }

  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = initial_partition->cluster_assignment[i];
    pairwise_assignment[i] = new int[nObs];
    for(int j = 0; j < nObs; j++){
      pairwise_assignment[i][j] = initial_partition->pairwise_assignment[i][j];
    }
  }
  return;
}

void Partition::Copy_Partition(LPPartition initial_partition){
  delete[] cluster_config;
  for(int k = 0; k < K; k++){
    delete[] clusters[k];
  }
  delete[] clusters;
  delete[] cluster_assignment;
  for(int i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; 
  }
  delete[] pairwise_assignment; 
  delete[] log_prior;

  indexing.resize(0);
  nObs = initial_partition->nObs;
  K = initial_partition->K;

  cluster_config = new int[K];
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  pairwise_assignment = new int*[nObs];
  log_prior = new double[K];
  indexing = std::vector<int>(initial_partition->indexing);
  
  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = initial_partition->cluster_assignment[i];
    pairwise_assignment[i] = new int[nObs];
    for(int j = 0 ; j < nObs; j++){
      pairwise_assignment[i][j] = initial_partition->pairwise_assignment[i][j];
    }
  }
  for(int k = 0; k < K; k++){
    cluster_config[k] = initial_partition->cluster_config[k];
    clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      clusters[k][i] = initial_partition->clusters[k][i];
    }
    log_prior[k] = initial_partition->log_prior[k];
  }
  return;
}

Partition::~Partition(){
  int i;
  delete[] cluster_config; cluster_config = NULL;
  for(i = 0; i < K; i++){
    delete[] clusters[i]; clusters[i] = NULL;
  }
  delete[] clusters; clusters = NULL;
  delete[] cluster_assignment; cluster_assignment = NULL;
  for(i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;
  delete[] log_prior; log_prior = NULL;
  indexing.resize(0);
  // delete[] log_det_Omegay; log_det_Omegay = NULL;
  // delete[] quad_forms; quad_forms = NULL;
  return;
}

void Partition::Initialize_Partition(int n, double eta){
  nObs = n;
  K = 1;
  cluster_config = new int[K];
  cluster_config[0] = nObs; // initial partition has 1 cluster of size nObs
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  indexing.resize(nObs);
  // log_det_Omegay = new double[K];
  // quad_forms = new double[K];
  log_prior = new double[K];

  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
  }
  // now since K = 1, we only have one cluster:
  for(int i = 0; i < nObs; i++){
    clusters[0][i] = i;
  }
  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = 0; // initial partition has 1 cluster
    indexing[i] = i;
  }
  get_pairwise();
  for(int k = 0; k < K; k++){
    log_pi_ep(k, eta);
  }

  // get_likelihood(A_or_B, 0);
  return;
}

void Partition::Initialize_Partition(int n, std::vector<int> ids, double eta){
  // ids contains the mainID (right?) so here we are putting mainIDs in clusters
  nObs = n;
  K = 1;
  cluster_config = new int[K];
  cluster_config[0] = nObs; // initial partition has 1 cluster of size nObs
  clusters = new int*[K];
  cluster_assignment = new int[nObs];

  // log_det_Omegay = new double[K];
  // quad_forms = new double[K];
  log_prior = new double[K];

  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
  }
  // now since K = 1, we only have one cluster:
  for(int i = 0; i < nObs; i++){
    clusters[0][i] = i;
  }
  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = 0; // initial partition has 1 cluster
  }
  indexing = std::vector<int>(ids);
  get_pairwise();
  for(int k = 0; k < K; k++){
    log_pi_ep(k, eta);
  }
  // get_likelihood(A_or_B, 0);
  return;
}

void Partition::Initialize_Partition(int n, std::vector<int> ids, double eta, int gamma_ptr[], int* K_init_ptr){
  nObs = n;
  K = *K_init_ptr;
  cluster_config = new int[K];
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  for(int k = 0; k < K; k++){
    cluster_config[k] = 0;
    for(int i = 0; i < n; i++){
      if(gamma_ptr[ids[i]] == (k+1)){
        cluster_config[k]++;
      }
    }
  }
  int count;
  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
    count = 0;
    for(int i = 0; i < n; i++){
      if(gamma_ptr[ids[i]] == (k+1)){
        clusters[k][count] = i;
        count++;
      }
    }
  }
  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = gamma_ptr[ids[i]]-1;
  }
  log_prior = new double[K];
  get_pairwise();
  for(int k = 0; k < K; k++){
    log_pi_ep(k, eta);
  }
  indexing = std::vector<int>(ids);
  return;
}

std::vector<int> Partition::Initialize_Partition(int n, std::vector<int> ids, double eta, int gamma_ptr[]){
  // returns the unique values for each cluster, in a vector long as the number of clusters of ids
  // works even if length of gamma_ptr is longer than length of ids, n is length of ids
  nObs = n;
  // Find the number of clusters, n_unique, and the unique values, gamma_unique:
  std::vector<int> gamma_unique;
  int n_unique = 0;
  for(int i = 0; i < n; i++){
    if(n_unique == 0){
      gamma_unique.push_back(gamma_ptr[ids[i]]);
      n_unique++;
    } else {
      bool newcl = true;
      for(int k = 0; k < n_unique; k++){
        if(gamma_ptr[ids[i]] == gamma_unique[k]){
          newcl = false;
        }
      }
      if(newcl){
        gamma_unique.push_back(gamma_ptr[ids[i]]);
        n_unique++;
      } 
    }
  }
  K = n_unique;

  cluster_config = new int[K];
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  // find cluster size, cluster_config
  for(int k = 0; k < K; k++){
    cluster_config[k] = 0;
    for(int i = 0; i < n; i++){
      if(gamma_ptr[ids[i]] == gamma_unique[k]){
        cluster_config[k]++;
      }
    }
  }
  // find clusters and cluster_assignment 
  int count;
  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
    count = 0;
    for(int i = 0; i < n; i++){
      if(gamma_ptr[ids[i]] == gamma_unique[k]){
        clusters[k][count] = i;
        cluster_assignment[i] = k;
        count++;
      }
    }
  }
  log_prior = new double[K];
  get_pairwise();
  for(int k = 0; k < K; k++){
    log_pi_ep(k, eta);
  }
  indexing = std::vector<int>(ids);
  return(gamma_unique);
}

void Partition::Initialize_Partition(int n, double eta, std::vector<int> gamma_ptr, std::vector<int> gamma_unique){
  nObs = n;
  K = (int)gamma_unique.size();
  cluster_config = new int[K];
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  for(int k = 0; k < K; k++){
    cluster_config[k] = 0;
    for(int i = 0; i < n; i++){
      if(gamma_ptr[i] == k){
        cluster_config[k]++;
      }
    }
  }
  int count;
  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
    count = 0;
    for(int i = 0; i < n; i++){
      if(gamma_ptr[i] == k){
        clusters[k][count] = i;
        cluster_assignment[i] = k;
        count++;
      }
    }
  }
  log_prior = new double[K];
  get_pairwise();
  for(int k = 0; k < K; k++){
    log_pi_ep(k, eta);
  }
  indexing.resize(nObs);
  for(int i = 0; i < nObs; i++){
    indexing[i] = i;
  }
  return;
}

void Partition::Initialize_CostDish(std::vector<int> nObsRest, std::vector<int> nTableRest, std::vector<LPPartition> costumers_tables, LPPartition tables_dishes){
  nObs = 0;
  for(int r = 0; r < nObsRest.size(); r++){
    nObs += nObsRest[r];
  }
  int x, y; // x, y are used for main ID
  int rest, table, cost; // rest, table, cost are used for restaurant ID
  K = tables_dishes->K;
  cluster_config = new int[K];
  cluster_assignment = new int[nObs];
  vector<vector<int> > clusters_v(K); // using a vector to make things easier
  for(int k = 0; k < K; k++){
    // let's form the kth cluster of costumer_dishes
    // in *cluster k of costumers_dishes* I will get the costumers at the tables in *cluster k of tables_dishes*
    cluster_config[k] = 0;
    for(int i = 0; i < tables_dishes->cluster_config[k]; i++){
      y = tables_dishes->clusters[k][i];
      // from the table's main id we get the restaurant id:
      Table_MainId_RestId(y, nTableRest, rest, table); // get rest and table
      // cout << "Element "<<i<<" in cluster "<<k <<" of tables_dishes. ";
      // cout << "Corresponds to " <<y<< " in the mainID. ";
      // cout << "Corresponds to rest " <<rest<< " and table " <<table << endl;      
      cluster_config[k] += costumers_tables[rest]->cluster_config[table]; 
      for(int ii = 0; ii < costumers_tables[rest]->cluster_config[table]; ii++){
        cost = costumers_tables[rest]->clusters[table][ii];
        // from the costumer's restaurant id we get the main id
        x = costumers_tables[rest]->indexing[cost];
        // x = Costumer_RestId_MainId(rest, cost, nObsRest); // need to transform into the MainID
        clusters_v[k].push_back( x ); // which people at that table  
        cluster_assignment[ x ] = k;
      }
    }
  }
  clusters = new int*[K];
  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
    for(int ii = 0; ii < cluster_config[k]; ii++){
      clusters[k][ii] = clusters_v[k][ii];
    }
  }
  log_prior = new double[K];
  get_pairwise();
  // The prior of this partition is not described by the ewens formula, so I don't know the prior!
  for(int k = 0; k < K; k++){
    // log_pi_ep(k);
    log_prior[k] = 0;
  }
  indexing = std::vector<int>(nObs);
  for(int i = 0; i < nObs; i++){
    indexing[i] = i;
  }
  return;
}

// void Partition::Initialize_Partition(int n, Rcpp::List gamma_init){
//   nObs = n;
//   K = gamma_init.size();
//   cluster_config = new int[K];
//   cluster_assignment = new int[nObs];
//   clusters = new int*[K];
//   // log_det_Omegay = new double[K];
//   // quad_forms = new double[K];
//   log_prior = new double[K];
  
//   Rcpp::NumericVector tmp_vec;
//   for(int k = 0; k < K; k++){
//     tmp_vec = gamma_init[k];
//     cluster_config[k] = tmp_vec.size();
//     clusters[k] = new int[cluster_config[k]];
//     for(int i = 0; i < cluster_config[k]; i++){
//       clusters[k][i] = tmp_vec[i] - 1; // R is 1-indexed and C++ is 0-indexed
//       cluster_assignment[clusters[k][i]] = k;
//     }
//     log_pi_ep(k);
//   }
//   get_pairwise();
//   // for(int k = 0; k < K; k++){
//   //   get_likelihood(A_or_B, k);
//   // }
//   return;
// }

void Partition::Initialize_Partition_nclusters(int n, double eta){
  nObs = n;
  K = n;
  cluster_config = new int[K];
  clusters = new int*[K];
  indexing.resize(nObs);
  for(int k = 0; k < K; k++){
    cluster_config[k] = 1;
    clusters[k] = new int[1];
    clusters[k][0] = k;
  }
  cluster_assignment = new int[nObs];
  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = i; // initial partition has 1 cluster
    indexing[i] = i;
  }
  // log_det_Omegay = new double[K];
  // quad_forms = new double[K];
  log_prior = new double[K];
  get_pairwise();
  for(int k = 0; k < K; k++){
    // get_likelihood(A_or_B, k);
    log_pi_ep(k, eta);
  }
  return;
}

void Partition::Initialize_Partition_nclusters(int n, std::vector<int> ids, double eta){
  nObs = n;
  K = n;
  cluster_config = new int[K];
  clusters = new int*[K];
  for(int k = 0; k < K; k++){
    cluster_config[k] = 1;
    clusters[k] = new int[1];
    clusters[k][0] = k;
  }
  cluster_assignment = new int[nObs];
  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = i; // initial partition has 1 cluster
  }
  // log_det_Omegay = new double[K];
  // quad_forms = new double[K];
  log_prior = new double[K];
  get_pairwise();
  indexing = std::vector<int>(ids);
  for(int k = 0; k < K; k++){
    // get_likelihood(A_or_B, k);
    log_pi_ep(k, eta);
  }
  return;
}


void Partition::get_pairwise(){
  pairwise_assignment = new int*[nObs];
  for(int i = 0; i < nObs; i++){
    pairwise_assignment[i] = new int[nObs];
  }
  for(int i = 0; i < nObs; i++){
    for(int j = i; j < nObs; j++){
      if(cluster_assignment[i] == cluster_assignment[j]){
        pairwise_assignment[i][j] = 1;
        pairwise_assignment[j][i] = 1;
      } else{
        pairwise_assignment[i][j] = 0;
        pairwise_assignment[j][i] = 0;
      }
    }
  }
  return;
}

void Partition::Print_Partition(){
  std::cout << "Number of clusters K: " << K << std::endl;
  std::cout << "Size of clusters:";
  for(int k = 0; k < K; k++){
    std::cout << cluster_config[k] << " ";
  }
  std::cout << std::endl;
  std::cout << "Clusters:" << std::endl;
  for(int k = 0; k < K; k++){
    std::cout << "Cluster " << k  << " : ";
    for(int j = 0; j < cluster_config[k]; j++){
      std::cout << clusters[k][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Log-prior:" ;
  for(int k = 0; k < K; k++){
    std::cout <<  log_prior[k] << " ";
  }
  std::cout << std::endl;
  std::cout << "Indexing:" ;
  for(int i = 0; i < indexing.size(); i++){
    std::cout <<  indexing[i] << " ";
  }
  std::cout << std::endl;
  return;
}
void Partition::Print_Partition_Short(){
  std::cout << "Number of clusters K: " << K << std::endl;
  std::cout << "Size of clusters: ";
  for(int k = 0; k < K; k++){
  std::cout << cluster_config[k] << " ";
  }
  std::cout << std::endl;
  return;
}
void Partition::Print_Partition_ToFile(string file){
  ofstream myfile;
  myfile.open(file.c_str(), std::ofstream::app);
  // myfile.open(file, std::ios_base::app);
  for(int i = 0; i < nObs-1; i++){
    for(int j = i+1; j < nObs; j++){
      myfile << pairwise_assignment[i][j] << ",";
    }
  }
  myfile << endl;
  myfile.close();
  return;
}
void Partition::Read_Partition_FromFile(string file, int n){
  ifstream filein ( file.c_str() );
  string value;
  int value_int;

  // mat mat_temp(n,n,fill::randu);
  int aa = 0;
  pairwise_assignment = new int*[nObs];
  for(int i = 0; i < nObs; i++){
    pairwise_assignment[i] = new int[nObs];
  }
  for(int i = 0; i < nObs-1; i++){
    pairwise_assignment[i][i] = 1;
    // mat_temp(i,i) = 1;
    for(int j = i+1; j < nObs; j++){
      getline ( filein, value, ',' );
      value_int = stoi(value);
      pairwise_assignment[j][i] = value_int;
      pairwise_assignment[i][j] = value_int;
      // mat_temp(i,j) = value_int;
      // mat_temp(j,i)f = value_int;
      aa++;
    }
  }  
  filein.close();
  return;
}
void Partition::Print_Y(bool LR){
  arma::mat Yloc;
  if(LR){
    Yloc = YLR;
  } else {
    Yloc = Y;
  }
  int cost_localID, cost_mainID;
  for(int k =0; k < K; k++){
    for(int t = 0; t < cluster_config[k]; t++){
      cost_localID = clusters[k][t];
      cost_mainID = indexing[cost_localID];
      cout << Yloc( cost_mainID ) << " ";
    }
    cout << endl;
  }
}

void Partition::log_pi_ep(int cluster_id, double eta){
  log_prior[cluster_id] = log(eta) + lgamma(cluster_config[cluster_id]);
}

double Partition::get_logprior(double eta){
  // cout << "get_logprior: ";
  double logpr = -lgamma(eta + nObs) + lgamma(eta);
  for(int k = 0; k < K; k++){
    // cout << cluster_config[k] << " ";
    logpr += log_prior[k];
  }
  // cout << endl;
  return logpr;
}

// double Partition::get_logprior2(double eta){
//   double logpr = -lgamma(eta + nObs) + lgamma(eta);
//   for(int k = 0; k < K; k++){
//     logpr += log(eta) + lgamma(cluster_config[k]);
//   }
//   return logpr;
// }

void Partition::AddTable_TableDish(int r, int newtable_id, int cluster_newtable, std::vector<int> nTableRest, double eta_TD){
  // since we added a table, we need to shift the indices of the other tables and change the whole partition structure
  int orig_nObs = nObs;
  nObs = orig_nObs+1;
  // cout << r << " " << newtable_id << " " << cluster_newtable << endl;
  if(cluster_newtable < K){ // it's a current dish
    // no need to change K or whole cluster_config, or whole log_prior
    // we need to change whole clusters, whole pairwise because the indexing has changed
    int** orig_clusters = new int*[K];
    int* orig_cluster_assignment = new int[orig_nObs];
    for(int k = 0; k < K; k++){
      orig_clusters[k] = new int[cluster_config[k]];
      for(int i = 0; i < cluster_config[k]; i++){
        orig_clusters[k][i] = clusters[k][i];
      }
    }
    for(int i = 0; i < orig_nObs; i++){
      orig_cluster_assignment[i] = cluster_assignment[i];
    }

    for(int k = 0; k < K; k++){
      delete[] clusters[k]; clusters[k] = NULL;
    }
    delete[] clusters; clusters = NULL;
    delete[] cluster_assignment; cluster_assignment = NULL;
    for(int i = 0; i < orig_nObs; i++){
      delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
    }
    delete[] pairwise_assignment; pairwise_assignment = NULL;


    cluster_config[cluster_newtable] += 1;
    clusters = new int*[K];
    cluster_assignment = new int[nObs];
    int y = Table_RestId_MainId(r, newtable_id, nTableRest);
    // cout << "y: " << y << endl;
    int yy;
    for(int k = 0; k < K; k++){
      clusters[k] = new int[cluster_config[k]];
      if(k == cluster_newtable){
        bool y_not_assigned = true;
        int count = 0;
        for(int i = 0; i < cluster_config[k]-1; i++){ // only up to ..-1 because we access the orig 
          yy = orig_clusters[k][i];
          // cout << "k, i, yy " << k << " " << i << " " << yy;
          if(yy < y){
            clusters[k][count] = yy;
            count++;
            // cout << " " << clusters[k][i] << endl;
          } else if(yy >= y){
            if(y_not_assigned){
              clusters[k][count] = y;
              count++;
              // cout << " y " << clusters[k][i]; 
              y_not_assigned = false;
              clusters[k][count] = yy + 1;
              count++;
              // cout << " " << clusters[k][i+1] << endl;
            } else {
              clusters[k][count] = yy + 1;
              count++;
              // cout << " " << clusters[k][i] << endl;
            }
          }
        }
        if(y_not_assigned){
          // clusters[k][cluster_config[k]-1] = y;
          clusters[k][count] = y;
          // cout << "y " << k << " " << cluster_config[k] -1 << " " << clusters[k][cluster_config[k]-1] << endl;
        }
      } else {
        for(int i = 0; i < cluster_config[k]; i++){ // only up to ..-1 because we access the orig 
          yy = orig_clusters[k][i];
          // cout << "k, i, yy " << k << " " << i << " " << yy;
          if(yy < y){
            clusters[k][i] = yy;
            // cout << " " << clusters[k][i] << endl;
          } else if(yy >= y){
            clusters[k][i] = yy + 1;
            // cout << " " << clusters[k][i] << endl;
          }
        }
      }
    }
    cluster_assignment[y] = cluster_newtable;
    for(int yy = 0; yy < orig_nObs; yy++){
      if(yy < y){
        cluster_assignment[yy] = orig_cluster_assignment[yy];
      } else {
        cluster_assignment[yy+1] = orig_cluster_assignment[yy];
      }
    }
    get_pairwise();
    log_pi_ep(cluster_newtable, eta_TD);
    for(int kk = 0; kk < K; kk++){
      delete[] orig_clusters[kk];
    }
    delete[] orig_clusters;
    delete[] orig_cluster_assignment;
  } else {
    int orig_K = K;
    int* orig_cluster_config = new int[orig_K];
    int** orig_clusters = new int*[orig_K];
    int* orig_cluster_assignment = new int[orig_nObs];
    double* orig_log_prior = new double[orig_K];
    for(int k = 0; k < orig_K; k++){
      orig_cluster_config[k] = cluster_config[k];
      orig_clusters[k] = new int[cluster_config[k]];
      for(int i = 0; i < cluster_config[k]; i++){
        orig_clusters[k][i] = clusters[k][i];
      }
      orig_log_prior[k] = log_prior[k];
    }
    for(int i = 0; i < orig_nObs; i++){
      orig_cluster_assignment[i] = cluster_assignment[i];
    }
    
    delete[] cluster_config; cluster_config = NULL;
    for(int k = 0; k < orig_K; k++){
      delete[] clusters[k]; clusters[k] = NULL;
    }
    delete[] clusters; clusters = NULL;
    delete[] cluster_assignment; cluster_assignment = NULL;
    delete[] log_prior; log_prior = NULL;
    for(int i = 0; i < orig_nObs; i++){
      delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
    }
    delete[] pairwise_assignment; pairwise_assignment = NULL;


    K = orig_K + 1;
    cluster_config = new int[K];
    clusters = new int*[K];
    cluster_assignment = new int[nObs];
    log_prior = new double[K];

    for(int k = 0; k < K; k++){
      if(k == K-1){
        cluster_config[k] = 1;
      } else {
        cluster_config[k] = orig_cluster_config[k];
      }
    }
    // How do the ID change? the new cluster has ID = y.
    // If your ID < y, then you don't change.
    // If your ID was >= y, then you become y+1
    int y = Table_RestId_MainId(r, newtable_id, nTableRest);
    int yy;
    for(int k = 0; k < K; k++){
      clusters[k] = new int[cluster_config[k]];
      if(k == K-1){
        clusters[k][0] = y;
      } else{
        for(int i = 0; i < cluster_config[k]; i++){
          yy = orig_clusters[k][i];
          if(yy < y){
            clusters[k][i] = yy;
          } else {
            clusters[k][i] = yy+1;
          }
        }
      }
    }
    cluster_assignment[y] = K-1;
    for(int yy = 0; yy < orig_nObs; yy++){
      if(yy < y){
        cluster_assignment[yy] = orig_cluster_assignment[yy];
      } else {
        cluster_assignment[yy+1] = orig_cluster_assignment[yy];
      }
    }
    get_pairwise();
    for(int k = 0; k < K; k++){
      if(k == K-1){ // need to re-compute
        log_pi_ep(k, eta_TD);
      } else {
        log_prior[k] = orig_log_prior[k];
      }
    }
    delete[] orig_cluster_config;
    for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length orig_K (K and not K+1.)
      delete[] orig_clusters[kk];
    }
    delete[] orig_clusters;
    delete[] orig_cluster_assignment;
    delete[] orig_log_prior;
  }
  indexing.resize(nObs);
  indexing[nObs-1] = nObs-1;
  return;
}

void Partition::MoveTables_TableDish(std::vector<int> tables, int starting_index){
  /*
  The tables here are divided as "before rmin", "rmin", "middle", "rmax", "after rmax"
  and we need to reorder them as "before rmin", "rmin", "rmax", "middle", "after rmax"
  thus we need to change the names of 
  - "middle" by shifting forward of rmax->K (= tsize) [new name is t + tsize]
  - "rmax" by shifting back, the new name starting_index + (t-tables[0])
    Note that (t-tables[0]) ranges from 0..K-1 so we shift that by starting_index.
  We need to change only clusters, cluster_assignment and indexing, (and pairwise)
  while nObs, K, cluster_config and log_prior do not change
  
  table 0 of tables becomes table starting_index, 
  and table starting_index becomes starting_index+tables.size()
  all the tables that before where after starting_index now need to get shifted
  */

  
  int t;
  int tsize = tables.size();
  for(int i = 0; i < K; i++){ // i is index for dish
    for(int j = 0; j < cluster_config[i]; j++){ // j is index for table
      t = clusters[i][j]; // t is the jth table in cluster i (mainID)
      if(t >= starting_index && t < tables[0]){ 
        // if it's in the middle between rmin and rmax
        clusters[i][j] = t + tsize;
      } else if(t >= tables[0] && t <= tables[tsize-1]){ 
        // if it's the old rmax
        clusters[i][j] = starting_index + (t - tables[0]);
      }
      // otherwise it does not change
    }
  }

  int* orig_cluster_assignment = new int[nObs];
  for(int i = 0; i < nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }
  std::vector<int> orig_indexing = std::vector<int>(indexing);

  for(int i = 0; i < nObs; i++){
    if(i >= starting_index && i < tables[0]){
      // if it's in the middle between rmin and rmax
      cluster_assignment[i + tsize] = orig_cluster_assignment[i];  
      indexing[i + tsize] = orig_indexing[i];
    } else if (i >= tables[0] && i <= tables[tsize-1]){
      // if it's the old rmax
      cluster_assignment[starting_index + (i - tables[0])] = orig_cluster_assignment[i];
      indexing[starting_index + (i - tables[0])] = orig_indexing[i];
    } 
  }
  for(int i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; 
  }
  delete[] pairwise_assignment; 
  get_pairwise();
  delete[] orig_cluster_assignment;
  return;
}

// void Partition::SplitTables_TableDish(int r1, int r2, std::vector<int> tables1, std::vector<int> tables2, std::vector<int> nTableRest, double eta_TD, bool local_print){
//   /*
//   We update tables_dishes after the split of a LR cluster, which splits some tables overlapping with the two splits
//   We need to change nObs, and the names of some tables.
//   We remove the tables that were initially in r1 (nTableRest[r1]) and add instead the new tables1. 
//   Then we have the ones which were after r1, which change name, and then we add r2 tables (tables2) at the end
//   Note: tables1/2 contain the localID in the old r1.
//   The number of dishes stays the same.
//   */
//   if(r2 != nTableRest.size()){
//     cout << "SplitTables_TableDish: Problem!" << endl;
//   }
//   int orig_nObs = nObs;
//   int orig_K = K;
//   int* orig_cluster_assignment = new int[orig_nObs];
//   for(int i = 0; i < orig_nObs; i++){
//     orig_cluster_assignment[i] = cluster_assignment[i];
//   }
  
//   delete[] cluster_config; cluster_config = NULL;
//   for(int k = 0; k < orig_K; k++){
//     delete[] clusters[k]; clusters[k] = NULL;
//   }
//   delete[] clusters; clusters = NULL;
//   delete[] cluster_assignment; cluster_assignment = NULL;
//   delete[] log_prior; log_prior = NULL;
//   for(int i = 0; i < orig_nObs; i++){
//     delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
//   }
//   delete[] pairwise_assignment; pairwise_assignment = NULL;

//   nObs = orig_nObs - nTableRest[r1] + tables1.size() + tables2.size();
//   // the shift is nTableRest[r1] - tables1.size()
//   // int shift = nTableRest[r1] - tables1.size(); // not used WHY?
//   // the shift should happen for the ID that are >= than cumsum(nTableRest)[r1]
//   // the change in names happens for the ID >= starting_index
//   int starting_index = nTableRest[r1];
//   for(int i = 0; i < r1; i++){
//     starting_index += nTableRest[i];
//   }

//   // we go through orig_cluster_assignment in the right order and reconstruct all the pieces
//   // cl here represents a dish - we consider the original dish and we keep it when we split the tables
//   std::vector<int> cluster_assignment_v, cluster_assignment_v2;
//   std::vector<std::vector<int> > clusters_v, clusters_v2;
//   std::vector<int> cluster_touched;
//   int cl, ind_oldcl, localID, mainID, iID;
//   cluster_assignment_v.reserve(nObs);
//   cluster_assignment_v2.reserve(nObs);
//   clusters_v2.resize(orig_K);
  
//   iID = 0; // new id for the tables (in the first loop coincides with i)
  
//   // we first look at the things before r1, which don't change
//   for(int i = 0; i < starting_index-nTableRest[r1]; i++){
//     cl = orig_cluster_assignment[i];
//     bool newcl = true;
//     for(int j = 0; j < cluster_touched.size(); j++){
//       if(cluster_touched[j] == cl){
//         ind_oldcl = j;
//         newcl = false;
//         break;
//       }
//     }
//     if(newcl){
//       cluster_touched.push_back(cl);
//       clusters_v.push_back(std::vector<int>(1,iID));
//       cluster_assignment_v.push_back(cluster_touched.size() - 1);
//     } else {
//       clusters_v[ind_oldcl].push_back(iID);
//       cluster_assignment_v.push_back(ind_oldcl);
//     }
//     if(local_print){
//       printf("Table %d (with iID %d) in cluster %d which is (%d) new.\n",i,iID,cl,newcl);
//     }
//     iID++;
//   }
//   // then we skip r1 and look at tables1
//   for(int i = 0; i < tables1.size(); i++){
//     localID = tables1[i];
//     mainID = starting_index-nTableRest[r1] + localID; // this is the old mainID
//     cl = orig_cluster_assignment[mainID];
//     bool newcl = true;
//     for(int j = 0; j < cluster_touched.size(); j++){
//       if(cluster_touched[j] == cl){
//         ind_oldcl = j;
//         newcl = false;
//         break;
//       }
//     }
//     if(newcl){
//       cluster_touched.push_back(cl);
//       clusters_v.push_back(std::vector<int>(1,iID));
//       cluster_assignment_v.push_back(cluster_touched.size() - 1);
//     } else {
//       clusters_v[ind_oldcl].push_back(iID);
//       cluster_assignment_v.push_back(ind_oldcl);
//     }
//     if(local_print){
//       printf("Table %d (with mainID %d and iID %d) in cluster %d which is (%d) new.\n",localID,mainID,iID,cl,newcl);
//     }
//     iID++;
//   }
//   // then we look at the things after r1, and before the original end, which don't change
//   for(int i = starting_index; i < orig_nObs; i++){
//     cl = orig_cluster_assignment[i];
//     bool newcl = true;
//     for(int j = 0; j < cluster_touched.size(); j++){
//       if(cluster_touched[j] == cl){
//         ind_oldcl = j;
//         newcl = false;
//         break;
//       }
//     }
//     if(newcl){
//       cluster_touched.push_back(cl);
//       clusters_v.push_back(std::vector<int>(1,iID));
//       cluster_assignment_v.push_back(cluster_touched.size() - 1);
//     } else {
//       clusters_v[ind_oldcl].push_back(iID);
//       cluster_assignment_v.push_back(ind_oldcl);
//     }
//     if(local_print){
//       printf("Table %d (with iID %d) in cluster %d which is (%d) new.\n",i,iID,cl,newcl);
//     }
//     iID++;
//   }
//   // finally we do tables2 with goes at the end
//   for(int i = 0; i < tables2.size(); i++){
//     localID = tables2[i];
//     mainID = starting_index-nTableRest[r1] + localID;
//     cl = orig_cluster_assignment[mainID];
//     bool newcl = true;
//     for(int j = 0; j < cluster_touched.size(); j++){
//       if(cluster_touched[j] == cl){
//         ind_oldcl = j;
//         newcl = false;
//         break;
//       }
//     }
//     if(newcl){
//       cluster_touched.push_back(cl);
//       clusters_v.push_back(std::vector<int>(1,iID));
//       cluster_assignment_v.push_back(cluster_touched.size() - 1);
//     } else {
//       clusters_v[ind_oldcl].push_back(iID);
//       cluster_assignment_v.push_back(ind_oldcl);
//     }
//     if(local_print){
//       printf("Table %d (with mainID %d and iID %d) in cluster %d which is (%d) new.\n",localID,mainID,iID,cl,newcl);
//     }
//     iID++;
//   }


//   if(local_print){
//     iID = 0;
    
//     // we first look at the things before r1, which don't change
//     for(int i = 0; i < starting_index-nTableRest[r1]; i++){
//       cl = orig_cluster_assignment[i];
//       clusters_v2[cl].push_back(iID);
//       cluster_assignment_v2.push_back(cl);
//       iID++;
//     }
//     // then we skip r1 and look at tables1
//     for(int i = 0; i < tables1.size(); i++){
//       localID = tables1[i];
//       mainID = starting_index-nTableRest[r1] + localID; // this is the old mainID
//       cl = orig_cluster_assignment[mainID];
//       clusters_v2[cl].push_back(iID);
//       cluster_assignment_v2.push_back(cl);
//       iID++;
//     }
//     // then we look at the things after r1, and before the original end, which don't change
//     for(int i = starting_index; i < orig_nObs; i++){
//       cl = orig_cluster_assignment[i];
//       clusters_v2[cl].push_back(iID);
//       cluster_assignment_v2.push_back(cl);
//       iID++;
//     }
//     // finally we do tables2 with goes at the end
//     for(int i = 0; i < tables2.size(); i++){
//       localID = tables2[i];
//       mainID = starting_index-nTableRest[r1] + localID;
//       cl = orig_cluster_assignment[mainID];
//       clusters_v2[cl].push_back(iID);
//       cluster_assignment_v2.push_back(cl);
//       iID++;
//     }
//     // cout << "Cluster assignment_v2: " << endl;
//     // for(int i = 0; i < nObs; i++){
//     //   cout << cluster_assignment_v2[i] << " ";
//     // }
//     // cout << "Cluster_v2: " << endl;
//     // for(int k = 0; k < clusters_v2.size(); k++){
//     //   cout << "Cluster " << k << ": ";
//     //   for(int j = 0; j < clusters_v2[k].size(); j++){
//     //     cout << clusters_v2[k][j] << " ";
//     //   } 
//     //   cout << endl;
//     // }
//   }
  

//   K = clusters_v.size();
//   cluster_config = new int[K];
//   clusters = new int*[K];
//   cluster_assignment = new int[nObs];
//   log_prior = new double[K];
//   for(int i = 0; i < K; i++){
//     cluster_config[i] = clusters_v[i].size();
//     clusters[i] = new int[cluster_config[i] ];
//     for(int j = 0; j < cluster_config[i]; j++){
//       clusters[i][j] = clusters_v[i][j];
//     }
//   }
//   if(local_print){
//     cout << "Cluster assignment: " << endl;
//   }
//   for(int i = 0; i < nObs; i++){
//     cluster_assignment[i] = cluster_assignment_v[i];
//     if(local_print){
//       cout << cluster_assignment_v[i] << " ";
//     }
//   }
//   if(local_print) cout << endl;
//   for(int k = 0; k < K; k++){
//     log_pi_ep(k, eta_TD);
//   }
//   get_pairwise();
//   indexing = std::vector<int>(nObs);
//   for(int i = 0; i < nObs; i++){
//     indexing[i] = i;
//   }

//   if(local_print){
//     Print_Partition();
//   }
//   // delete[] orig_cluster_config; 
//   // for(int k = 0; k < orig_K; k++){
//   //   delete[] orig_clusters[k]; 
//   // }
//   // delete[] orig_clusters; 
//   delete[] orig_cluster_assignment; 
//   // delete[] orig_log_prior; 
//   return;
// }

// this version removes the change of order in the dishes (keep the same order and names as before)
void Partition::SplitTables_TableDish(int r1, int r2, std::vector<int> tables1, std::vector<int> tables2, std::vector<int> nTableRest, double eta_TD, bool local_print){
  /*
  We update tables_dishes after the split of a LR cluster, which splits some tables overlapping with the two splits
  We need to change nObs, and the names of some tables.
  We remove the tables that were initially in r1 (nTableRest[r1]) and add instead the new tables1. 
  Then we have the ones which were after r1, which change name, and then we add r2 tables (tables2) at the end
  Note: tables1/2 contain the localID in the old r1.
  The number of dishes stays the same.
  DOES NOT CHANGE THE NAME/ORDER OF THE DISHES
  */
  if(r2 != nTableRest.size()){
    cout << "SplitTables_TableDish: Problem!" << endl;
  }
  int orig_nObs = nObs;
  int orig_K = K;
  int* orig_cluster_assignment = new int[orig_nObs];
  for(int i = 0; i < orig_nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }
  
  delete[] cluster_config; cluster_config = NULL;
  for(int k = 0; k < orig_K; k++){
    delete[] clusters[k]; clusters[k] = NULL;
  }
  delete[] clusters; clusters = NULL;
  delete[] cluster_assignment; cluster_assignment = NULL;
  delete[] log_prior; log_prior = NULL;
  for(int i = 0; i < orig_nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;

  nObs = orig_nObs - nTableRest[r1] + tables1.size() + tables2.size();
  // the shift is nTableRest[r1] - tables1.size()
  // int shift = nTableRest[r1] - tables1.size(); // not used WHY?
  // the shift should happen for the ID that are >= than cumsum(nTableRest)[r1]
  int starting_index = nTableRest[r1];
  for(int i = 0; i < r1; i++){
    starting_index += nTableRest[i];
  }

  // we go through orig_cluster_assignment in the right order and reconstruct all the pieces
  // cl here represents a dish - we conside the original dish and we keep it when we split the tables
  std::vector<int> cluster_assignment_v;
  std::vector<std::vector<int> > clusters_v;
  std::vector<int> cluster_touched;
  int cl, localID, mainID, iID;
  cluster_assignment_v.reserve(nObs);
  clusters_v.resize(orig_K);
  
  iID = 0; // new id for the tables (in the first loop coincides with i)
  
  // we first look at the things before r1, which don't change
  for(int i = 0; i < starting_index-nTableRest[r1]; i++){
    cl = orig_cluster_assignment[i];
    clusters_v[cl].push_back(iID);
    cluster_assignment_v.push_back(cl);
    iID++;
  }
  // then we skip r1 and look at tables1
  for(int i = 0; i < tables1.size(); i++){
    localID = tables1[i];
    mainID = starting_index-nTableRest[r1] + localID; // this is the old mainID
    cl = orig_cluster_assignment[mainID];
    clusters_v[cl].push_back(iID);
    cluster_assignment_v.push_back(cl);
    iID++;
  }
  // then we look at the things after r1, and before the original end, which don't change
  for(int i = starting_index; i < orig_nObs; i++){
    cl = orig_cluster_assignment[i];
    clusters_v[cl].push_back(iID);
    cluster_assignment_v.push_back(cl);
    iID++;
  }
  // finally we do tables2 with goes at the end
  for(int i = 0; i < tables2.size(); i++){
    localID = tables2[i];
    mainID = starting_index-nTableRest[r1] + localID;
    cl = orig_cluster_assignment[mainID];
    clusters_v[cl].push_back(iID);
    cluster_assignment_v.push_back(cl);
    iID++;
  }
  // cout << "Cluster assignment_v: " << endl;
  // for(int i = 0; i < nObs; i++){
  //   cout << cluster_assignment_v[i] << " ";
  // }
  // cout << "Cluster_v: " << endl;
  // for(int k = 0; k < clusters_v.size(); k++){
  //   cout << "Cluster " << k << ": ";
  //   for(int j = 0; j < clusters_v[k].size(); j++){
  //     cout << clusters_v[k][j] << " ";
  //   } 
  //   cout << endl;
  // }

  K = clusters_v.size();
  cluster_config = new int[K];
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  log_prior = new double[K];
  for(int i = 0; i < K; i++){
    cluster_config[i] = clusters_v[i].size();
    clusters[i] = new int[cluster_config[i] ];
    for(int j = 0; j < cluster_config[i]; j++){
      clusters[i][j] = clusters_v[i][j];
    }
  }
  // if(local_print){
  //   cout << "Cluster assignment: " << endl;
  // }
  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = cluster_assignment_v[i];
    // if(local_print){
    //   cout << cluster_assignment_v[i] << " ";
    // }
  }
  // if(local_print) cout << endl;
  for(int k = 0; k < K; k++){
    log_pi_ep(k, eta_TD);
  }
  get_pairwise();
  indexing = std::vector<int>(nObs);
  for(int i = 0; i < nObs; i++){
    indexing[i] = i;
  }

  // if(local_print){
  //   Print_Partition();
  // }
  // delete[] orig_cluster_config; 
  // for(int k = 0; k < orig_K; k++){
  //   delete[] orig_clusters[k]; 
  // }
  // delete[] orig_clusters; 
  delete[] orig_cluster_assignment; 
  // delete[] orig_log_prior; 
  return;
}


// void Partition::AddTables_TableDish(int r, std::vector<int> newtables_id, std::vector<int> cluster_newtables, std::vector<int> nTableRest){
//   int num_splits = newtables_id.size() + 1;
//   int orig_nObs = nObs;
//   int orig_K = K;
//   nObs = orig_nObs + num_splits - 1;

//   int reindexing_threshold = 0;
//   for(int rr = 0; rr < r; rr++){
//     reindexing_threshold += nTableRest[rr];
//   }
//   reindexing_threshold += nTableRest[r] - num_splits + 1;
//   int reindexing_shift = num_splits - 1;

//   // figure out which new tables go into new clusters and which ones into old clusters
//   std::vector<int> new_clusters;
//   std::vector<int> old_clusters(orig_K,0);
//   std::vector<std::vector<int> > in_old_clusters;
//   std::vector<std::vector<int> > in_new_clusters;
//   bool new_cluster_flag, new_old_cluster_flag;
//   for(int h = 0; h < num_splits - 1; h++){
//     if(cluster_newtables[h] >= orig_K){
//       new_cluster_flag = true;
//       for(int j = 0; j < new_clusters.size(); j++){
//         if(cluster_newtables[h] == new_clusters[j]){
//           new_cluster_flag = false;
//           in_new_clusters[j].push_back( h ); // or newtables_id[h]
//         }
//       }
//       if(new_cluster_flag){
//         new_clusters.push_back(cluster_newtables[h]);
//         in_new_clusters.push_back(std::vector<int>(1, h )); // or newtables_id[h]
//       }
//     }
//   }

//   in_old_clusters.resize(orig_K);
//   for(int h = 0; h < num_splits - 1; h++){
//     if(cluster_newtables[h] < orig_K){
//       old_clusters[ cluster_newtables[h] ] += 1;
//       in_old_clusters[ cluster_newtables[h] ].push_back( h );
//     }
//   }
//   // the problem about this is that if I have an old or new cluster, I don't know where it is.
//   // it could be better to have something that is as long as the number of old clusters


//   int* orig_cluster_config = new int[orig_K];
//   int** orig_clusters = new int*[orig_K];
//   int* orig_cluster_assignment = new int[orig_nObs];
//   double* orig_log_prior = new double[orig_K];
//   for(int k = 0; k < orig_K; k++){
//     orig_cluster_config[k] = cluster_config[k];
//     orig_clusters[k] = new int[cluster_config[k]];
//     for(int i = 0; i < cluster_config[k]; i++){
//       orig_clusters[k][i] = clusters[k][i];
//     }
//     orig_log_prior[k] = log_prior[k];
//   }
//   for(int i = 0; i < orig_nObs; i++){
//     orig_cluster_assignment[i] = cluster_assignment[i];
//   }

//   delete[] cluster_config; cluster_config = NULL;
//   for(int k = 0; k < orig_K; k++){
//     delete[] clusters[k]; clusters[k] = NULL;
//   }
//   delete[] clusters; clusters = NULL;
//   delete[] cluster_assignment; cluster_assignment = NULL;
//   delete[] log_prior; log_prior = NULL;
//   for(int i = 0; i < orig_nObs; i++){
//     delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
//   }
//   delete[] pairwise_assignment; pairwise_assignment = NULL;

//   K += new_clusters.size();
//   cluster_config = new int[K];
//   clusters = new int*[K];
//   cluster_assignment = new int[nObs];
//   log_prior = new double[K];

//   for(int k = 0; k < K; k++){
//     if(k < orig_K){
//       cluster_config[k] = orig_cluster_config[k];
//     } else {
//       cluster_config[k] = 0;
//     }
//   }
//   int kk;
//   for(int k = 0; k < num_splits - 1; k++){
//     kk = cluster_newtables[k];
//     if(kk >= K){
//       cout << "Partition::AddTables_TableDish ERROR, invalid cluster_newtables" << endl;
//       return;
//       // this test should be placed at the beginning
//     }
//     cluster_config[kk] += 1;
//   }

//   std::vector<int> ys(num_splits - 1);
//   for(int k = 0; k < num_splits - 1; k++){ // ys should be in order
//     ys[k] = Table_RestId_MainId(r, newtables_id[k], nTableRest);
//     // cout << k << " " << ys[k] << endl;
//   }
//   // cout << "after ys" << endl;
//   int y, yy;
//   for(int k = 0; k < K; k++){
//     // cout << "allocating clusters "<< k << endl; 
//     clusters[k] = new int[cluster_config[k]];
//     if(k < orig_K){
//       int i = 0;
//       int ior = 0;
//       bool orig_done = false;
//       yy = orig_clusters[k][ior];
//       if(yy >= reindexing_threshold){
//         yy += reindexing_shift;
//       }
//       ior += 1;
//       // cout << "yy " << yy << endl;
//       for(int j = 0; j < in_old_clusters[k].size(); j++){
//         y = ys[ in_old_clusters[k][j] ];
//         // cout << "y " << y << endl;
//         while((yy < y) & (!orig_done)){ // check i < orig_cluster_config[k]
//           // if done = true, tutto e' falso
//           // if done = false, continuo normale
//           clusters[k][i] = yy;
//           // cout << i << " " << yy << endl;
//           i += 1;
//           if(ior < orig_cluster_config[k]){
//             yy = orig_clusters[k][ior];
//             if(yy >= reindexing_threshold){
//               yy += reindexing_shift;
//             }
//             ior += 1;
//             // cout << "yy " << yy << endl;
//           } else {
//             orig_done = true;
//           }
//         }
//         clusters[k][i] = y;
//         // cout << i << " " << y << endl;
//         i += 1;
//       }
//       while(!orig_done){
//         clusters[k][i] = yy; 
//         // cout << i << " " << yy << endl;
//         i += 1;
//         if(ior < orig_cluster_config[k]){
//           yy = orig_clusters[k][ior];
//           if(yy >= reindexing_threshold){
//             yy += reindexing_shift;
//           }
//           ior += 1;
//           // cout << "yy " << yy << endl;
//         } else {
//           orig_done = true;
//         } 
//       }
//     }
//     if(k >= orig_K){
//       // cout << "k >= orig_K" << endl;
//       for(int h = 0; h < new_clusters.size(); h++){
//         // cout << "new_clusters[h] " << new_clusters[h] << endl;
//         if( new_clusters[h] == k ){
//           // cout << "if( new_clusters[h] == k )" << in_new_clusters[h].size() << endl;
//           for(int j = 0; j < in_new_clusters[h].size(); j++){
//             // cout << "for loop j " << j << ";";
//             // cout << in_new_clusters[h][j] << " ";
//             y = ys[ in_new_clusters[h][j] ];
//             // cout << h << " " << new_clusters[h] << " " << j << " " << y << endl;
//             clusters[k][j] = y;
//             // cout << y << " ";
//           }
//           // cout << endl;
//         }
//       }
//     }
//   }
  
//   for(int k = 0; k < num_splits - 1; k++){
//     cluster_assignment[ys[k]] = cluster_newtables[k];
//   }
//   for(int yy = 0; yy < orig_nObs; yy++){
//     if(yy < reindexing_threshold){
//       cluster_assignment[yy] = orig_cluster_assignment[yy];
//     } else {
//       cluster_assignment[yy + reindexing_shift] = orig_cluster_assignment[yy];
//     }
//   }
  
//   get_pairwise();

//   for(int k = 0; k < K; k++){
//     if(k < orig_K){
//       if(old_clusters[k] == 0){
//         log_prior[k] = orig_log_prior[k];
//       } else {
//         log_pi_ep(k);
//       }
//     } else {
//       log_pi_ep(k);
//     }
//   }
//   delete[] orig_cluster_config;
//   for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length orig_K (K and not K+1.)
//     delete[] orig_clusters[kk];
//   }
//   delete[] orig_clusters;
//   delete[] orig_cluster_assignment;
//   delete[] orig_log_prior;
//   return;
// }

void Partition::RemoveTable_TableDish(int r, int oldtable_mainid, double eta_TD){
  /*
  RemoveTable_TableDish is used by Merge_CostTable (which is called by get_tables_merge)
    in that case the tables merged will belong to the same dish, so it leads to no empty dishes
  it is also used by ProposalSM_CostTable where instead any two tables within a restaurant can be merged
    so we need to be careful that we are not left with empty dishes
  */
  // we need to do reindexing
  // the number of elements of the partition is changed so we need to change a lot of other things
  int orig_nObs = nObs;
  nObs = orig_nObs-1;

  int orig_K = K;
  int* orig_cluster_config = new int[orig_K];
  int** orig_clusters = new int*[orig_K];
  int* orig_cluster_assignment = new int[orig_nObs];
  double* orig_log_prior = new double[orig_K];
  for(int k = 0; k < orig_K; k++){
    orig_cluster_config[k] = cluster_config[k];
    orig_clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      orig_clusters[k][i] = clusters[k][i];
    }
    orig_log_prior[k] = log_prior[k];
  }
  for(int i = 0; i < orig_nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }
  
  delete[] cluster_config; cluster_config = NULL;
  for(int k = 0; k < orig_K; k++){
    delete[] clusters[k]; clusters[k] = NULL;
  }
  delete[] clusters; clusters = NULL;
  delete[] cluster_assignment; cluster_assignment = NULL;
  delete[] log_prior; log_prior = NULL;
  for(int i = 0; i < orig_nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;

  int oldtable_cluster = orig_cluster_assignment[oldtable_mainid];
  if(orig_cluster_config[oldtable_cluster] == 1){
    // then we have one less cluster
    K = orig_K - 1;
    cluster_config = new int[K];
    clusters = new int*[K];
    cluster_assignment = new int[nObs];
    log_prior = new double[K];
    for(int k = 0; k < K; k++){
      if(k < oldtable_cluster){
        cluster_config[k] = orig_cluster_config[k];  
      } else {
        cluster_config[k] = orig_cluster_config[k+1];  
      }
    }
    int yy;
    for(int k = 0; k < K; k++){
      clusters[k] = new int[cluster_config[k]];
      if(k < oldtable_cluster){
        for(int i = 0; i < cluster_config[k]; i++){
          yy = orig_clusters[k][i];
          if(yy < oldtable_mainid){
            clusters[k][i] = yy;
            cluster_assignment[yy] = k;
          } else {
            clusters[k][i] = yy-1;
            cluster_assignment[yy-1] = k;
          }
        }
      } else{
        for(int i = 0; i < cluster_config[k]; i++){
          yy = orig_clusters[k+1][i];
          if(yy < oldtable_mainid){
            clusters[k][i] = yy;
            cluster_assignment[yy] = k;
          } else {
            clusters[k][i] = yy-1;
            cluster_assignment[yy-1] = k;
          }
        }
      }
    }
    get_pairwise();
    for(int k = 0; k < K; k++){
      if(k < oldtable_cluster){ // need to re-compute
        log_prior[k] = orig_log_prior[k];
      } else {
        log_prior[k] = orig_log_prior[k+1];
      }
    }
  } else {
    // the number of clusters stays the same
    K = orig_K;
    cluster_config = new int[K];
    clusters = new int*[K];
    cluster_assignment = new int[nObs];
    log_prior = new double[K];
    for(int k = 0; k < K; k++){
      if(k < oldtable_cluster){
        cluster_config[k] = orig_cluster_config[k];  
      } else if(k == oldtable_cluster){
        cluster_config[k] = orig_cluster_config[k]-1; 
      } else {
        cluster_config[k] = orig_cluster_config[k];  
      }
    }
    int yy;
    for(int k = 0; k < K; k++){
      clusters[k] = new int[cluster_config[k]];
      if(k < oldtable_cluster){
        for(int i = 0; i < cluster_config[k]; i++){
          yy = orig_clusters[k][i];
          if(yy < oldtable_mainid){
            clusters[k][i] = yy;
            cluster_assignment[yy] = k;
          } else {
            clusters[k][i] = yy-1;
            cluster_assignment[yy-1] = k;
          }
        }
      } else if(k == oldtable_cluster){
        int count = 0;
        for(int i = 0; i < orig_cluster_config[k]; i++){
          yy = orig_clusters[k][i];
          // if yy == oldtable_mainid we don't do anything
          if(yy < oldtable_mainid){
            clusters[k][count] = yy;
            count++;
            cluster_assignment[yy] = k;
          } else if(yy > oldtable_mainid){
            clusters[k][count] = yy-1; 
            count++;
            cluster_assignment[yy-1] = k;
          }
        }
      } else {
        for(int i = 0; i < cluster_config[k]; i++){
          yy = orig_clusters[k][i];
          if(yy < oldtable_mainid){
            clusters[k][i] = yy;
            cluster_assignment[yy] = k;
          } else{
            clusters[k][i] = yy-1;
            cluster_assignment[yy-1] = k;
          }
        }
      }
    }
    get_pairwise();
    for(int k = 0; k < K; k++){
      if(k < oldtable_cluster){ // need to re-compute
        log_prior[k] = orig_log_prior[k];
      } else if(k == oldtable_cluster){
        log_pi_ep(k, eta_TD);
      } else {
        log_prior[k] = orig_log_prior[k];
      }
    }
  }
  indexing.resize(nObs);
  delete[] orig_cluster_config;
  for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length orig_K (K and not K+1.)
    delete[] orig_clusters[kk];
  }
  delete[] orig_clusters;
  delete[] orig_cluster_assignment;
  delete[] orig_log_prior;
}
// Function will split the cluster into two parts
// new cluster1 will still be called cluster split_k
// new cluster2 will be called cluster K+1
void Partition::Split(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, double eta){
  // create a temporary copy of the main attributes
  int orig_K = K;
  int* orig_cluster_config = new int[orig_K];
  int** orig_clusters = new int*[orig_K];
  int* orig_cluster_assignment = new int[nObs];
  double* orig_log_prior = new double[orig_K];
  // double* orig_log_det_Omegay = new double[orig_K];
  // double* orig_quad_forms = new double[orig_K];
  
  for(int k = 0; k < orig_K; k++){
    orig_cluster_config[k] = cluster_config[k];
    orig_clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      orig_clusters[k][i] = clusters[k][i];
    }
    orig_log_prior[k] = log_prior[k];
    // orig_log_det_Omegay[k] = log_det_Omegay[k];
    // orig_quad_forms[k] = quad_forms[k];
  }
  for(int i = 0; i < nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }
  // clear up the memory from the original values and re-initialize
  // no need to delete pairwise_assignment or cluster_assignment; the sizes are fixed
  K = orig_K+1; // we now have K+1 clusters
  delete[] cluster_config; cluster_config = NULL;
  cluster_config = new int[K];

  for(int k = 0; k < K-1; k++){
    delete[] clusters[k]; clusters[k] = NULL;
  }
  delete[] clusters; clusters = NULL;
  clusters = new int*[K];
  delete[] log_prior; log_prior = NULL;
  // delete[] quad_forms; quad_forms = NULL;
  // delete[] log_det_Omegay; log_det_Omegay = NULL;
  log_prior = new double[K];
  // quad_forms = new double[K];
  // log_det_Omegay = new double[K];
  for(int i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;

  for(int k = 0; k < K; k++){
    if(k == split_k){
      cluster_config[k] = size1;
    } else if(k == K-1){
      cluster_config[k] = size2;
    } else{
      cluster_config[k] = orig_cluster_config[k];
    }
  }
  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
    if(k == split_k){
      for(int i = 0; i < size1; i++){
        clusters[k][i] = new_cluster1[i];
      }

    } else if(k == K - 1){ // remember the 0-indexing...
      for(int i = 0; i < size2;i++){
        clusters[k][i] = new_cluster2[i];
      }
    } else{
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = orig_clusters[k][i];
      }
    }
  }

  for(int i = 0; i < nObs; i++){
    if(orig_cluster_assignment[i] != split_k){
      cluster_assignment[i] = orig_cluster_assignment[i];
    }
  }
  for(int ii = 0; ii < size1; ii++){
    cluster_assignment[new_cluster1[ii]] = split_k;
  }
  for(int ii = 0; ii < size2; ii++){
    cluster_assignment[new_cluster2[ii]] = K - 1; // remember, cluster labels go from 0 to K-1.
  }
  get_pairwise();
  for(int k = 0; k < K; k++){
    if(k == split_k){ // need to re-compute
      log_pi_ep(k, eta);
      // get_likelihood(A_or_B, k);
    } else if(k == K-1){
      log_pi_ep(k, eta);
      // get_likelihood(A_or_B, k);
    } else{
      log_prior[k] = orig_log_prior[k];
      // quad_forms[k] = orig_quad_forms[k];
      // log_det_Omegay[k] = orig_log_det_Omegay[k];
    }
  }
  // free up memory by deleting the local copies
  delete[] orig_cluster_config;
  for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length orig_K (K and not K+1.)
    delete[] orig_clusters[kk];
  }
  delete[] orig_clusters;
  delete[] orig_cluster_assignment;
  delete[] orig_log_prior;
  // delete[] orig_quad_forms;
  // delete[] orig_log_det_Omegay;
  return;
}

void Partition::KSplit(int split_k, int num_splits, std::vector<std::vector<int> > indices, std::vector<int> ns, double eta){
  int orig_K = K;
  int* orig_cluster_config = new int[orig_K];
  int** orig_clusters = new int*[orig_K];
  int* orig_cluster_assignment = new int[nObs];
  double* orig_log_prior = new double[orig_K];
  // double* orig_log_det_Omegay = new double[orig_K];
  // double* orig_quad_forms = new double[orig_K];

  for(int k = 0; k < orig_K; k++){
    orig_cluster_config[k] = cluster_config[k];
    orig_clusters[k] = new int[cluster_config[k]];
    for(int i = 0; i < cluster_config[k]; i++){
      orig_clusters[k][i] = clusters[k][i];
    }
    orig_log_prior[k] = log_prior[k];
    // orig_log_det_Omegay[k] = log_det_Omegay[k];
    // orig_quad_forms[k] = quad_forms[k];
  }
  for(int i = 0; i < nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }

  /* Now we free the space */
  delete[] cluster_config; cluster_config = NULL;
  for(int i = 0; i < orig_K; i++){
    delete[] clusters[i]; clusters[i] = NULL;
  }
  delete[] clusters; clusters = NULL;
  delete[] log_prior; log_prior = NULL;
  // delete[] quad_forms; quad_forms = NULL;
  // delete[] log_det_Omegay; log_det_Omegay = NULL;
  for(int i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;
  // no need to delete cluster_assignment; the sizes are fixed

  /* Here we update each of them */
  K = orig_K - 1 + num_splits;
  cluster_config = new int[K];
  clusters = new int*[K];
  for(int k = 0; k < K; k++){
    if(k == split_k){
      cluster_config[k] = ns[0];
    } else if(k < orig_K){
      cluster_config[k] = orig_cluster_config[k];
    } else{
      cluster_config[k] = ns[k - orig_K + 1];
    }
  }
  for(int k = 0; k < K; k++){
    clusters[k] = new int[cluster_config[k]];
    if(k == split_k){
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = indices[0][i];
      }
    } else if(k < orig_K){
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = orig_clusters[k][i];
      }
    } else{
      for(int i = 0; i < cluster_config[k]; i++){
        clusters[k][i] = indices[k - orig_K + 1][i];
      }
    }
  }
  for(int i = 0; i < nObs; i++){
    if(orig_cluster_assignment[i] != split_k){
      cluster_assignment[i] = orig_cluster_assignment[i];
    }
  }
  for(int kk = 0; kk < num_splits; kk++){
    if(kk == 0){
      for(int i = 0; i < ns[kk]; i++){
        cluster_assignment[ indices[kk][i] ] = split_k;
      }
    } else {
      for(int i = 0; i < ns[kk]; i++){
        cluster_assignment[ indices[kk][i] ] = kk + orig_K - 1;
      }  
    }
  }
  get_pairwise();
  log_prior = new double[K];
  // quad_forms = new double[K];
  // log_det_Omegay = new double[K];
  for(int k = 0; k < K; k++){
    if(k == split_k){ // need to re-compute
      log_pi_ep(k, eta);
      // get_likelihood(A_or_B, k);
    } else if(k < orig_K){
      log_prior[k] = orig_log_prior[k];
      // quad_forms[k] = orig_quad_forms[k];
      // log_det_Omegay[k] = orig_log_det_Omegay[k];
    } else{
      log_pi_ep(k, eta);
      // get_likelihood(A_or_B, k);
    }
  }
  
  // if(opt_method == 0){
  //   // just do what we were doing before
  //   for(int k = 0; k < K; k++){
  //     if(A_or_B){
  //       get_parameter(A_or_B, k);
  //     } else {
  //       get_parameter(A_or_B, k);
  //     }
  //   }
  // } else if(opt_method == 1){
  //   // optimize alpha_beta given the two partitions
  //   get_alpha_beta();
  // } else if(opt_method == 2){
  //   // optimize alpha_beta given the two partitions, by doing coordinate ascend
  //   get_alpha_beta_iter();
  // }
  delete[] orig_cluster_config;
  for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length orig_K (K and not K+1.)
    delete[] orig_clusters[kk];
  }
  delete[] orig_clusters;
  delete[] orig_cluster_assignment;
  delete[] orig_log_prior;
  // delete[] orig_quad_forms;
  // delete[] orig_log_det_Omegay;
  return;
}

void Partition::Merge(int k_1, int k_2, double eta){
  // k1, k2 get merged into k_min
  int k_max = max(k_1, k_2);
  int k_min = min(k_1, k_2);
  int new_cluster_size = cluster_config[k_min] + cluster_config[k_max];

    // make a pointer to the new merged cluster
  int* new_merged_cluster = new int[new_cluster_size];
  for(int i = 0; i < cluster_config[k_min]; i++){
    new_merged_cluster[i] = clusters[k_min][i];
  }
  for(int i = 0; i < cluster_config[k_max]; i++){
    new_merged_cluster[cluster_config[k_min] + i] = clusters[k_max][i];
  }
  // Update cluster_assignment
  // for original cluster k_max: this now becomes k_min
  // for clusters with original label greater than k_max, we need to decrement by 1
  int tmp_assignment = 0;
  for(int i = 0; i < nObs; i++){
    if(cluster_assignment[i] > k_max){ // we decrement cluster label by 1
      tmp_assignment = cluster_assignment[i];
      cluster_assignment[i] = tmp_assignment - 1;
    } else if(cluster_assignment[i] == k_max){
      cluster_assignment[i] = k_min;
    }
  }
  // make a temporary copy of clusters and cluster_config
  int orig_K = K;
  int* orig_cluster_config = new int[orig_K];
  int** orig_clusters = new int*[orig_K];
  int* orig_cluster_assignment = new int[nObs];
  double* orig_log_prior = new double[orig_K];
  // double* orig_log_det_Omegay = new double[orig_K];
  // double* orig_quad_forms = new double[orig_K];
  
  for(int kk = 0; kk < orig_K; kk++){
    orig_cluster_config[kk] = cluster_config[kk];
    orig_clusters[kk] = new int[cluster_config[kk]];
    for(int i = 0; i < cluster_config[kk]; i++){
      orig_clusters[kk][i] = clusters[kk][i];
    }
    orig_log_prior[kk] = log_prior[kk];
    // orig_log_det_Omegay[kk] = log_det_Omegay[kk];
    // orig_quad_forms[kk] = quad_forms[kk];
  }
  for(int i = 0; i < nObs; i++){
    orig_cluster_assignment[i] = cluster_assignment[i];
  }
  delete[] cluster_config;
  cluster_config = new int[K-1];

  for(int kk = 0; kk < orig_K; kk++){
    delete[] clusters[kk];
  }
  delete[] clusters;
  clusters = new int*[K-1];
  delete[] log_prior;
  log_prior = new double[K-1];
  // delete[] quad_forms; quad_forms = NULL;
  // delete[] log_det_Omegay; log_det_Omegay = NULL;
  // quad_forms = new double[K-1];
  // log_det_Omegay = new double[K-1];
  // looping over the OLD labels
  // remember the labels don't change until we get to k_max
  // this loop visits every cluster EXCEPT k_max
  for(int kk = 0; kk < orig_K; kk++){
    if(kk == k_min){
      cluster_config[kk] = new_cluster_size;
      clusters[kk] = new int[cluster_config[kk]];
      for(int i = 0; i < cluster_config[kk]; i++){
        clusters[kk][i] = new_merged_cluster[i];
      }
      log_pi_ep(kk, eta);
      // get_likelihood(A_or_B, kk);
    } else if(kk < k_max){
      cluster_config[kk] = orig_cluster_config[kk];
      clusters[kk] = new int[cluster_config[kk]];
      for(int i = 0; i < cluster_config[kk]; i++){
        clusters[kk][i] = orig_clusters[kk][i];
      }
      log_prior[kk] = orig_log_prior[kk];
      // quad_forms[kk] = orig_quad_forms[kk];
      // log_det_Omegay[kk] = orig_log_det_Omegay[kk];
    } else if(kk > k_max){
      cluster_config[kk-1] = orig_cluster_config[kk];
      clusters[kk-1] = new int[cluster_config[kk-1]];// ceci changed cluster_config[kk] into cluster_config[kk-1]
      for(int i = 0; i < cluster_config[kk-1]; i++){ // ceci changed cluster_config[kk] into cluster_config[kk-1]
        clusters[kk-1][i] = orig_clusters[kk][i];
      }
      log_prior[kk-1] = orig_log_prior[kk];
      // quad_forms[kk-1] = orig_quad_forms[kk];
      // log_det_Omegay[kk-1] = orig_log_det_Omegay[kk];
    }
  }
  for(int i = 0; i < nObs; i++){
    delete[] pairwise_assignment[i]; pairwise_assignment[i] = NULL;
  }
  delete[] pairwise_assignment; pairwise_assignment = NULL;
  
  // update K
  K = orig_K - 1;

  get_pairwise();

  
  // clean-up the memory
  delete[] orig_cluster_config;
  for(int kk = 0; kk < orig_K; kk++){ //original_clusters has length K and not K+1.
    delete[] orig_clusters[kk];
  }
  delete[] orig_clusters;
  delete[] orig_cluster_assignment;
  delete[] orig_log_prior;
  // delete[] orig_quad_forms;
  // delete[] orig_log_det_Omegay;
  delete[] new_merged_cluster;
}

void Partition::Split_and_Merge(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int k_star_1, int k_star_2, double eta){
  // the idea is that newcluster1 is merged with kstar1 and newcluster2 is merged with kstar2.
  if(k_star_1 != k_star_2){
    Split(split_k, new_cluster1, new_cluster2, size1, size2, eta);
    if((split_k == k_star_1) & (split_k != k_star_2)){ // leave new_cluster1 alone and just attempt to merge new_cluster2 with k_star_2
      Merge(K-1, k_star_2, eta); // remember new_cluster2's label is K-1
    } else if( (split_k != k_star_1) & (split_k == k_star_2)){
      // leave new_clsuter2 by itself and merge k_star_1 with new_cluster1
      // remember new_cluster1's label is still split_k;
      Merge(split_k, k_star_1, eta);
    } else if((split_k != k_star_1) & (split_k != k_star_2) & (k_star_2 > max(split_k, k_star_1))){
      // We need to perform two merges
      Merge(split_k, k_star_1, eta);
      // k_star_2's label is now decremented by 1
      // new_cluster2's label is now K
      Merge(K, k_star_2 - 1, eta);
    } else if((split_k != k_star_1) & (split_k != k_star_2) & (k_star_2 < max(split_k, k_star_1))){
      Merge(split_k, k_star_1, eta);
      Merge(K, k_star_2, eta);
    }
  }
  return;
}

// Function to modify a partition if any of the clusters are disconnected
// void Partition::Modify(int cl_ind){
//   // cl_ind is the index of the cluster that needs to be modified
//   if(cluster_config[cl_ind] > 1)
//   {
//     int n = cluster_config[cl_ind];
//     int *index = clusters[cl_ind];
//     mat A_cluster = Submatrix(A_block, n, n, index, index);
//     int *components = new int[n];
//     int *count = new int;
//     Connected_Components(A_cluster, n, components, count);
//     if(*count > 1)
//     {
//       //I will use split iteratively, first split 0 from !=0
//       //then among the remaining, split 1 from !=1, etc
//       *count = *count - 1;
//       int *new_components, *index1, *index2, *i_ind, n1, n2;
//       for(int tmp = 0; tmp < *count; tmp++)
//       {
//         index1 = new int[n];
//         index2 = new int[n];
//         i_ind = new int[n];
//         n1 = 0;
//         n2 = 0;
//         for(int i = 0; i < n; i++){
//           if(components[i] == tmp){
//             index1[n1] = index[i];
//             n1++;
//           } else {
//             index2[n2] = index[i];
//             i_ind[n2] = i;
//             n2++;
//           }
//         }
//         Split(cl_ind, index1, index2, n1, n2);
//         if(tmp > 0)
//           delete[] index;
//         index = index2;
//         n = n2;
//         new_components = new int[n2];
//         for(int j = 0; j < n2; j++)
//           new_components[j] = components[i_ind[j]];
//         delete[] components;
//         components = new_components;
//         cl_ind = K-1;
//         delete[] index1;
//       }
//       delete[] index2;
//     }
//     delete[] components;
//     delete count;
//   }
// }
// 


std::vector<int> Partition::SampleSM(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR){
  int pars_size = 4;
  std::vector<int> pars(pars_size);
  std::vector<int> ret(pars_size+1);

  LPPartition proposal = new Partition(this);
  double lqratio, lpratio, llikratio;
  proposal->ProposalSM(this, ngibbs, lqratio, lpratio, llikratio, eng, pars, eta_LR);
  double a = lqratio + temp * (lpratio + llikratio);
  
  // WHAT ARE PARS? -- what do we need them for?
  // pars[0] = i1;
  // pars[1] = i2;
  // pars[2] = c1;
  // pars[3] = c2;

  a = exp(a);
  a = min(1.0,a);
  uniform_real_distribution<double> distribution(0.0,1.0);
  if(distribution(eng) < a){
    Copy_Partition(proposal);
    // returned_value = true;
    ret[0] = 1;
  } else {
    // returned_value = false;
    ret[0] = 0;
  }
  for(int i = 0; i < pars_size; i++){
    ret[i+1] = pars[i];
  }
  delete proposal;
  // return returned_value;
  return ret;
} 

void Partition::ProposalSM(LPPartition startingP, int ngibbs, double& lqratio, double& lpratio, double& llikratio, 
                           std::default_random_engine& myRandomEngine, std::vector<int>& pars, double eta_LR){
  // Should modify proposal, not startingP, update the values of lqratio, lpratio, llikratio
  // Create a random device and use it to generate a random seed
  // std::random_device myRandomDevice;
  // unsigned seed = myRandomDevice();
  // std::default_random_engine myRandomEngine(seed);
  std::vector<int> i12 = sample_without_replacement(2, startingP->nObs, myRandomEngine);
  int i1 = i12[0];
  int i2 = i12[1];
  int c1 = startingP->cluster_assignment[i1];
  int c2 = startingP->cluster_assignment[i2];
  if(c2 < c1){
    // c1 should be smaller than c2
    // change the labels i1 and i2 so that c1 < c2
    int tmp = i1;
    i1 = i2;
    i2 = tmp;
    tmp = c1;
    c1 = c2;
    c2 = tmp;
  } 
  if((c1 == c2) && (i2 < i1)){
    int tmp = i1;
    i1 = i2;
    i2 = tmp;
  }
  pars[0] = i1;
  pars[1] = i2;
  pars[2] = c1;
  pars[3] = c2;
  // samecl_index contains the other people that belong to c1 or c2 (but not i1, i2) localID
  // if c1 == c2, then it's the other people in that cluster
  // if c1 != c2, then it's the other people in the two clusters
  std::vector<int> samecl_index = startingP->get_samecl_index(i1, i2, c1, c2);
  
  // launch will NOT modify the current partition, but update new_c1_ind, new_c2_ind,samecl_newcl
  std::list<int> new_c1_ind,new_c2_ind; // lists are optimized to remove elements in any position (use "remove")
  std::vector<int> samecl_newcl;
  
  LaunchSM(startingP, i1, i2, samecl_index, c1, c2, new_c1_ind, new_c2_ind, samecl_newcl, ngibbs, myRandomEngine);
  if(c1 == c2){
    // need to compute llikratio using a function that computes the log likelihood of one cluster 
    // everyone together
    llikratio = -startingP->cluster_llikelihood(c1, true);
    lqratio = 0;
    // cluster labels of the two new clusters
    int new_c1 = c1;
    int new_c2 = startingP->K;
    // int nc1 = 0;
    // int nc2 = 0;
    ProposalSM_iter(true, startingP, samecl_index, new_c1, new_c2, new_c1_ind, new_c2_ind,samecl_newcl, lqratio, myRandomEngine);
    int nc1 = new_c1_ind.size();
    int nc2 = new_c2_ind.size();
    lpratio = log(eta_LR) + lbeta(nc1,nc2);
    // update llikratio with the cluster a, and cluster b
    llikratio += startingP->cluster_llikelihood(new_c1_ind, true) + startingP->cluster_llikelihood(new_c2_ind, true);

    // NOW WE UPDATE THE CURRENT PARTITION, BY SPLITTING
    new_c1_ind.sort();
    new_c2_ind.sort();
    int *index1, *index2;
    index1 = new int[nc1];
    index2 = new int[nc2];
    for(int i = 0; i < nc1; i++){
      // if this does not work we can use the iterator as in cluster_llikelihood
      index1[i] = new_c1_ind.front();
      new_c1_ind.pop_front();
    }
    for(int i = 0; i < nc2; i++){
      index2[i] = new_c2_ind.front();
      new_c2_ind.pop_front();
    }
    Split(c1, index1, index2, nc1, nc2, eta_LR);
    delete[] index1;
    delete[] index2;
  } else {
    // maybe at the end: merge everyone into one cluster
    lpratio = - log(eta_LR) - lbeta(startingP->cluster_config[c1], startingP->cluster_config[c2]);
    lqratio = 0;
    
    llikratio = - startingP->cluster_llikelihood(c1, true) - startingP->cluster_llikelihood(c2, true);
    
    ProposalSM_iter(false, startingP, samecl_index, c1, c2, new_c1_ind, new_c2_ind,samecl_newcl, lqratio, myRandomEngine);
    // MERGE EVERYONE INTO ONE CLUSTER
    Merge(c1,c2, eta_LR);
    llikratio += cluster_llikelihood(c1, true);
  }
  return;
}

void Partition::ProposalSM_iter(bool SplitMerge, LPPartition startingP, std::vector<int>& samecl_index, int new_c1, int new_c2,
  std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, double& lqratio, std::default_random_engine& myRandomEngine){
  // does not modify the current partition
  arma::vec P(2);
  arma::vec l(2);
  arma::vec logpr, pr;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  for(int i = 0; i < samecl_index.size(); i++){
    if(samecl_newcl[i] == new_c1){
      P[0] = new_c1_ind.size()-1;
      P[1] = new_c2_ind.size();
      l[0] = startingP->cluster_llikelihood(new_c1_ind, true) + startingP->cluster_llikelihood(new_c2_ind, true);
      new_c1_ind.remove(samecl_index[i]);
      new_c2_ind.push_back(samecl_index[i]);
      // update nTableRest (does not change)
      l[1] = startingP->cluster_llikelihood(new_c1_ind, true) + startingP->cluster_llikelihood(new_c2_ind, true);
      new_c2_ind.pop_back();
      new_c1_ind.push_back(samecl_index[i]);
    } else {
      P[0] = new_c1_ind.size();
      P[1] = new_c2_ind.size()-1;
      l[1] = startingP->cluster_llikelihood(new_c1_ind, true) + startingP->cluster_llikelihood(new_c2_ind, true);
      new_c2_ind.remove(samecl_index[i]);
      new_c1_ind.push_back(samecl_index[i]);
      l[0] = startingP->cluster_llikelihood(new_c1_ind, true) + startingP->cluster_llikelihood(new_c2_ind, true);
      new_c1_ind.pop_back();
      new_c2_ind.push_back(samecl_index[i]);
    }
    logpr = log(P)+l;
    logpr = logpr - max(logpr);
    logpr = logpr - log(sum(exp(logpr)));
    pr = exp(logpr);
    // update the membership of samecl_index[i] in the current partition
    if(SplitMerge){
      if(distribution(myRandomEngine) < pr[0]){ // do 0
        if(samecl_newcl[i] == new_c1){
          // then we do noting
        } else {
          new_c2_ind.pop_back();
          new_c1_ind.push_back(samecl_index[i]);
          samecl_newcl[i] = new_c1;
        }
        lqratio += -logpr[0];
      } else { // do 1
        if(samecl_newcl[i] == new_c1){
          new_c1_ind.pop_back();
          new_c2_ind.push_back(samecl_index[i]);
          samecl_newcl[i] = new_c2;
        } else {
          // then we do noting
        }
        lqratio += -logpr[1];
      } // end sampling of new_c1 or new_c2
    } else {
      if(startingP->cluster_assignment[ samecl_index[i] ] == new_c1){
        lqratio += logpr[0];
      } else if(startingP->cluster_assignment[ samecl_index[i] ] == new_c2) {
        lqratio += logpr[1];
      } else {
        cout << "Partition::ProposalSM_iter Problem!";
        Print_Partition();
      }
    }
  }
  return;
}


// need to write the function that computes the likelihood for a cluster
// maybe it's easier if it takes a int* instead of vector

// in the split case, it makes sense to change the current partition (proposal) because that will be furtherly updated
// in the merge case, it does not make sense to change the partition, because then it will be merged again.
// So probably instead of modifying the partititon we can return the two lists.

void Partition::LaunchSM(LPPartition startingP, int i1, int i2, std::vector<int> samecl_index, 
  int c1, int c2, std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, 
  int ngibbs, std::default_random_engine& gen){
  // will NOT modify the current partition
  // samecl_index contain the indices (from 1 to n)
  // new_c1_ind, new_c2_ind contain the indices (from 1 to n) divided into the two clusters
  // samecl_newcl contains either new_c1 or new_c2 for the elements of samecl_index, with the same order

  // // Create a random device and use it to generate a random seed
  // std::random_device myRandomDevice;
  // unsigned seed = myRandomDevice();
  // // Initialize a default_random_engine with the seed
  // std::default_random_engine gen(seed);


  int new_c1, new_c2;
  if(c1 == c2){
    new_c1 = c1;
    new_c2 = startingP->K; // 0-indexing
  } else {
    new_c1 = c1;
    new_c2 = c2;
  }
  new_c1_ind.push_back(i1);
  new_c2_ind.push_back(i2);

  if(samecl_index.size() > 0){
    
    std::bernoulli_distribution bern(0.5);
    
    samecl_newcl.resize(samecl_index.size());
    // randomly assign each element of samecl into one cluster or the other
    for(int i = 0; i < samecl_index.size(); i++){
      if(bern(gen)){  // bern(gen) returns bool
        new_c1_ind.push_back(samecl_index[i]); // becomes new_c1
        samecl_newcl[i] = new_c1;
      } else {
        new_c2_ind.push_back(samecl_index[i]); // new_c2
        samecl_newcl[i] = new_c2;
      }
    }
    double lqratio_ignored = 0;
    for(int t = 0; t < ngibbs; t++){
      ProposalSM_iter(true, startingP, samecl_index, new_c1, new_c2, new_c1_ind, new_c2_ind, samecl_newcl, lqratio_ignored, gen);
    } // end loop t
    // we are done changing the lists
  } 
  return;
}


double Partition::cluster_llikelihood(std::list<int> index, bool LR){
  arma::mat Yloc;
  double alpha_0, beta_0, mu_0, k_0;
  if(LR){
    Yloc = YLR;
    alpha_0 = alpha_L;
    beta_0 = beta_L;
    mu_0 = mu_L;
    k_0 = k_L;
  } else {
    Yloc = Y;
    alpha_0 = alpha_H;
    beta_0 = beta_H;
    mu_0 = mu_H;
    k_0 = k_H;
  }
  // pretend this list was a cluster
  int n = index.size();
  arma::vec y_cl(n);
  int i = 0;
  for(std::list<int>::iterator it = index.begin(); it != index.end(); ++it){
    y_cl[i] = Yloc[ *it ];
    i++;
  }
  double ybar = arma::mean(y_cl);
  double k_n = k_0 + n;
  // double mu_n = (k_0*mu_0 + n*ybar)/k_n;
  double alpha_n = alpha_0 + n/2;
  double beta_n = beta_0 + sum(pow(y_cl-ybar,2)) + k_0 * n * pow(ybar-mu_0,2)/(2*k_n);
  double loglike = lgamma(alpha_n) - lgamma(alpha_0) + alpha_0 * log(beta_0) - alpha_n * log(beta_n);
  loglike += 0.5*(log(k_0)-log(k_n)) -n/2*log(2*M_PI);
  return(loglike);
}

double Partition::get_logposterior(double eta_LR){
  double res = 0.0;
  res += get_logprior(eta_LR);
  res += get_loglike(true);
  return res;
}

double Partition::cluster_llikelihood(int cl, bool LR){
  arma::mat Yloc;
  double alpha_0, beta_0, mu_0, k_0;
  if(LR){
    Yloc = YLR;
    alpha_0 = alpha_L;
    beta_0 = beta_L;
    mu_0 = mu_L;
    k_0 = k_L;
  } else {
    Yloc = Y;
    alpha_0 = alpha_H;
    beta_0 = beta_H;
    mu_0 = mu_H;
    k_0 = k_H;
  }
  // pretend this list was a cluster
  int n = cluster_config[cl];
  arma::vec y_cl(n);
  for(int i = 0; i < n; i++){
    y_cl[i] = Yloc[ clusters[cl][i] ];
  }
  double ybar = arma::mean(y_cl);
  double k_n = k_0 + n;
  // double mu_n = (k_0*mu_0 + n*ybar)/k_n;
  double alpha_n = alpha_0 + n/2;
  double beta_n = beta_0 + sum(pow(y_cl-ybar,2)) + k_0 * n * pow(ybar-mu_0,2)/(2*k_n);
  double loglike = lgamma(alpha_n) - lgamma(alpha_0) + alpha_0 * log(beta_0) - alpha_n * log(beta_n);
  loglike += 0.5*(log(k_0)-log(k_n)) -n/2*log(2*M_PI);
  return(loglike);
}
// I think the likelihood function should take the partition and a vector (not a pointer, or a cluster) of elements 
// now I am feeding it lists/
// Eventually we can just put this inside of update_particle

double Partition::get_loglike(bool LR){
  double res = 0.0;
  for(int k = 0; k < K; k++){
    res += cluster_llikelihood(k, LR);
  }
  return res;
}

std::vector<int> Partition::get_samecl_index(int i1, int i2, int c1, int c2){
  // create samecl_index, which contains the other people in those cluster (but not i1, i2)
  std::vector<int> samecl_index;
  int indtmp;
  if(c1 == c2){
    // If I resize, it will add 0's to the vector!
    samecl_index.reserve(cluster_config[c1]-2); 
    for(int i = 0; i < cluster_config[c1]; i++){
      indtmp = clusters[c1][i];
      if((indtmp != i1) && (indtmp != i2)){
        samecl_index.push_back(indtmp);
      }
    }
  } else {
    samecl_index.reserve(cluster_config[c1] + cluster_config[c2] - 2); 
    for(int i = 0; i < cluster_config[c1]; i++){
      indtmp = clusters[c1][i];
      if(indtmp != i1){
        samecl_index.push_back(indtmp);
      }
    }
    for(int i = 0; i < cluster_config[c2]; i++){
      indtmp = clusters[c2][i];
      if(indtmp != i2){
        samecl_index.push_back(indtmp);
      }
    }
  }
  return samecl_index;
}

bool Partition::FindKMSplit(int split_k, arma::vec y, double eta){
  int n_cl = cluster_config[split_k];
  arma::mat U = arma::zeros<mat>(n_cl, 1); // holds the data passed to k-means
  U.set_size(n_cl, 1);
  for(int i = 0; i < n_cl; i++){
    U(i,0)  = y[i];
  }
  arma::mat means;

  std::vector<std::vector<int> > indices(2);
  std::vector<int> ns(2);
  
  // double min_score;
  int * membership;
  bool status = false;
  while(!status){
    status = arma::kmeans(means, U.t(), 2, random_subset, 10, false);
    // status = kmeans_repeat(means, U, 2, 10, min_score);
    if(status == false){
      cout << "clustering failed" << endl;
      // go directly to next iteration: i.e. do kmeans again
    }
  }
  membership = which_is_nearest_k(means, U.t());

  // indices has the actual indices from 0 to nObs-1
  for(int i = 0; i < n_cl; i++){
    indices[membership[i]].push_back( clusters[split_k][i] );
    (ns[membership[i]])++;
  }
  delete[] membership;

  int size1 = indices[0].size();
  int size2 = indices[1].size();
  int* new_cluster1 = new int[size1];
  int* new_cluster2 = new int[size2];
  for(int i = 0; i < size1; i++){
    new_cluster1[i] = indices[0][i];
  }
  for(int i = 0; i < size2; i++){
    new_cluster2[i] = indices[1][i];
  }
  Split(split_k, new_cluster1, new_cluster2, size1, size2, eta);
  delete[] new_cluster1;
  delete[] new_cluster2;
  return true;
}

void Partition::Combine_Restaurants(LPPartition part1, LPPartition part2){
  /*
  Combine_Restaurants creates the merged partition of costumers_tables[rmin], costumers_tables[rmax]
  It considers the tables from the two partitions and simply puts them in the same restaurant
  without merging them. The number of tables of the new merged restaurant is K1 + K2
  The names of the tables for part2 are shifted and start by part1->K up to part2->K-1
  and are incorporated in the same order of the tables names (first tables 0..K1-1 then K1..(K1+K2-1))
  The names of the costumers of part2 are also shifted by part1->nObs
  */
  nObs = part1->nObs + part2->nObs;
  K = part1->K + part2->K;
  cluster_config = new int[K];
  clusters = new int*[K];
  for(int i = 0; i < part1->K; i++){
    cluster_config[i] = part1->cluster_config[i];
    clusters[i] = new int[part1->cluster_config[i]];
    for(int j = 0; j < part1->cluster_config[i]; j++){
      clusters[i][j] = part1->clusters[i][j];
    }
  }
  for(int i = 0; i < part2->K; i++){
    cluster_config[part1->K + i] = part2->cluster_config[i];
    clusters[part1->K + i] = new int[part2->cluster_config[i]];
    for(int j = 0; j < part2->cluster_config[i]; j++){
      clusters[part1->K + i][j] = part1->nObs + part2->clusters[i][j];
    }
  }
  cluster_assignment = new int[nObs];
  for(int i = 0; i < part1->nObs; i++){
    cluster_assignment[i] = part1->cluster_assignment[i];
  }
  for(int i = 0; i < part2->nObs; i++){
    cluster_assignment[part1->nObs + i] = part1->K + part2->cluster_assignment[i];
  }
  get_pairwise();
  log_prior = new double[K];
  for(int i = 0; i < part1->K; i++){
    log_prior[i] = part1->log_prior[i];
  }
  for(int i = 0; i < part2->K; i++){
    log_prior[part1->K + i] = part2->log_prior[i];
  }
  indexing.resize(nObs);
  for(int i = 0; i < part1->nObs; i++){
    indexing[i] = part1->indexing[i];
  }
  for(int i = 0; i < part2->nObs; i++){
    indexing[part1->nObs + i] = part2->indexing[i];
  }
  return;
}

void Partition::Divide_Restaurants(LPPartition newpart, std::vector<int> index1, std::vector<int> index2, double eta_CT, bool local_print){
  /*
  Modifies Partition which represents costumers_tables[r1] and creates newpart, which represents costumers_tables[r2].
  We change the observations and the clusters within each restaurant. We go through each index (which contains localID of the costumers in the two splits)
  and create clusters in the order of appearance in index1/2 (first cluster is the cluster of index1[0], second cluster is the first new cluster appearing in index1, etc)
  It is important that SplitTables_TableDish considers the same order when creating the tables (it does now - Jan 7, 2021)
  Then goes through index2 and creates newpart.
  */
  // if(local_print) cout << endl << "Divide_Restaurants" << endl;
  // updates this and newpart
  int nObs_orig = nObs;
  int K_orig = K;
  nObs = index1.size();
  newpart->nObs = index2.size();

  std::vector<int> scl1,scl2;
  std::vector<std::vector<int> > clusters1,clusters2;
  std::vector<int> cluster_assignment1(nObs);
  std::vector<int> cluster_assignment2(newpart->nObs);
  std::vector<int> indexing1;
  std::vector<int> indexing2;
  // std::pair<std::set<int>::iterator,bool> ret;
  int cl, ind_oldcl;
  bool newcl;
  for(int i = 0; i < index1.size(); i++){
    indexing1.push_back(indexing[index1[i]]);

    cl = cluster_assignment[ index1[i] ];
    newcl = true;
    for(int c = 0; c < scl1.size(); c++){
      if(scl1[c] == cl){
        newcl = false;
        ind_oldcl = c;
        break;
      }
    }
    if(newcl){
      scl1.push_back(cl);
      clusters1.push_back(std::vector<int>(1,i));
      cluster_assignment1[i] = scl1.size()-1;
    } else {
      clusters1[ind_oldcl].push_back(i);
      cluster_assignment1[i] = ind_oldcl;
    }
    // if(local_print){
    //   printf("i %d, index1[i] %d, indexing %d, in cluster %d which is (%d) new and saved as %d.\n",i,index1[i],indexing[index1[i]],cl,newcl,cluster_assignment1[i]);
    // }
  }
  for(int i = 0; i < index2.size(); i++){
    indexing2.push_back(indexing[index2[i]]);

    cl = cluster_assignment[ index2[i] ];
    newcl = true;
    for(int c = 0; c < scl2.size(); c++){
      if(scl2[c] == cl){
        newcl = false;
        ind_oldcl = c;
        break;
      }
    }
    if(newcl){
      scl2.push_back(cl);
      clusters2.push_back(std::vector<int>(1,i));
      cluster_assignment2[i] = scl2.size()-1;
    } else {
      clusters2[ind_oldcl].push_back(i);
      cluster_assignment2[i] = ind_oldcl;
    }
  }

  K = clusters1.size();
  newpart->K = clusters2.size();

  delete[] cluster_config;
  for(int i = 0; i < K_orig; i++){
    delete[] clusters[i];
  }
  delete[] clusters; 
  delete[] cluster_assignment;
  for(int i = 0; i < nObs_orig; i++){
    delete[] pairwise_assignment[i]; 
  }
  delete[] pairwise_assignment; 
  delete[] log_prior;

  cluster_config = new int[K];
  clusters = new int*[K];
  cluster_assignment = new int[nObs];
  log_prior = new double[K];

  newpart->cluster_config = new int[newpart->K];
  newpart->clusters = new int*[newpart->K];
  newpart->cluster_assignment = new int[newpart->nObs];
  newpart->log_prior = new double[newpart->K];

  for(int i = 0; i < clusters1.size(); i++){
    cluster_config[i] = clusters1[i].size();
    clusters[i] = new int[cluster_config[i]];
    for(int j = 0; j < clusters1[i].size(); j++){
      clusters[i][j] = clusters1[i][j];
    }
  }

  for(int i = 0; i < clusters2.size(); i++){
    newpart->cluster_config[i] = clusters2[i].size();
    newpart->clusters[i] = new int[newpart->cluster_config[i]];
    for(int j = 0; j < clusters2[i].size(); j++){
      newpart->clusters[i][j] = clusters2[i][j];
    }
  }

  for(int i = 0; i < nObs; i++){
    cluster_assignment[i] = cluster_assignment1[i];
  }
  for(int i = 0; i < newpart->nObs; i++){
    newpart->cluster_assignment[i] = cluster_assignment2[i];
  }
  for(int k = 0; k < K; k++){
    log_pi_ep(k, eta_CT);
  }
  for(int k = 0; k < newpart->K; k++){
    newpart->log_pi_ep(k, eta_CT);
  }
  get_pairwise();
  newpart->get_pairwise();
  indexing = std::vector<int>(indexing1);
  newpart->indexing = std::vector<int>(indexing2);
  return;
}

void Partition::format_partition(double* vector, int starting_index){
  for(int i = 0; i < nObs; i++){
    vector[starting_index + indexing[i]] = cluster_assignment[i] + 1; // Remember that R is 1-indexed
  }
  return;
}

void Partition::get_postmean(double* vector, int starting_index, bool LR){
  arma::mat Yloc;
  if(LR){
    Yloc = YLR;
  } else {
    Yloc = Y;
  }
  int cost_localID, cost_mainID;
  double ybar;
  for(int k =0; k < K; k++){
    arma::vec y_cl(cluster_config[k]);
    for(int t = 0; t < cluster_config[k]; t++){
      cost_localID = clusters[k][t];
      cost_mainID = indexing[cost_localID];
      y_cl[t] =  Yloc[ cost_mainID ];
    }
    ybar = arma::mean(y_cl);
    for(int t = 0; t < cluster_config[k]; t++){
      cost_localID = clusters[k][t];
      cost_mainID = indexing[cost_localID];
      vector[starting_index + cost_mainID] = ybar;
    }
  }
  return;
}

// Rcpp::List Partition::format_partition()
// {
//   Rcpp::List output_list;
//   Rcpp::NumericVector output_vec;
//   Rcpp::List tmp_list;

//   arma::vec tmp_vec = arma::zeros<vec>(nObs);   
//   for(int i = 0; i < nObs; i++){
//     tmp_vec[ indexing[i] ] = cluster_assignment[i] + 1; // Remember that R is 1-indexed
//   }
//   output_vec = Rcpp::wrap(tmp_vec);
//   output_vec.attr("dim") = R_NilValue;
//   output_list["cl_vec"] = output_vec;

//   tmp_list = Rcpp::List(K);
//   for(int k = 0; k < K; k++){
//     tmp_vec = arma::zeros<vec>(cluster_config[k]);
//     for(int i = 0; i < cluster_config[k]; i++){
//       tmp_vec(i) = indexing[ clusters[k][i] ] + 1; // Remember that R is 1-indexed
//     }
//     output_vec = Rcpp::wrap(arma::sort(tmp_vec, "ascend"));
//     output_vec.attr("dim") = R_NilValue;
//     tmp_list[k] = output_vec;
//   }
//   output_list["cl_list"] = tmp_list;
  
//   return output_list;  
// }

// bool Particle::Find_KM_Splits(bool A_or_B, int split_k, int num_splits, std::vector<std::vector<int> >& indices, std::vector<int>& ns){
//   // indices now contains numbers between 0 and nObs-1
  
//     // int SCiter = 2;
//     int n_cl = partition->cluster_config[split_k];
//     arma::mat U = arma::zeros<mat>(n_cl, 1); // holds the data passed to k-means
//     U.set_size(n_cl, 1);
//     for(int i = 0; i < n_cl; i++){
//       U(i,0)  = parameter[partition->clusters[split_k][i]];
//     }
//     arma::mat means;

//     indices.clear();
//     indices.resize(num_splits);
//     ns.clear();
//     ns.resize(num_splits);

//     double min_score;
//     int * membership;
//     bool status = false;
//     while(!status){
//       // status = arma::kmeans(means, U.t(), num_splits, random_subset, 10, false);
//       status = kmeans_repeat(means, U, num_splits, 10, min_score);
//       if(status == false){
//         cout << "clustering failed" << endl;
//         // go directly to next iteration: i.e. do kmeans again
//       }
//     }
//     // membership is long n_cl and has elements between 0 and num_splits-1
//     membership = which_is_nearest_k(means, U.t());
    
//     for(int j = 0; j<num_splits; j++){
//       ns[j] = 0;
//     }

//     // indices has the actual indices from 0 to nObs-1
//     for(int i = 0; i < n_cl; i++){
//       indices[membership[i]].push_back( partition->clusters[split_k][i] );
//       (ns[membership[i]])++;
//     }

//     int *components;
//     int *count;
//     for(int j = 0; j < num_splits; j++){
//       mat A_cluster = Submatrix(A_block, ns[j], ns[j], indices[j], indices[j]);
//       components = new int[ ns[j] ];
//       count = new int;
//       Connected_Components(A_cluster, ns[j],components,count);
//       if((*count)!=1){
//         // instead of iterating let's change indices:
//         // save the connected components of indices[j] as a vector
//         std::vector<std::vector<int> > new_indices;
//         new_indices.resize(*count);
//         for(int i = 0; i < ns[j]; i++){
//           new_indices[ components[i] ].push_back(indices[j][i]);
//         }
//         // we can eliminate the current element of indices[j] and ns[j]
//         indices[j].clear();
//         // now we update the previous element indices[j] and add new ones
//         for(int c = 0; c < (*count); c++){
//           if(c == 0){
//             indices[j] = new_indices[0];
//             ns[j] = new_indices[0].size();
//           } else {
//             indices.push_back(new_indices[c]);
//             ns.push_back(new_indices[c].size());
//           }
//         }
//       }
//       delete[] components;
//       delete count;
//     }
//     delete[] membership;

//     Partition_KSplit(A_or_B, split_k, indices, ns);
//     return true;
  
// }
