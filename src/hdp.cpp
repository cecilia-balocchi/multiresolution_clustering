#include <vector>  
#include <iostream>
#include <set>
#include <armadillo>
#include <algorithm>
#include <random>
#include <cmath>
#include "partition.h"
#include "hdp.h"
#include "various_functions.h"
#include "objects_functions.h"
using namespace std;
using namespace arma;
 
extern arma::mat Y;
extern arma::mat YLR;
extern arma::imat Mapping;
extern arma::imat InvMapping;
// extern double eta;
// extern double eta_LR;
// extern double eta_CT;
// extern double eta_TD;
extern double mu_H;
extern double mu_L;
extern double k_H;
extern double k_L;
extern double alpha_H;
extern double beta_H;
extern double alpha_L;
extern double beta_L;
extern double alpha1;
extern double alpha2;

// int nObs; // number of costumer over all the restaurantes
// int nRest; // number of restaurants
// std::vector<int> nObsRest; // vector containing the number of people in each restaurant
// std::vector<int> nTableRest; // vector containing the number of people in each restaurant
// std::vector<LPPartition> costumers_tables; // how the costumers are divided into tables, one for each restaurant
// LPPartition tables_dishes; // how the tables are divided into dishes
// LPPartition costumers_dishes; // overall partition
// int K; // this could be the size of the overall partition costumers_dishes


// Imagine the data to be clustered is represented in a vector, and we have a sort of index to understand 
// which y_i corresponds to a certain group.

// In general, we have a certain number of census tracts, but after we cluster them the real number of groups is
// reduced. 
// We could have N CTs, and one vector/pointer saying in which cluster they belong
// For each cluster, we know which CTs belong to that cluster, and the size of the cluster
// For each CT, we know the number of BGs and we have a vector containing the indices of BGs
// For each cluster of CTs, we can have the indices of all the BGs under the CTs in that cluster

// Y is a vector


Hdp::Hdp(){
  nObs = 0;
  nRest = 0;
  nObsRest = std::vector<int>();
  nTableRest = std::vector<int>();
  costumers_tables = std::vector<LPPartition>();
  tables_dishes = NULL;
  // K = 0;
  CTinRest = std::vector< std::vector<int> >();
  dishes_thetas = std::vector<double>();
  return;
}

Hdp::Hdp(LPHdp initial_hdp){
  nObs = initial_hdp->nObs;
  nRest = initial_hdp->nRest;
  nObsRest = std::vector<int>(nRest,0);
  for(int r = 0; r < nRest; r++){
    nObsRest[r] = initial_hdp->nObsRest[r];
  }
  nTableRest = std::vector<int>(nRest,0);
  for(int r = 0; r < nRest; r++){
    nTableRest[r] = initial_hdp->nTableRest[r];
  }
  // costumers_tables = std::vector<LPPartition>( initial_hdp->costumers_tables ); // this would copy, and the pointers would be pointing at the same objects
  costumers_tables = std::vector<LPPartition>( nRest );
  for(int r = 0; r < nRest; r++){
    costumers_tables[r] = new Partition( initial_hdp->costumers_tables[r] ); // new pointers that point at copies of the original objects
  }
  tables_dishes = new Partition(initial_hdp->tables_dishes);
  CTinRest = std::vector<std::vector<int> >(initial_hdp->CTinRest);
  dishes_thetas = std::vector<double>(initial_hdp->dishes_thetas);
  return;
}

Hdp::~Hdp(){
  for(int r = 0; r < nRest; r++){
    delete costumers_tables[r];
  }
  delete tables_dishes;
  // costumers_tables = NULL;
  tables_dishes = NULL;
  // nObsRest = NULL;
  nObs = 0;
  nRest = 0;
  // K = 0;
  return;
}

void Hdp::Copy_Hdp(LPHdp initial_hdp){
  for(int r = 0; r < nRest; r++){
    delete costumers_tables[r];
  }
  delete tables_dishes;
  tables_dishes = NULL;
  
  nObs = initial_hdp->nObs;
  nRest = initial_hdp->nRest;
  // K = initial_hdp->K;

  nObsRest.resize(nRest);
  for(int r = 0; r < nRest; r++){
    nObsRest[r] = initial_hdp->nObsRest[r];
  }
  nTableRest.resize(nRest);
  for(int r = 0; r < nRest; r++){
    nTableRest[r] = initial_hdp->nTableRest[r];
  }
  costumers_tables.resize(nRest);
  for(int r = 0; r < nRest; r++){
    // new pointers that point at copies of the original objects
    costumers_tables[r] = new Partition( initial_hdp->costumers_tables[r] ); 
  }
  tables_dishes = new Partition(initial_hdp->tables_dishes);
  CTinRest = std::vector<std::vector<int> >(initial_hdp->CTinRest);
  dishes_thetas = std::vector<double>(initial_hdp->dishes_thetas);
  return;
}

void Hdp::Merge_Rest(int r1, int r2){
  /*
  Merge_Rest needs to adjust the Hdp object: 
  - update the simple objects like nRest, nObsRest, nTableRest and CTinRest
  - update the partitions: costumers_tables and tables_dishes
  Combine_Restaurants creates the merged partition of costumers_tables[rmin], costumers_tables[rmax]
    without merging any tables.
  MoveTables_TableDish updates tables_dishes by renaming the tables: 
    the ones between rmin and rmax, and the rmax ones need to change their name in tables_dishes
  */
  int rmax = max(r1, r2);
  int rmin = min(r1, r2);
  // int rnew = rmin;
  int nRest_orig = nRest;
  nRest = nRest_orig-1;

  std::vector<int> nObsRest_orig(nObsRest);
  std::vector<int> nTableRest_orig(nTableRest);
  std::vector<LPPartition> costumers_tables_orig(costumers_tables); 

  nObsRest.resize(nRest);
  for(int i = rmax; i < nRest; i++){
    nObsRest[i] = nObsRest_orig[i+1];
  }
  nObsRest[rmin] += nObsRest_orig[rmax];

  // before changing nTableRest, we need these two for MoveTables_TableDish:
  // tables_rmax contains the mainID for the tables that from rmax are merged into rmin
  std::vector<int> tables_rmax(costumers_tables_orig[rmax]->K);
  for(int i = 0; i < costumers_tables_orig[rmax]->K; i++){
    tables_rmax[i] = Table_RestId_MainId(rmax, i, nTableRest);
  }
  // starting index will be the one right after the tables of rmin
  int starting_index = costumers_tables_orig[rmin]->K;
  for(int r = 0; r < rmin; r++){
    starting_index += costumers_tables_orig[r]->K;
  }

  nTableRest.resize(nRest);
  for(int i = rmax; i < nRest; i++){
    nTableRest[i] = nTableRest_orig[i+1];
  }
  nTableRest[rmin] += nTableRest_orig[rmax];

  // this makes a copy of the pointers, not of the objects
  nTableRest.resize(nRest);
  for(int i = rmax; i < nRest; i++){ // we are shifting the pointers
    costumers_tables[i] = costumers_tables_orig[i+1];
  }
  LPPartition newpartition = new Partition();
  newpartition->Combine_Restaurants(costumers_tables_orig[rmin], costumers_tables_orig[rmax]);
  delete costumers_tables_orig[rmax];
  delete costumers_tables_orig[rmin];
  costumers_tables[rmin] = newpartition;
  // now we need to fix tables_dishes
  tables_dishes->MoveTables_TableDish(tables_rmax, starting_index);
  // finally CTinRest (easy)
  std::vector<std::vector<int> > CTinRest_orig(CTinRest);
  for(int i = 0; i < CTinRest_orig[rmax].size();i++){
    CTinRest[rmin].push_back(CTinRest_orig[rmax][i]);
  }
  for(int i = rmax; i < nRest; i++){
    CTinRest[i].resize(0);
    CTinRest[i] = std::vector<int>(CTinRest_orig[i+1]);
  }
  // DO NOT delete newpartition
  return;
}

void Hdp::Split_Rest(int r1, int r2, std::vector<int> index1, std::vector<int> index2, double eta_CT, double eta_TD, double local_print){
  /*
  Split_Rest needs to adjust the Hdp object: 
  - update the simple objects like nRest, nObsRest, nTableRest and CTinRest
  - update the partitions: costumers_tables and tables_dishes
  cost_index1/2 contain the localID of the costumers which are in CT in index1/2
  tables1/2 contain the unique tables (localID) which costumers in cost_index1/2 occupy 
    (in order of appearance in cost_index1/2)
  SplitTables_TableDish updates tables_dishes: [uses tables1/2]
    removes the previous tables and adds the ones that were created when splitting
  Divide_Restaurants updates costumers_tables: [uses cost_index1/2] 
    specifically: costumers_tables[r1] and the new partition which gets added at the end of costumers_tables
  Remember: here the splits happens in a unique way
  */
  int nRest_orig = nRest;
  nRest = nRest_orig+1;
  if(r2 != nRest-1){
    cout << "Split_Rest: Problem!" << endl;
  }
  int new_r2 = nRest-1; // this should be equal to r2

  std::vector<int> nObsRest_orig(nObsRest);
  std::vector<int> nTableRest_orig(nTableRest);
  std::vector<LPPartition> costumers_tables_orig(costumers_tables); 

  // does costumer belong to a group in index1 or index2?
  std::vector<int> cost_index1, cost_index2;
  for(int i = 0; i < costumers_tables[r1]->nObs; i++){
    int r = InvMapping(costumers_tables[r1]->indexing[i],1);
    // if(local_print){
    //   printf("i: %d, ids[i]: %d, r: %d;",i, costumers_tables[r1]->indexing[i],r);
    // }
    bool r_in_index1 = false;
    for(int j = 0; j < index1.size(); j++){
      if(index1[j] == r){
        r_in_index1 = true;
        break;
      }
    }
    if(r_in_index1 == true){
      // if(local_print){
      //   cout << " is in sub-restaurant 1." << endl;
      // }
      cost_index1.push_back(i);
    } else {
      // if(local_print){
      //   cout << " is in sub-restaurant 2." << endl;
      // }
      cost_index2.push_back(i);
    }
  }
  nObsRest.resize(nRest);
  nObsRest[r1] = cost_index1.size();
  nObsRest[new_r2] = cost_index2.size();

  std::vector<int> tables1,tables2; // unique set of tables, but with right order
  int tab;
  bool newtab;
  for(int i = 0; i < cost_index1.size(); i++){
    tab = costumers_tables[r1]->cluster_assignment[cost_index1[i]];
    newtab = true;
    for(int c = 0; c < tables1.size(); c++){
      if(tables1[c] == tab){
        newtab = false;
      }
    }
    if(newtab){
      tables1.push_back(tab);
    }
  }
  for(int i = 0; i < cost_index2.size(); i++){
    tab = costumers_tables[r1]->cluster_assignment[cost_index2[i]];
    newtab = true;
    for(int c = 0; c < tables2.size(); c++){
      if(tables2[c] == tab){
        newtab = false;
      }
    }
    if(newtab){
      tables2.push_back(tab);
    }
  }
  // if(local_print){
  //   cout << "tables1 ";
  //   for(int c = 0; c < tables1.size(); c++) 
  //     cout << tables1[c] << " ";
  //   cout << endl;
  //   cout << "tables2 ";
  //   for(int c = 0; c < tables2.size(); c++) 
  //     cout << tables2[c] << " ";
  //   cout << endl;
  // }

  // update tables_dishes by removing the previous tables and adding the ones that were created when splitting
  tables_dishes->SplitTables_TableDish(r1, new_r2, tables1, tables2, nTableRest, eta_TD, local_print);
  nTableRest.resize(nRest);
  nTableRest[r1] = tables1.size();
  nTableRest[new_r2] = tables2.size();

  LPPartition newpart = new Partition();
  costumers_tables.resize(nRest);
  costumers_tables[r1]->Divide_Restaurants(newpart, cost_index1, cost_index2, eta_CT, local_print);
  costumers_tables[new_r2] = newpart;
  
  // update CTinRest: index1 contains the CT in rest r1 and same for index2
  CTinRest[r1].resize(0);
  CTinRest[r1] = std::vector<int>(index1);
  CTinRest.resize(nRest);
  CTinRest[new_r2] = std::vector<int>(index2);
  return;
}

void Hdp::Change_Rest(proposal_info* info, double eta_CT, double eta_TD, double local_print){
  /*
  - update the simple objects like nRest, nObsRest, nTableRest and CTinRest
  - update costumers_tables
  Note if info->kLR_new == nRest_orig we are just splitting.
  */
  // if(info->kLR_new == nRest_orig) // just split

  if(local_print){
    cout << info->jLR << " " << info->kLR_old << " " << info->kLR_new << endl << "cost_locID_j" << endl;
    // cout << nTableRest.size() << endl;
    // for(int r = 0; r < nTableRest.size(); r++){
    //   cout << nTableRest[r] << " ";
    // }
    // cout << endl;
    for(int i = 0; i < info->n_cost_j; i++){
      cout << info->cost_locID_j[i] << " ";
    }
    cout << endl << "cost_indexing_j" << endl;
    for(int i = 0; i < info->n_cost_j; i++){
      cout << info->cost_indexing_j[i] << " ";
    }
    cout << endl << "oldtables" << endl;
    for(int i = 0; i < info->n_cost_j; i++){
      cout << info->oldtables[i] << " ";
    }
    cout << endl << "new_tables" << endl;
    for(int i = 0; i < info->n_cost_j; i++){
      cout << info->new_tables[i] << " ";
    }
    cout << endl << "new_dishes" << endl;
    for(int i = 0; i < info->n_cost_j; i++){
      cout << info->new_dishes[i] << " ";
    }
    cout << endl;
  }

  bool add_rest = false;
  bool left_zero = false;

  int nRest_orig = nRest;
  if(info->kLR_new == nRest_orig){
    nRest += 1;  
    add_rest = true;
    // cout << "add_rest" << endl;
  } 
  if(info->cost_locID_NOTj.size() == 0){
    nRest -= 1;
    left_zero = true;
    // cout << "left_zero" << endl;
  }
  
  int r1 = info->kLR_old;
  int new_r2 = info->kLR_new;
  int min_r, max_r;
  if(r1 < new_r2){
    min_r = r1;
    max_r = new_r2;
  } else {
    min_r = new_r2;
    max_r = r1;
  }

  std::vector<int> nObsRest_orig(nObsRest);
  std::vector<int> nTableRest_orig(nTableRest); // TODO: not used?

  /* ~~~ update nObsRest ~~~ */ 

  if(add_rest && left_zero){
    // do nothing
  } 
  if(!add_rest && left_zero){
    nObsRest.resize(nRest);
    nObsRest[min_r] = nObsRest_orig[r1] + nObsRest_orig[new_r2];
    for(int r = max_r+1; r < nRest_orig; r++){
      nObsRest[r-1] = nObsRest_orig[r];
    }
  }
  if(!left_zero){
    if(add_rest){
      nObsRest.push_back(0);
    }
    nObsRest[r1] = info->cost_locID_NOTj.size(); 
    nObsRest[new_r2] += info->n_cost_j;
  }

  /* ~~~ update CTinRest ~~~ */ 

  std::vector<int> temp_vec;
  if(add_rest && left_zero){
    // does not change
  }
  if(!add_rest && left_zero){
    temp_vec = std::vector<int>(CTinRest[new_r2]);
    temp_vec.push_back(info->jLR_index);

    std::vector<std::vector<int> > orig_CTinRest(CTinRest);
    CTinRest.resize(nRest);

    CTinRest[min_r] = temp_vec;
    for(int r = max_r+1; r < nRest_orig;r++){
      CTinRest[r-1] = orig_CTinRest[r];
    }
  } 
  if(!left_zero){
    temp_vec = std::vector<int>(CTinRest[r1]);
    CTinRest[r1].resize(0);
    for(int i = 0; i < temp_vec.size();i++){
      if(temp_vec[i] != info->jLR_index){
        CTinRest[r1].push_back(temp_vec[i]);
      }
    }
    if(add_rest){
      CTinRest.push_back(std::vector<int>());
    }
    CTinRest[new_r2].push_back(info->jLR_index);
  }

  int n_tables_r2 = 0; // number of tables in new_r2 before changing it.
  if(!add_rest){
    // n_tables_r2 = tables_dishes->cluster_config[new_r2];
    n_tables_r2 = costumers_tables[new_r2]->K;
  }
  if(local_print) cout << "n_tables_r2 " << n_tables_r2 << endl;

  int temp;
  bool temp_bool;
  int orig_K, orig_nObs;
  if(add_rest && left_zero){
    // the costumers in the restaurant don't change, but tables and dishes do
    
    /* ~~~ let's change costumers_tables[r1] ~~~ */ 
    // costumers do not change
    // indexing does not change
    // K, cluster_config, clusters, cluster_assignment, log_prior change

    orig_K = costumers_tables[r1]->K;
    delete[] costumers_tables[r1]->cluster_config;
    for(int k = 0; k < orig_K; k++){
      delete[] costumers_tables[r1]->clusters[k];
    }
    delete[] costumers_tables[r1]->clusters;
    delete[] costumers_tables[r1]->log_prior;

    std::vector<int> seen_tables;
    std::vector< std::vector<int> > new_clusters;
    for(int i = 0; i < info->n_cost_j; i++){
      costumers_tables[r1]->cluster_assignment[i] = info->new_tables[i];

      temp_bool = false; // is info->new_tables[i] in seen_tables?
      for(int j = 0; j < seen_tables.size(); j++){
        if(info->new_tables[i] == seen_tables[j]){
          temp_bool = true;
          new_clusters[j].push_back( i );
          break;
        }
      }
      if(temp_bool == false){
        seen_tables.push_back(info->new_tables[i]);
        new_clusters.push_back(std::vector<int>(1, i));
      }
    }
    costumers_tables[r1]->K = new_clusters.size();

    // we could have a little check here:
    // if seen_tables[j] != j there could be some problems
    for(int j = 0; j < seen_tables.size(); j++){
      if(seen_tables[j] != j){
        cout << "Change_Rest problem, add_rest && left_zero, seen_tables[j] != j for j="<<j << endl;
      }
    }
    // let's also add a check that all the cluster_assignments are < K
    for(int i = 0; i < costumers_tables[r1]->nObs; i++){
      if(costumers_tables[r1]->cluster_assignment[i] >= costumers_tables[r1]->K){
        cout << "Change_Rest problem, add_rest && left_zero, cluster_assignment >= K for i="<<i << endl;
      }
    }

    costumers_tables[r1]->cluster_config = new int[costumers_tables[r1]->K];
    costumers_tables[r1]->clusters = new int*[costumers_tables[r1]->K];
    costumers_tables[r1]->log_prior = new double[costumers_tables[r1]->K];

    for(int k = 0; k < costumers_tables[r1]->K; k++){
      costumers_tables[r1]->cluster_config[k] = new_clusters[k].size();
      costumers_tables[r1]->clusters[k] = new int[costumers_tables[r1]->cluster_config[k]];
      for(int j = 0; j < costumers_tables[r1]->cluster_config[k]; j++){
        costumers_tables[r1]->clusters[k][j] = new_clusters[k][j];
      }
      costumers_tables[r1]->log_pi_ep(k, eta_CT);
    }
  }
  if(!add_rest && left_zero){
    // we merge the two restaurants in new_r2 and delete r1. 
    // The new partition is saved in min_r, max_r is removed and we shift everyone else.
    // In other words:
    // create a new partition by copying the pointer to new_r2, to which we'll add r1, which will represent the merge. 
    // this object will get save in min_r. The other partitions will get shifted in costumers_tables
    // the new costumers from r1 are added after the ones originally in new_r2
    // and so are the new tables created by the costumers in r1.

    // what happens if r1 == new_r2? Not possible because left zero means there was zero probability to choose r1.

    std::vector<LPPartition> orig_costumers_tables( nRest_orig );
    for(int r = 0; r < nRest_orig; r++){
      orig_costumers_tables[r] = costumers_tables[r]; // this should save only the pointers
    }

    LPPartition new_rest = orig_costumers_tables[new_r2];
    orig_nObs = new_rest->nObs;
    orig_K = new_rest->K;
    
    std::vector<int> new_cluster_assignment;
    std::vector<std::vector<int>> new_clusters;
    std::vector<int> seen_tables;

    new_rest->nObs += info->n_cost_j;
    
    for(int i = 0; i < orig_nObs; i++){
      new_cluster_assignment.push_back(new_rest->cluster_assignment[i]);
    }
    for(int k = 0; k < orig_K; k++){
      new_clusters.push_back( std::vector<int>() );
      for(int j = 0; j < new_rest->cluster_config[k]; j++){
        new_clusters[k].push_back(new_rest->clusters[k][j]);
      }
      seen_tables.push_back(k);
    }

    for(int i = 0; i < info->n_cost_j; i++){
      new_rest->indexing.push_back( info->cost_indexing_j[i] );
      new_cluster_assignment.push_back( info->new_tables[i] );
      
      temp_bool = false; // is info->new_tables[i] in seen_tables?
      for(int j = 0; j < seen_tables.size(); j++){
        if(info->new_tables[i] == seen_tables[j]){
          temp_bool = true;
          new_clusters[j].push_back( orig_nObs+i );
          break;
        }
      }
      if(temp_bool == false){
        seen_tables.push_back(info->new_tables[i]);
        temp = new_clusters.size();
        new_clusters.push_back(std::vector<int>(1, orig_nObs+i));
        if(new_clusters.size() == temp){
          cout << "Change_Rest problem, !add_rest && left_zero, new_clusters size did not increase" << endl;
        }
      }
    }
    new_rest->K = new_clusters.size();

    // we could have a little check here:
    // if seen_tables[j] != j there could be some problems
    for(int j = 0; j < seen_tables.size(); j++){
      if(seen_tables[j] != j){
        cout << "Change_Rest problem, !add_rest && left_zero, seen_tables[j] != j for j="<<j << endl;
      }
    }
    // let's also add a check that all the cluster_assignments are < K
    for(int i = 0; i < new_cluster_assignment.size(); i++){
      if(new_cluster_assignment[i] >= new_rest->K){
        cout << "Change_Rest problem, !add_rest && left_zero, cluster_assignment >= K r1,r2:" << r1 << " " << new_r2 << " for i="<<i << endl;
      }
    }
    if(new_clusters.size() != seen_tables.size()){
      cout << "Change_Rest problem, !add_rest && left_zero, new_clusters.size() != seen_tables.size()" << endl;
    }

    for(int k = 0; k < orig_K; k++){
      delete[] new_rest->clusters[k];
    }
    delete[] new_rest->cluster_config;
    delete[] new_rest->clusters;
    delete[] new_rest->cluster_assignment;
    delete[] new_rest->log_prior;
    for(int i = 0; i < orig_nObs; i++){
      delete[] new_rest->pairwise_assignment[i];
    }
    delete[] new_rest->pairwise_assignment;
    new_rest->cluster_assignment = new int[new_rest->nObs];
    new_rest->cluster_config = new int[new_rest->K];
    new_rest->clusters = new int*[new_rest->K];
    new_rest->log_prior = new double[new_rest->K];
    for(int k = 0; k < new_rest->K; k++){
      new_rest->cluster_config[k] = new_clusters[k].size();
      new_rest->clusters[k] = new int[new_rest->cluster_config[k]];
      for(int j = 0; j < new_rest->cluster_config[k]; j++){
        new_rest->clusters[k][j] = new_clusters[k][j];
      }
      new_rest->log_pi_ep(k, eta_CT);
    }
    for(int i = 0; i < new_rest->nObs; i++){
      new_rest->cluster_assignment[i] = new_cluster_assignment[i];
    }
    new_rest->get_pairwise();

    delete costumers_tables[r1];

    costumers_tables.resize(nRest);
    costumers_tables[min_r] = new_rest;
    for(int r = max_r+1; r < nRest_orig;r++){
      costumers_tables[r-1] = orig_costumers_tables[r];
    }
  }
  if(!left_zero){
    // the costumers in the restaurant change, together with tables and dishes

    /* ~~~ let's change costumers_tables[r1] ~~~ */
    
    orig_nObs = costumers_tables[r1]->nObs;
    orig_K = costumers_tables[r1]->K;
    int* orig_cluster_config = new int[orig_K];
    int* orig_cluster_assignment = new int[orig_nObs];
    double* orig_log_prior = new double[orig_K];
    for(int k = 0; k < orig_K; k++){
      orig_cluster_config[k] = costumers_tables[r1]->cluster_config[k];
      orig_log_prior[k] = costumers_tables[r1]->log_prior[k];
      delete[] costumers_tables[r1]->clusters[k];
    }
    for(int i = 0; i < orig_nObs; i++){
      orig_cluster_assignment[i] = costumers_tables[r1]->cluster_assignment[i];
    }
    std::vector<int> orig_indexing = std::vector<int>(costumers_tables[r1]->indexing);
    delete[] costumers_tables[r1]->cluster_config;
    delete[] costumers_tables[r1]->clusters;
    delete[] costumers_tables[r1]->log_prior;
    delete[] costumers_tables[r1]->cluster_assignment;

    std::vector<int> new_cluster_config;
    std::vector< std::vector<int> > new_clusters;
    std::vector<double> new_log_prior;
    std::vector<int> new_cluster_assignment;  // the order of the indices is the same as info->cost_locID_NOTj 
                                              // which is the original order minus cost_locID_j

    costumers_tables[r1]->nObs = info->cost_locID_NOTj.size();
    new_cluster_assignment.resize(costumers_tables[r1]->nObs);
    int k_newname = 0;
    for(int k = 0; k < orig_K; k++){
      if(info->remove_table[k] == 1){
        costumers_tables[r1]->K--;
      } else {          
        temp = orig_cluster_config[k];
        for(int i = 0; i < info->oldtables.size(); i++){
          if(info->oldtables[i] == k)
            temp--;
        }
        new_cluster_config.push_back(temp);
        if(temp == orig_cluster_config[k]){
          new_log_prior.push_back(orig_log_prior[k]);
        } else {
          new_log_prior.push_back( log(eta_CT) + lgamma(temp) );
        }

        temp_vec.resize(0);
        temp_vec.reserve(temp);
        for(int i = 0; i < info->cost_locID_NOTj.size(); i++){
          if(orig_cluster_assignment[ info->cost_locID_NOTj[i] ] == k){
            // note: info->cost_locID_NOTj[i] was the old ID, the new one is i
            temp_vec.push_back(i);
            new_cluster_assignment[i] = k_newname;
          }
        }
        if(local_print){
          if(temp_vec.size() != temp){
            cout << "temp and temp_vec.size not matching for k=" << k << " and k_newname=" << k_newname << endl;
          }
        }
        new_clusters.push_back(temp_vec);

        if(new_r2 == r1){
          // we change the name of the tables in info->new_tables
          for(int i = 0; i < info->n_cost_j; i ++){
            if(info->new_tables[i] == k){
              info->new_tables[i] = k_newname;
            }
          }
        }

        k_newname++;
      }
    }

    // let's add a check that all the cluster_assignments are < K
    for(int i = 0; i < costumers_tables[r1]->nObs; i++){
      if(new_cluster_assignment[i] >= costumers_tables[r1]->K){
        cout << "Change_Rest problem, !left_zero, cluster_assignment in r1 >= K r1,r2:" << r1 << " " << new_r2 << " for i="<<i << endl;
      }
    }
    if(k_newname != costumers_tables[r1]->K){
      cout << "Change_Rest problem, !left_zero, k_newname != K: " << k_newname << " " << costumers_tables[r1]->K << endl;
    }

    costumers_tables[r1]->cluster_config = new int[costumers_tables[r1]->K];
    costumers_tables[r1]->clusters = new int*[costumers_tables[r1]->K];
    costumers_tables[r1]->log_prior = new double[costumers_tables[r1]->K];
    costumers_tables[r1]->cluster_assignment = new int[costumers_tables[r1]->nObs];

    for(int k = 0; k < costumers_tables[r1]->K; k++){
      costumers_tables[r1]->cluster_config[k] = new_cluster_config[k];
      costumers_tables[r1]->clusters[k] = new int[ costumers_tables[r1]->cluster_config[k] ];
      for(int j = 0; j < new_clusters[k].size(); j++){
        costumers_tables[r1]->clusters[k][j] = new_clusters[k][j];
      }
      costumers_tables[r1]->log_prior[k] = new_log_prior[k];
    }
    for(int i = 0; i < costumers_tables[r1]->nObs; i++){
      costumers_tables[r1]->cluster_assignment[i] = new_cluster_assignment[i];
    }
    costumers_tables[r1]->indexing.resize(costumers_tables[r1]->nObs);
    for(int i = 0; i < costumers_tables[r1]->nObs; i++){
      costumers_tables[r1]->indexing[i] = orig_indexing[ info->cost_locID_NOTj[i] ];
    }

    for(int i = 0; i < orig_nObs; i++){
      delete[] costumers_tables[r1]->pairwise_assignment[i]; 
    }
    delete[] costumers_tables[r1]->pairwise_assignment; 
    costumers_tables[r1]->get_pairwise();

    delete[] orig_cluster_config;
    delete[] orig_log_prior;
    delete[] orig_cluster_assignment;
  
    /* ~~~ let's create costumers_tables[new_r2] ~~~ */
    // there are some differences depending if add_rest or not

    LPPartition new_rest;
    if(add_rest){
      new_rest = new Partition();
    } else {
      new_rest = costumers_tables[new_r2];
    }
    
    orig_nObs = new_rest->nObs;
    orig_K = new_rest->K;
    new_cluster_assignment.resize(0);
    new_cluster_config.resize(0);
    new_clusters.resize(0);
    std::vector<int> seen_tables;

    new_rest->nObs += info->n_cost_j;
    if(!add_rest){
      for(int i = 0; i < orig_nObs; i++){
        new_cluster_assignment.push_back(new_rest->cluster_assignment[i]);
      }
      for(int k = 0; k < orig_K; k++){
        new_cluster_config.push_back(new_rest->cluster_config[k]); // **********
        new_clusters.push_back( std::vector<int>() );
        for(int j = 0; j < new_rest->cluster_config[k]; j++){
          new_clusters[k].push_back(new_rest->clusters[k][j]);
        }
        seen_tables.push_back(k);
      }
    }
    
    for(int i = 0; i < info->n_cost_j; i++){
      new_rest->indexing.push_back( orig_indexing[info->cost_locID_j[i]] ); // NOTE: orig_indexing is from r1. why not using info->cost_indexing_j[i]??
      
      temp_bool = false; // is info->new_tables[i] in seen_tables?
      for(int j = 0; j < seen_tables.size(); j++){
        if(info->new_tables[i] == seen_tables[j]){
          temp_bool = true;
          new_cluster_config[j]++;
          new_clusters[j].push_back( orig_nObs+i );
          new_cluster_assignment.push_back( j );
          break;
        }
      }
      if(temp_bool == false){
        seen_tables.push_back(info->new_tables[i]);
        new_cluster_config.push_back(1);
        new_clusters.push_back(std::vector<int>(1, orig_nObs+i));
        new_cluster_assignment.push_back( new_clusters.size()-1 );
      }

    }
    new_rest->K = new_cluster_config.size();

    // we could have a little check here:
    // if seen_tables[j] != j there could be some problems
    // temp_bool = false;
    // for(int j = 0; j < seen_tables.size(); j++){
    //   if(seen_tables[j] != j){
    //     cout << "Change_Rest problem, !left_zero, in new_r2 seen_tables[j] != j for j="<<j << endl;
    //     temp_bool = true;
    //   }
    // }
    // if(temp_bool){
    //   cout << "seen_tables: ";
    //   for(int j = 0; j < seen_tables.size(); j++){
    //     cout << seen_tables[j] << " ";
    //   }
    //   cout << endl;
    // }
    // let's also add a check that all the cluster_assignments are < K
    temp_bool = false;
    for(int i = 0; i < new_rest->nObs; i++){
      if(new_cluster_assignment[i] >= new_rest->K){
        cout << "Change_Rest problem, !left_zero, cluster_assignment in new_r2 >= Kr1,r2:" << r1 << " " << new_r2 << " for i="<<i << endl;
        temp_bool = true;
      }
    }
    if(temp_bool){
      cout << "old cluster_assignment: ";
      for(int i = 0; i < orig_nObs; i++){
        cout << new_rest->cluster_assignment[i] << " ";
      }
      cout << endl << "new_tables: ";
      for(int i = 0; i < info->n_cost_j; i++){
        cout << info->new_tables[i] << " ";
      }
      cout << endl;
    }
    

    if(!add_rest){
      for(int k = 0; k < orig_K; k++){
        delete[] new_rest->clusters[k];
      }
      delete[] new_rest->cluster_config;
      delete[] new_rest->clusters;
      delete[] new_rest->cluster_assignment;
      delete[] new_rest->log_prior;
      for(int i = 0; i < orig_nObs; i++){
        delete[] new_rest->pairwise_assignment[i]; 
      }
      delete[] new_rest->pairwise_assignment; 
    }
    new_rest->cluster_assignment = new int[new_rest->nObs];
    new_rest->cluster_config = new int[new_rest->K];
    new_rest->clusters = new int*[new_rest->K];
    new_rest->log_prior = new double[new_rest->K];
    for(int k = 0; k < new_rest->K; k++){
      new_rest->cluster_config[k] = new_cluster_config[k];
      new_rest->clusters[k] = new int[new_cluster_config[k]];
      for(int j = 0; j < new_cluster_config[k]; j++){
        new_rest->clusters[k][j] = new_clusters[k][j];
      }
      new_rest->log_pi_ep(k, eta_CT);
    }
    for(int i = 0; i < new_rest->nObs; i++){
      new_rest->cluster_assignment[i] = new_cluster_assignment[i];
    }
    new_rest->get_pairwise();

    if(add_rest){
      costumers_tables.push_back(new_rest);
    }    
  }


  /* ~~~ let's change tables_dishes ~~~ */ 
  
  orig_nObs = tables_dishes->nObs;
  orig_K = tables_dishes->K;

  std::vector<int> newtables,newdishes; 
  for(int i = 0; i < info->n_cost_j; i++){
    if(info->new_tables[i] >= n_tables_r2){ // info->new_tables[i] was already in r2
      temp_bool = false; // is info->new_tables[i] in newtables?
      for(int j = 0; j < newtables.size(); j++){
        if(newtables[j] == info->new_tables[i]){
          temp_bool = true;
          break;
        }
      }
      if(temp_bool == false){
        newtables.push_back(info->new_tables[i]);
      }
    } 
    if(info->new_dishes[i] >= orig_K){
      temp_bool = false; // is info->new_tables[i] in newdishes?
      for(int j = 0; j < newdishes.size(); j++){
        if(newdishes[j] == info->new_dishes[i]){
          temp_bool = true;
          break;
        }
      }
      if(temp_bool == false){
        newdishes.push_back(info->new_dishes[i]);
      }
    }
  }
  int less_tables = 0;
  int less_dishes = 0;
  std::vector<int> new_dishes_names(orig_K,-1);
  if(local_print) cout << "Starting new_dishes_names" << endl;
  int shift = 0;
  for(int k = 0; k < orig_K; k++){
    if(info->remove_dish[k] == 1){
      shift++;
      less_dishes++;
    } else {
      new_dishes_names[k] = k-shift;
    }
    if(local_print) cout << k << " " << new_dishes_names[k] << endl;
  }
  std::vector<int> new_tables_names(orig_nObs,-1);
  if(local_print) cout << "Starting new_tables_names" << endl;
  int starting_index = 0;
  int shift_newtables = 0;
  shift = 0;
  for(int r = 0; r < nRest_orig; r++){
    if(local_print) cout << "Restaurant " << r << endl;
    if(r == r1){
      if(local_print) cout << "Restaurant " << r << " is r1." << endl;
      for(int i = 0; i < nTableRest[r]; i++){
        if(info->remove_table[i] == 1){
          shift++;
          less_tables++;
        } else {
          new_tables_names[starting_index+i] = starting_index + i - shift;
        }
        if(local_print) cout << i << " " << starting_index+i << " " << new_tables_names[starting_index+i] << endl;
      }

      if((!add_rest && left_zero) && (new_r2 == max_r)){
        if(local_print) cout << "flag1" << endl;
        shift -= nTableRest[new_r2];
        shift -= newtables.size();
        shift_newtables = starting_index + nTableRest[new_r2]; 
      }
    } else if(r == new_r2){
      if(local_print) cout << "Restaurant " << r << " is r2." << endl;
      if((!add_rest && left_zero) && (new_r2 == max_r)){
        if(local_print) cout << "flag2" << endl;
        for(int i = 0; i < nTableRest[r]; i++){
          new_tables_names[starting_index + i] = shift_newtables - nTableRest[new_r2] + i;
          if(local_print) cout << i << " " << starting_index+i << " " << new_tables_names[starting_index+i] << endl;
        }
        shift += nTableRest[new_r2];
      } else { // same as the other r
        if(local_print) cout << "flag3" << endl;
        for(int i = 0; i < nTableRest[r]; i++){
          new_tables_names[starting_index + i] = starting_index + i - shift;
          if(local_print) cout << i << " " << starting_index+i << " " << new_tables_names[starting_index+i] << endl;
        }
      }
    } else {
      for(int i = 0; i < nTableRest[r]; i++){
        new_tables_names[starting_index + i] = starting_index + i - shift;
        if(local_print) cout << i << " " << starting_index+i << " " << new_tables_names[starting_index+i] << endl;
      }
    } 
    starting_index += nTableRest[r];
    if(r == new_r2){
      if( !(!add_rest && left_zero) || !(new_r2 == max_r) ){ // = !((!add_rest && left_zero) && (new_r2 == max_r))
        if(local_print) cout << "flag4" << endl;
        shift_newtables = starting_index - shift;
        shift -= newtables.size();
      }
    }
  }
  if(add_rest){ // before we never reached new_r2 and now we need to initialize
    shift_newtables = starting_index - shift;
  }
  // for(int r = 0; r < nRest_orig; r++){
  //   if(r == r1){
  //     for(int i = 0; i < nTableRest[r]; i++){
  //       if(info->remove_table[i] == 1){
  //         shift++;
  //         less_tables++;
  //       } else {
  //         new_tables_names[starting_index+i] = starting_index + i - shift;
  //       }
  //     }
  //   } else {
  //     for(int i = 0; i < nTableRest[r]; i++){
  //       new_tables_names[starting_index + i] = starting_index + i - shift;
  //     }
  //   } 
  //   starting_index += nTableRest[r];
  //   if(r == new_r2){
  //     shift -= newtables.size();
  //     shift_newtables = starting_index;
  //   }
  // }
  
  tables_dishes->nObs = orig_nObs - less_tables + newtables.size();
  tables_dishes->K = orig_K - less_dishes + newdishes.size();

  if(local_print){
    cout << "tables_dishes->nObs " << tables_dishes->nObs << " " << orig_nObs << " " << less_tables << " " << newtables.size() << endl;
    cout << "tables_dishes->K " << tables_dishes->K << " " << orig_K << " " << less_dishes << " " << newdishes.size() << endl;
  }
  // i nuovi dishes (nuovi clusters) vengono aggiunti alla fine.
  
  std::vector<int> cluster_assignment_v;
  std::vector<std::vector<int> > clusters_v;
  cluster_assignment_v.resize(tables_dishes->nObs);
  clusters_v.resize(tables_dishes->K);
  
  int cl;
  for(int i=0; i < orig_nObs; i++){
    cl = tables_dishes->cluster_assignment[i];
    if(new_tables_names[i] >= 0){ // we don't want removed tables
      clusters_v[ new_dishes_names[cl] ].push_back( new_tables_names[i] ); ///////
      cluster_assignment_v[ new_tables_names[i] ] = new_dishes_names[cl];
    }
  }
  for(int i=0; i < newtables.size(); i++){
    for(int j = 0; j < info->n_cost_j; j++){
      if(info->new_tables[j] == newtables[i]){
        cl = info->new_dishes[j];
        break;
      }
    }
    if(cl < orig_K){
      cl = new_dishes_names[cl];
    } else {
      for(int j = 0; j < newdishes.size(); j++){
        if(newdishes[j] == cl){
          cl = orig_K - less_dishes + j;
          break;
        }
      }
    }
    clusters_v[ cl ].push_back( shift_newtables+i );
    cluster_assignment_v[ shift_newtables+i ] = cl;
  }

  delete[] tables_dishes->cluster_config;
  for(int k = 0; k < orig_K; k++){
    delete[] tables_dishes->clusters[k];
  }
  delete[] tables_dishes->clusters;
  delete[] tables_dishes->cluster_assignment;
  delete[] tables_dishes->log_prior;
  for(int i = 0; i < orig_nObs; i++){
    delete[] tables_dishes->pairwise_assignment[i]; 
  }
  delete[] tables_dishes->pairwise_assignment; 

  tables_dishes->cluster_config = new int[tables_dishes->K];
  tables_dishes->clusters = new int*[tables_dishes->K];
  tables_dishes->cluster_assignment = new int[tables_dishes->nObs];
  tables_dishes->log_prior = new double[tables_dishes->K];
  for(int i = 0; i < tables_dishes->K; i++){
    tables_dishes->cluster_config[i] = clusters_v[i].size();
    tables_dishes->clusters[i] = new int[tables_dishes->cluster_config[i] ];
    for(int j = 0; j < tables_dishes->cluster_config[i]; j++){
      tables_dishes->clusters[i][j] = clusters_v[i][j];
    }
  }
  for(int i = 0; i < tables_dishes->nObs; i++){
    tables_dishes->cluster_assignment[i] = cluster_assignment_v[i];
  }
  for(int k = 0; k < tables_dishes->K; k++){
    tables_dishes->log_pi_ep(k, eta_TD);
  }
  tables_dishes->get_pairwise();
  tables_dishes->indexing = std::vector<int>(tables_dishes->nObs);
  for(int i = 0; i < tables_dishes->nObs; i++){
    tables_dishes->indexing[i] = i;
  }

  /* ~~~ now nTableRest ~~~ */

  if(add_rest && left_zero){
    nTableRest[r1] -= less_tables;
    nTableRest[r1] += newtables.size();
  }
  if(!add_rest && left_zero){
    nTableRest[min_r] = nTableRest[new_r2] + newtables.size();
    std::vector<int> orig_nTableRest(nTableRest);
    nTableRest.resize(nRest);
    for(int r = max_r+1; r < nRest_orig;r++){
      nTableRest[r-1] = orig_nTableRest[r];
    }
  }
  if(!left_zero){
    nTableRest[r1] -= less_tables;
    if(add_rest){
      nTableRest.push_back(0);
    }
    nTableRest[new_r2] += newtables.size();
  }
  
  // update dishes_thetas
  dishes_thetas = get_dishes_thetas();

  return;
}

void Hdp::Initialize_Hdp(int n, vector<vector<int> > ctinrest, vector<vector<int> > rests, double eta_CT, double eta_TD){
  // ctinrest gives the Census tracts in each restaurant. Potentially from it we could find rests
  // rests is a vector of vector of the same length as the number of restaurants
  // each element of rests contains the mainID of the costumers for each restaurant
  // Initialize_Hdp puts every person in each restaurant in one cluster and each restaurant in one cluster
  nObs = n;
  nRest = rests.size();
  nObsRest = vector<int>(nRest,0);
  nTableRest = vector<int>(nRest,0);
  costumers_tables = vector<LPPartition>( nRest );
  for(int r = 0; r < nRest; r++){
    nObsRest[r] = rests[r].size();
    nTableRest[r] = 1;
    costumers_tables[r] = new Partition;
    costumers_tables[r]->Initialize_Partition( nObsRest[r], rests[r], eta_CT ); // one table in each restaurant
  }
  tables_dishes = new Partition; 
  tables_dishes->Initialize_Partition( nRest, eta_TD ); // all the tables serve the same dish
  CTinRest = std::vector<std::vector<int> >(ctinrest);
  dishes_thetas.resize(1);
  dishes_thetas[0] = arma::mean(arma::mean(Y));
  return;
}

void Hdp::Initialize_Hdp_nRestDish(int n, vector<vector<int> > ctinrest, vector<vector<int> > rests, double eta_CT, double eta_TD){
  // nj is a vector of the same length as the number of restaurants
  // each element of nj is the number of costumers for each restaurant
  // Initialize_Hdp puts every person in each restaurant in one cluster and each restaurant in one cluster
  nObs = n;
  nRest = rests.size();
  nObsRest = vector<int>(nRest,0);
  nTableRest = vector<int>(nRest,0);
  costumers_tables = vector<LPPartition>( nRest );
  dishes_thetas.resize(nRest);
  double temp;
  for(int r = 0; r < nRest; r++){
    nObsRest[r] = rests[r].size();
    nTableRest[r] = 1;
    costumers_tables[r] = new Partition;
    costumers_tables[r]->Initialize_Partition( nObsRest[r], rests[r], eta_CT ); // one table in each restaurant
    temp = 0.0;
    for(int i = 0; i < rests[r].size(); i++){
      temp += Y(rests[r][i]);
    }
    dishes_thetas[r] = temp/rests[r].size();
  }
  tables_dishes = new Partition; 
  tables_dishes->Initialize_Partition_nclusters( nRest, eta_TD ); // all the tables serve different dishes
  CTinRest = std::vector<std::vector<int> >(ctinrest);
  return;
}

void Hdp::Initialize_Hdp_nclusters(int n, vector<vector<int> > ctinrest, vector<vector<int> >  rests, double eta_CT, double eta_TD){
  nObs = n;
  nRest = rests.size();
  nObsRest = vector<int>(nRest,0);
  nTableRest = vector<int>(nRest,0);
  costumers_tables = vector<LPPartition>( nRest );
  int nDish = 0;
  for(int r = 0; r < nRest; r++){
    nObsRest[r] = rests[r].size();
    nTableRest[r] = nObsRest[r];
    nDish += nTableRest[r]; // each restaurant serves different dishes
    costumers_tables[r] = new Partition;
    costumers_tables[r]->Initialize_Partition_nclusters( nObsRest[r], rests[r], eta_CT ); // nj[r] tables in each restaurant
    for(int i = 0; i < rests[r].size(); i++){
      dishes_thetas.push_back(Y(rests[r][i]));
    }
  }
  tables_dishes = new Partition;
  tables_dishes->Initialize_Partition_nclusters( nDish, eta_TD ); // all the tables serve different dishes
  CTinRest = std::vector<std::vector<int> >(ctinrest);
  return;
}

void Hdp::Initialize_Hdp_part(int n, vector<vector<int> > ctinrest, vector<vector<int> >  rests, int part_init_ptr[], int* K_init_ptr, double eta_CT, double eta_TD){
  // this only works if rests has length 1
  nObs = n;
  nRest = rests.size();
  if(nRest != 1){
    cout << "Error in Initialize_Hdp_part: nRest should have length 1." << endl;
  }
  nObsRest = vector<int>(nRest,0);
  nTableRest = vector<int>(nRest,0);
  costumers_tables = vector<LPPartition>( nRest );
  int nDish = *K_init_ptr; // each table serves different dishes
  dishes_thetas.resize(nDish);
  // no for loop (only r = 0)
  nObsRest[0] = rests[0].size();
  nTableRest[0] = *K_init_ptr;
  costumers_tables[0] = new Partition;
  costumers_tables[0]->Initialize_Partition( nObsRest[0], rests[0], eta_CT, part_init_ptr, K_init_ptr ); // nj[0] tables in each restaurant
  // no for loop  
  tables_dishes = new Partition;
  tables_dishes->Initialize_Partition_nclusters( nDish, eta_TD ); // all the tables serve different dishes
  CTinRest = std::vector<std::vector<int> >(ctinrest);
  double temp;
  int ntemp;
  for(int k = 0; k < nDish; k++){
    temp = 0.0;
    ntemp = 0;
    for(int i = 0; i < rests[0].size(); i++){
      if(part_init_ptr[rests[0][i]] == (k+1)){
        temp += Y(rests[0][i]);
        ntemp++;
      }
    }
    dishes_thetas[k] = temp/ntemp;
  }
  return;
}

void Hdp::Initialize_Hdp_part_nrest(int n, vector<vector<int> > ctinrest, vector<vector<int> >  rests, int part_init_ptr[], double eta_CT, double eta_TD){
  // I think this can work even if rests is not = Rests (the one containing the BG in each CT)
  nObs = n;
  nRest = rests.size();
  nObsRest = vector<int>(nRest,0);
  nTableRest = vector<int>(nRest,0);
  costumers_tables = vector<LPPartition>( nRest );
  int nDish = 0;
  std::vector<int> dishes;           // it's like a gamma_unique (with the unique values of the cluster labels)
  std::vector<int> tables_in_dishes; // it's like a gamma (a vector with the cluster assignments)
  for(int r = 0; r < nRest; r++){
    nObsRest[r] = rests[r].size();
    costumers_tables[r] = new Partition;
    std::vector<int> local_dishes;
    local_dishes = costumers_tables[r]->Initialize_Partition( nObsRest[r], rests[r], eta_CT, part_init_ptr); // nj[r] tables in each restaurant
    if(nDish == 0){
      nDish += local_dishes.size();
      for(int k = 0; k < local_dishes.size(); k++){
        dishes.push_back(local_dishes[k]);
        tables_in_dishes.push_back(k);
      }
    } else {
      for(int k = 0; k < local_dishes.size(); k++){
        bool new_dish = true;
        for(int h = 0; h < dishes.size(); h++){
          if(local_dishes[k] == dishes[h]){
            new_dish = false;
            tables_in_dishes.push_back(h);
            break;
          }
        }
        if(new_dish == true){
          dishes.push_back(local_dishes[k]);
          tables_in_dishes.push_back(nDish);
          nDish++;
        }
      }
    }
    nTableRest[r] = costumers_tables[r]->K;
    // nDish += nTableRest[r]; // each restaurant serves different dishes
  }
  tables_dishes = new Partition;
  // for(int i = 0; i < tables_in_dishes.size(); i++){
  //   cout << tables_in_dishes[i] << " ";
  // }
  // cout << endl;
  // for(int i = 0; i < dishes.size(); i++){
  //   cout << dishes[i] << " ";
  // }
  // cout << endl;
  tables_dishes->Initialize_Partition( (int)tables_in_dishes.size(), eta_TD, tables_in_dishes, dishes ); // all the tables serve different dishes
  CTinRest = std::vector<std::vector<int> >(ctinrest);
  dishes_thetas = get_dishes_thetas();
  return;
}

void Hdp::Print_Hdp_Short(){
  cout << "Total number of costumers: " << nObs << endl;
  cout << "Number of restaurants: " << nRest << endl;
  cout << "Costumers in each restaurant: ";
  for(int r = 0; r < nRest; r++){
    cout << nObsRest[r] << " ";
  }
  cout << endl;
  cout << "Tables in each restaurant: ";
  for(int r = 0; r < nRest; r++){
    cout << nTableRest[r] << " ";
  }
  cout << endl;
  cout << "Total number of tables: " << tables_dishes->nObs << endl;
  cout << "Total number of dishes: " << tables_dishes->K << endl;
  return;
}

void Hdp::Print_Hdp(){
  cout << "Total number of costumers: " << nObs << endl;
  cout << "Number of restaurants: " << nRest << endl;
  cout << "Costumers in each restaurant: ";
  for(int r = 0; r < nRest; r++){
    cout << nObsRest[r] << " ";
  }
  cout << endl;
  cout << "Tables in each restaurant: ";
  for(int r = 0; r < nRest; r++){
    cout << nTableRest[r] << " ";
  }
  cout << endl;
  cout << "Costumers in tables in each restaurant: ";
  for(int r = 0; r < nRest; r++){
    for(int k = 0; k < costumers_tables[r]->K; k++){
      cout << costumers_tables[r]->cluster_config[k] << " ";
    }
    cout << "-- ";
  }
  cout << endl;
  // cout << "Total number of dishes: " << K << endl;
  cout << "- Partition of tables into dishes: " << endl;
  tables_dishes->Print_Partition();
  cout << "Theta values for each dish: (" << dishes_thetas.size() << "): " << endl;
  for(int i = 0; i < dishes_thetas.size(); i++)
    cout << dishes_thetas[i] << "   ";
  cout << endl;
  return;
}

void Hdp::Print_Y(){
  LPPartition part = tables_dishes;
  int table_mainID, trest, table_restID, cost_restID, cost_mainID;
  for(int k =0; k < part->K; k++){
    for(int t = 0; t < part->cluster_config[k]; t++){
      table_mainID = part->clusters[k][t];
      Table_MainId_RestId(table_mainID, nTableRest, trest, table_restID);
      for(int c = 0; c < costumers_tables[trest]->cluster_config[table_restID]; c++){
        cost_restID = costumers_tables[trest]->clusters[table_restID][c];
        // cost_mainID = Costumer_RestId_MainId(trest, cost_restID, nObsRest);
        cost_mainID = costumers_tables[trest]->indexing[cost_restID];
        cout << Y( cost_mainID ) << " ";
      }
    }
    cout << endl;
  }
}

std::vector<bool> Hdp::SampleSM(int ngibbs, double temp, std::default_random_engine& eng, double eta_CT, double eta_TD, bool local_print){
  // not sure if I need to change costumers_dishes here. NO
  // in fact CostTable does not care if there is an additional dish
  
  bool tmp;
  // std::vector<bool> changed(nRest+1, false); // this was to record all the acceptance, for all rest
  std::vector<bool> changed(2, false); 
  for(int r = 0; r < nRest; r++){
    if(costumers_tables[r]->nObs > 1){
      tmp = SampleSM_CostTable(r, ngibbs, temp, eng, eta_CT, eta_TD, local_print);
      if(tmp){
        changed[0] = true;
      }
    }
  }
  if(tables_dishes->nObs > 1){
    changed[1] = SampleSM_TableDishes(ngibbs, temp, eng, eta_TD); 
  }
  return changed;
} 

bool Hdp::SampleSM_TableDishes(int ngibbs, double temp, std::default_random_engine& eng, double eta_TD){
  // do we really need to copy this into startingP? can we just pass this? is startingP modified into ProposalSM
  // LPHdp startingP = new Hdp(this);

  LPHdp proposal = new Hdp(this);
  double lqratio, lpratio, llikratio;
  // proposal->ProposalSM(startingP, ngibbs, lqratio, lpratio, llikratio);
  proposal->ProposalSM_TableDishes(this, ngibbs, lqratio, lpratio, llikratio, eng, eta_TD);
  double a = lqratio + temp * (lpratio + llikratio);
  a = exp(a);
  a = min(1.0,a);
  
  // cout << lqratio << " " << lpratio << " " << llikratio << endl;
  // cout << exp(lqratio) << " " << exp(lpratio) << " " << exp(llikratio) << " ";
  // cout << a << endl;
  
  // default_random_engine generator(time(0));
  uniform_real_distribution<double> distribution(0.0,1.0);

  bool returned_value;
  if(distribution(eng) < a){ // distribution(generator) < a
    // change the current partition into proposal
    // cout << "**Tables_dishes changed." << endl;
    // cout << "**BEFORE" << endl;
    // tables_dishes->Print_Partition();
    Copy_Hdp(proposal);
    // cout << "**AFTER" << endl;
    // tables_dishes->Print_Partition();
    returned_value = true;
  } else {
    returned_value = false;
  }
  // delete startingP;
  delete proposal;
  return returned_value;
}

void Hdp::ProposalSM_TableDishes(LPHdp startingH, int ngibbs, double& lqratio, double& lpratio, double& llikratio, std::default_random_engine& myRandomEngine, double eta_TD){
  // Should modify proposal, not startingP, update the values of lqratio, lpratio, llikratio
  // It's updading proposal's tables_dishes, BUT NOT costumers_dishes
  
  // std::random_device myRandomDevice;
  // unsigned seed = myRandomDevice();
  // std::default_random_engine myRandomEngine(seed);

  LPPartition startingP = startingH->tables_dishes;
  int i1, i2, c1, c2;
  std::vector<int> i12 = sample_without_replacement(2, startingP->nObs, myRandomEngine);
  i1 = i12[0];
  i2 = i12[1];
  c1 = startingP->cluster_assignment[i1];
  c2 = startingP->cluster_assignment[i2];
  if(c2 < c1){
    int tmp = i1;
    i1 = i2;
    i2 = tmp;
    tmp = c1;
    c1 = c2;
    c2 = tmp;
  }

  // cout << "Sampled tables: " << i1 << " " << i2 << endl;
  // cout << "Corresponding clusters: " << c1 << " " << c2 << endl;
  // samecl_index contains the other people in those cluster (but not i1, i2) TablesMainID
  std::vector<int> samecl_index = startingP->get_samecl_index(i1, i2, c1, c2);

  // launch will NOT modify the current partition, but update new_c1_ind, new_c2_ind,samecl_newcl
  std::list<int> new_c1_ind,new_c2_ind; // lists are optimized to remove elements in any position (use "remove")
  std::vector<int> samecl_newcl;
  LaunchSM_TableDishes(startingH, i1, i2, samecl_index, c1, c2, new_c1_ind, new_c2_ind, samecl_newcl, ngibbs, myRandomEngine);
  
  if(c1 == c2){
    llikratio = alpha1 * (- startingH->dish_llikelihood(c1)); // startingH is not changed
    lqratio = 0;
    int new_c1 = c1;
    int new_c2 = startingP->K;
    ProposalSM_TableDishes_iter(true, startingH, samecl_index, new_c1, new_c2, new_c1_ind, new_c2_ind,samecl_newcl, lqratio, myRandomEngine);
    int nc1 = new_c1_ind.size();
    int nc2 = new_c2_ind.size();
    lpratio = log(eta_TD) + lbeta(nc1,nc2);
    llikratio += alpha1 * (startingH->dish_llikelihood(new_c1_ind) + startingH->dish_llikelihood(new_c2_ind));

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

    tables_dishes->Split(c1, index1, index2, nc1, nc2, eta_TD);
    delete[] index1;
    delete[] index2;
    // K = K + 1;
  } else {
    // maybe at the end: merge everyone into one cluster
    lpratio = - log(eta_TD) - lbeta(startingP->cluster_config[c1], startingP->cluster_config[c2]);
    lqratio = 0;
    
    llikratio = alpha1 * (- startingH->dish_llikelihood(c1) - startingH->dish_llikelihood(c2));

    ProposalSM_TableDishes_iter(false, startingH, samecl_index, c1, c2, new_c1_ind, new_c2_ind,samecl_newcl, lqratio, myRandomEngine);
    tables_dishes->Merge(c1,c2, eta_TD);
    // K = K - 1;
    
    llikratio += alpha1 * (dish_llikelihood(c1));
  }
  dishes_thetas = get_dishes_thetas();
  return;
}

void Hdp::ProposalSM_TableDishes_iter(bool SplitMerge, LPHdp startingH, std::vector<int>& samecl_index, int new_c1, int new_c2,
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
      l[0] = startingH->dish_llikelihood(new_c1_ind) + startingH->dish_llikelihood(new_c2_ind);
      new_c1_ind.remove(samecl_index[i]);
      new_c2_ind.push_back(samecl_index[i]);
      // update nTableRest (does not change)
      l[1] = startingH->dish_llikelihood(new_c1_ind) + startingH->dish_llikelihood(new_c2_ind);
      new_c2_ind.pop_back();
      new_c1_ind.push_back(samecl_index[i]);
    } else {
      P[0] = new_c1_ind.size();
      P[1] = new_c2_ind.size()-1;
      l[1] = startingH->dish_llikelihood(new_c1_ind) + startingH->dish_llikelihood(new_c2_ind);
      new_c2_ind.remove(samecl_index[i]);
      new_c1_ind.push_back(samecl_index[i]);
      l[0] = startingH->dish_llikelihood(new_c1_ind) + startingH->dish_llikelihood(new_c2_ind);
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
      if(startingH->tables_dishes->cluster_assignment[ samecl_index[i] ] == new_c1){
        lqratio += logpr[0];
      } else if(startingH->tables_dishes->cluster_assignment[ samecl_index[i] ] == new_c2) {
        lqratio += logpr[1];
      } else {
        cout << "ProposalSM_TableDishes_iter: Problem!";
        tables_dishes->Print_Partition();
      }
    }
  }
  return;
}

void Hdp::LaunchSM_TableDishes(LPHdp startingH, int i1, int i2, std::vector<int> samecl_index, int c1, int c2, std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, int ngibbs, std::default_random_engine& gen){
  // will NOT modify the current partition

  // samecl_index contain the indices (from 1 to n)
  // new_c1_ind, new_c2_ind contain the indices (from 1 to n) divided into the two clusters
  // samecl_newcl contains either new_c1 or new_c2 for the elements of samecl_index, with the same order

  // std::random_device myRandomDevice;
  // unsigned seed = myRandomDevice();
  // std::default_random_engine gen(seed);
  // gen(seed);

  int new_c1, new_c2;
  if(c1 == c2){
    new_c1 = c1;
    new_c2 = startingH->tables_dishes->K; // 0-indexing
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
      ProposalSM_TableDishes_iter(true, startingH, samecl_index, new_c1, new_c2, new_c1_ind, new_c2_ind, samecl_newcl, lqratio_ignored, gen);
    } // end loop t
  } 
  return;
}

bool Hdp::SampleSM_CostTable(int r, int ngibbs, double temp, std::default_random_engine& eng, double eta_CT, double eta_TD, bool local_print){
  // do we really need to copy this into startingP? can we just pass this? is startingP modified into ProposalSM
  // LPHdp startingP = new Hdp(this);

  LPHdp proposal = new Hdp(this);
  double lqratio, lpratio, llikratio;
  proposal->ProposalSM_CostTable(this, r, ngibbs, lqratio, lpratio, llikratio, eng, eta_CT, eta_TD, local_print);
  double a = lqratio + temp * (lpratio + llikratio);
  a = exp(a);
  // if(r == 1){
  //   proposal->costumers_tables[r]->Print_Partition();
  //   proposal->costumers_tables[r]->Print_Y();
  //   cout << lqratio << " " << lpratio << " " << llikratio << endl;
  //   cout << exp(lqratio) << " " << exp(lpratio) << " " << exp(llikratio) << " "; 
  //   cout << a << endl; 
  // }
  
  a = min(1.0,a);
  uniform_real_distribution<double> distribution(0.0,1.0);

  bool accepted;
  if(distribution(eng) < a){ // distribution(generator) < a
    // change the current partition into proposal
    Copy_Hdp(proposal);
    accepted = true;
    // cout << "accepted!" << endl;
  } else {
    accepted = false;
    // cout << "rejected" << endl;
  }
  delete proposal;

  return accepted;
}

void Hdp::ProposalSM_CostTable(LPHdp startingH, int rest, int ngibbs, double& lqratio, double& lpratio, double& llikratio, std::default_random_engine& myRandomEngine, double eta_CT, double eta_TD, bool local_print){
  // Should modify proposal, not startingP, update the values of lqratio, lpratio, llikratio
  // It's updading proposal's tables_dishes, 
  // Proposal at the beginning is equal to startingH
  lqratio = 0;
  LPPartition startingP = startingH->costumers_tables[rest];

  if(startingP->nObs == 1){
    return;
  }
  
  std::vector<int> i12 = sample_without_replacement(2, startingP->nObs, myRandomEngine);
  int i1 = i12[0];
  int i2 = i12[1];
  int c1 = startingP->cluster_assignment[i1];
  int c2 = startingP->cluster_assignment[i2];
  if(c2 < c1){
    // this is temporary because not really correct
    int tmp;
    tmp = i1;
    i1 = i2;
    i2 = tmp;
    tmp = c1;
    c1 = c2;
    c2 = tmp;
  }
  if(local_print) cout << rest << " " << i1 << " " << i2 << " " << c1 << " " << c2 << endl;
  int dish1, dish2, new_c1, new_c2, new_dish1, new_dish2;
  if(c1 == c2){
    int c1MainID = Table_RestId_MainId(rest, c1, startingH->nTableRest);
    if(local_print) cout << rest << " " << c1MainID << " " << startingH->tables_dishes->nObs << endl;
    dish1 = startingH->tables_dishes->cluster_assignment[c1MainID];
    dish2 = dish1;
    // i2 will go in new cluster
  } else {
    int c1MainID = Table_RestId_MainId(rest, c1, startingH->nTableRest);
    int c2MainID = Table_RestId_MainId(rest, c2, startingH->nTableRest);
    if(local_print) cout << rest << " " << c1MainID << " " << c2MainID << " " << startingH->tables_dishes->nObs << endl;
    dish1 = startingH->tables_dishes->cluster_assignment[c1MainID];
    dish2 = startingH->tables_dishes->cluster_assignment[c2MainID];
    // everyone will go into dish1
  }
  if(local_print) cout << "done mainID" << endl;
  std::vector<double> pr;
  pr = startingH->get_pr(i2, rest);
  // std::vector<double> pr(startingH->tables_dishes->K+1);
  // for(int d = 0; d < startingH->tables_dishes->K; d++){
  //   pr[d] = 1; //(double)1/(startingH->tables_dishes->K+1);
  // }
  // pr[startingH->tables_dishes->K] = 1; //(double)1/(startingH->tables_dishes->K+1);

  if(c1 == c2){
    new_c1 = c1;
    new_c2 = startingH->costumers_tables[rest]->K;
    new_dish1 = dish1;
    std::discrete_distribution<int> distribution(pr.begin(), pr.end());
    new_dish2 = distribution(myRandomEngine);
    lqratio += -log(pr[new_dish2]);
  } else {
    new_c1 = c1;
    new_c2 = c1; // on purpose for LaunchSM and iter
    new_dish1 = dish1;
    new_dish2 = dish1;
    lqratio += +log(pr[dish2]);
  }
  // cout << "Sampled costumers: " << i1 << " " << i2 << endl;
  // cout << "MainID: " << startingP->indexing[i1] << " " << startingP->indexing[i2] << endl;
  // cout << "Corresponding clusters: " << c1 << " " << c2 << endl;
  // cout << "Corresponding ys: " << Y( startingP->indexing[i1] ) << " " << Y( startingP->indexing[i2] ) << endl;
  // cout << "New clusters: " << new_c1 << " " << new_c2 << endl;
  // cout << "New dishes: " << new_dish1 << " " << new_dish2 << endl;
  // cout << endl << endl;

  // samecl_index contains the other people in those cluster (but not i1, i2) localID
  std::vector<int> samecl_index = startingP->get_samecl_index(i1, i2, c1, c2);

  // launch will NOT modify the current partition, but update new_c1_ind, new_c2_ind,samecl_newcl
  std::list<int> new_c1_ind,new_c2_ind; // lists are optimized to remove elements in any position (use "remove")
  std::vector<int> samecl_newcl;
  int params[10] = {i1, i2, c1, c2, dish1, dish2, new_c1, new_c2, new_dish1, new_dish1};
  // LaunchSM_CostTable(startingH, rest, i1, i2, samecl_index, c1, c2, dish1, dish2, new_c1_ind, new_c2_ind, samecl_newcl, ngibbs, myRandomEngine);
  LaunchSM_CostTable(startingH, rest, params, samecl_index, new_c1_ind, new_c2_ind, samecl_newcl, ngibbs, myRandomEngine);
  
  if(c1 == c2){
    ProposalSM_CostTable_iter(true, startingH, rest, params, samecl_index, new_c1_ind, new_c2_ind, samecl_newcl, lqratio, myRandomEngine);
    int nc1 = new_c1_ind.size();
    int nc2 = new_c2_ind.size();
    lpratio = log(eta_CT) + lbeta(nc1,nc2);
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
    costumers_tables[rest]->Split(c1, index1, index2, nc1, nc2, eta_CT);
    nTableRest[rest] += 1;
    delete[] index1;
    delete[] index2;
    // int i1_mainID = Costumer_RestId_MainId(rest, i1, nObsRest);
    // int cluster_newtable = costumers_dishes->cluster_assignment[i1_mainID]; // AHHHHH Not supposed to do this
    int newtable_id = nTableRest[rest] -1; // because of 0-indexing

    llikratio = alpha1 * (- tables_dishes->get_logprior(eta_TD));
    tables_dishes->AddTable_TableDish(rest, newtable_id, new_dish2, nTableRest, eta_TD);
    llikratio += alpha1 * (tables_dishes->get_logprior(eta_TD));
    
    if(new_dish2 == new_dish1){
      llikratio += 0;
    } else {
      if(new_dish2 == startingH->tables_dishes->K){
        llikratio += alpha1 * (- startingH->dish_llikelihood(dish1)); // startingH is not changed
        llikratio += alpha1 * (dish_llikelihood(new_dish1) + dish_llikelihood(new_dish2));
      } else {
        llikratio += alpha1 * (- startingH->dish_llikelihood(dish1) - startingH->dish_llikelihood(new_dish2)); // startingH is not changed
        llikratio += alpha1 * (dish_llikelihood(new_dish1) + dish_llikelihood(new_dish2));
      }
    }
  } else {
    // maybe at the end: merge everyone into one cluster
    lpratio = - log(eta_CT) - lbeta(startingP->cluster_config[c1], startingP->cluster_config[c2]);

    ProposalSM_CostTable_iter(false, startingH, rest, params, samecl_index, new_c1_ind, new_c2_ind, samecl_newcl, lqratio, myRandomEngine);

    costumers_tables[rest]->Merge(c1,c2, eta_CT);
    nTableRest[rest] += -1;
    int c2_mainID = Table_RestId_MainId(rest, c2, nTableRest);
    bool flag = false;
    if(tables_dishes->cluster_config[dish2] == 1){
      // then K is reduced and the new_dish names shift
      flag = true;
      if(dish1 > dish2){
        new_dish1 = dish1 - 1;
      }
    }
    llikratio = alpha1 * (- tables_dishes->get_logprior(eta_TD));
    tables_dishes->RemoveTable_TableDish(rest, c2_mainID, eta_TD);    
    llikratio += alpha1 * (tables_dishes->get_logprior(eta_TD));

    if(dish1 == dish2){
      llikratio += 0;
    } else {
      llikratio += alpha1 * (- startingH->dish_llikelihood(dish1) - startingH->dish_llikelihood(dish2)); // startingH is not changed
      // what if the name of the new_dish change because we removed one element?
      if(flag){
        // cout << "flag" << endl;
        llikratio += alpha1 * (dish_llikelihood(new_dish1));
      } else {
        llikratio += alpha1 * (dish_llikelihood(new_dish1) + dish_llikelihood(dish2));
      }
    }
  }
  dishes_thetas = get_dishes_thetas();
  return;
}

void Hdp::LaunchSM_CostTable(LPHdp startingH, int rest, int params[], std::vector<int> samecl_index, 
  std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, int ngibbs, std::default_random_engine& gen){
  // will NOT modify the current partition
  int i1, i2, c1, c2, dish1, dish2, new_c1, new_c2;//, new_dish1, new_dish2;
  i1 = params[0];
  i2 = params[1];
  c1 = params[2];
  c2 = params[3];
  dish1 = params[4];
  dish2 = params[5];
  if(c1 == c2){
    new_c1 = params[6];
    new_c2 = params[7];
    // new_dish1 = params[8];
    // new_dish2 = params[9];
  } else {
    new_c1 = c1;
    new_c2 = c2;
    // new_dish1 = dish1;
    // new_dish2 = dish2;
  }
  // samecl_index contain the indices (from 1 to n)
  // new_c1_ind, new_c2_ind contain the indices (from 1 to n - localID) divided into the two clusters
  // samecl_newcl contains either new_c1 or new_c2 for the elements of samecl_index, with the same order

  new_c1_ind.push_back(i1); // new_c1 is the future cluster of i1
  new_c2_ind.push_back(i2);

  if(samecl_index.size() > 0){
    std::bernoulli_distribution bern(0.5);
    samecl_newcl.resize(samecl_index.size());
    // randomly assign each element of samecl into one cluster or the other
    for(int i = 0; i < samecl_index.size(); i++){
      if(bern(gen)){  // bern(gen) returns bool
        new_c1_ind.push_back(samecl_index[i]); // becomes new_c1
        samecl_newcl[i] = new_c1;
        // cout << "(" << samecl_index[i] << "," << new_c1 << ") ";
      } else {
        new_c2_ind.push_back(samecl_index[i]); // new_c2
        samecl_newcl[i] = new_c2;
        // cout << "(" << samecl_index[i] << "," << new_c2 << ") ";
      }
    }
    // double y1 = Y(Costumer_RestId_MainId(rest, i1, nObsRest));
    // double y2 = Y(Costumer_RestId_MainId(rest, i2, nObsRest));
    // for(int i = 0; i < samecl_index.size(); i++){
    //   if( Y(Costumer_RestId_MainId(rest, samecl_index[i], nObsRest)) == y1){
    //     new_c1_ind.push_back(samecl_index[i]); // becomes new_c1
    //     samecl_newcl[i] = new_c1;
    //   } else {
    //     if(Y(Costumer_RestId_MainId(rest, samecl_index[i], nObsRest)) == y2){
    //       new_c2_ind.push_back(samecl_index[i]); // new_c2
    //       samecl_newcl[i] = new_c2;
    //     } else {
    //       cout << "ops" << endl;
    //     }
    //   }
    // }
    // cout << endl;
    double lqratio_ignored = 0;
    for(int t = 0; t < ngibbs; t++){
      ProposalSM_CostTable_iter(true, startingH, rest, params, samecl_index, new_c1_ind, new_c2_ind, samecl_newcl, lqratio_ignored, gen);
    } // end loop t
  } 
  return;
}

// this is better when the different restaurant are mixed together (in the wrong way) in tables_dishes
void Hdp::ProposalSM_CostTable_iter(bool SplitMerge, LPHdp startingH, int rest, int params[], std::vector<int> samecl_index,
  std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, double& lqratio, std::default_random_engine& myRandomEngine){
  int c1, c2, dish1, dish2; // i1, i2, 
  int new_c1, new_c2; //, new_dish1, new_dish2;
  // i1 = params[0];
  // i2 = params[1];
  c1 = params[2];
  c2 = params[3];
  dish1 = params[4];
  dish2 = params[5];
  if(c1 == c2){
    new_c1 = params[6];
    new_c2 = params[7];
    // new_dish1 = params[8];
    // new_dish2 = params[9];
  } else {
    new_c1 = c1;
    new_c2 = c2;
    // new_dish1 = dish1;
    // new_dish2 = dish2;
  }

  arma::vec P(2);
  arma::vec l(2);
  arma::vec logpr, pr;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  for(int i = 0; i < samecl_index.size(); i++){
    if(samecl_newcl[i] == new_c1){
      P[0] = new_c1_ind.size()-1;
      P[1] = new_c2_ind.size();
      l[0] = table_llikelihood(new_c1_ind, rest) + table_llikelihood(new_c2_ind, rest);
      new_c1_ind.remove(samecl_index[i]);
      new_c2_ind.push_back(samecl_index[i]);
      l[1] = table_llikelihood(new_c1_ind, rest) + table_llikelihood(new_c2_ind, rest);
      new_c2_ind.pop_back();
      new_c1_ind.push_back(samecl_index[i]);
    } else {
      P[0] = new_c1_ind.size();
      P[1] = new_c2_ind.size()-1;
      l[1] = table_llikelihood(new_c1_ind, rest) + table_llikelihood(new_c2_ind, rest);
      new_c2_ind.remove(samecl_index[i]);
      new_c1_ind.push_back(samecl_index[i]);
      l[0] = table_llikelihood(new_c1_ind, rest) + table_llikelihood(new_c2_ind, rest);
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
      if(startingH->costumers_tables[rest]->cluster_assignment[ samecl_index[i] ] == new_c1){
        lqratio += logpr[0];
      } else if(startingH->costumers_tables[rest]->cluster_assignment[ samecl_index[i] ] == new_c2) {
        lqratio += logpr[1];
      } else {
        cout << "ProposalSM_CostTable_iter Problem!"; //////////
      }
    }
  }
  return;
}

// void Hdp::SampleGibbs_CostTable(int r, std::default_random_engine& eng){
//   for(int i = 0; i < costumers_tables[r]->nObs; i++){
//     ConditionalGibbs_CostTable(r, i, eng);
//   }
//   return;
// }

// void Hdp::ConditionalGibbs_CostTable(int r, int i, std::default_random_engine& eng){
//   LPPartition restaurant = costumers_tables[r];
//   std::vector<double> newdish_pr;
//   arma::vec P, l, logpr;
//   int dish, tMainID, iMainID, itable, itableMainID, idish, ncl;
//   iMainID = Costumer_RestId_MainId(r, i, nObsRest);
//   itable = restaurant->cluster_assignment[i];
//   itableMainID = Table_RestId_MainId(r, itable, nTableRest);
//   idish = tables_dishes->cluster_assignment[ itableMainID ];
//   ncl = restaurant->cluster_config[itable];
//   if(ncl == 1){
//     P = zeros<vec>(restaurant->K);
//     l = zeros<vec>(restaurant->K);
//     for(int t = 0; t < restaurant->K;t++){
//       if(t < itable){
//         tMainID = Table_RestId_MainId(r, t, nTableRest);
//         dish = tables_dishes->cluster_assignment[ tMainID ];
//         double pp,ll;
//         pp = restaurant->cluster_config[t]; 
//         if(dish == idish){  
//           // I know there are other tables in this dish (because I am not doing mine), so I should be fine
//           ll = dish_llikelihood(dish)-dish_llikelihood_minusCost(dish,iMainID);
//         } else {
//           ll = dish_llikelihood_plusCost(dish,iMainID)-dish_llikelihood(dish);
//         }
//         P[t] = pp;
//         l[t] = ll;
//       } else if(t > itable){
//         tMainID = Table_RestId_MainId(r, t, nTableRest);
//         dish = tables_dishes->cluster_assignment[ tMainID ];
//         double pp,ll;
//         pp = restaurant->cluster_config[t]; 
//         if(dish == idish){  
//           // I know there are other tables in this dish (because I am not doing mine), so I should be fine
//           ll = dish_llikelihood(dish)-dish_llikelihood_minusCost(dish,iMainID);
//         } else {
//           ll = dish_llikelihood_plusCost(dish,iMainID)-dish_llikelihood(dish);
//         }
//         P[t-1] = pp;
//         l[t-1] = ll;
//       }
//     }
//     P[restaurant->K -1] = eta;
//     l[restaurant->K -1] = newdish_llikelihood(iMainID, idish, newdish_pr);
    
//     logpr = log(P) + l;
//     logpr = logpr - max(logpr);
//     logpr = logpr - log(sum(exp(logpr)));

//     std::vector<double> pr(logpr.n_rows);
//     for(int k = 0; k < logpr.n_rows;k++){
//       pr[k] = exp(logpr[k]);
//       // cout << pr[k] << " ";
//     }
//     // cout << endl;
//     std::discrete_distribution<int> distribution(pr.begin(), pr.end());
//     int sampled = distribution(eng);

//     if(sampled < restaurant->K){
//       if(sampled >= itable){
//         sampled = sampled + 1;
//       } 
//       restaurant->Merge(itable, sampled); // merge will fix which one is larger
//       // Since i was alone, we eliminate the table
//       tables_dishes->RemoveTable_TableDish(r, itableMainID);  
//       nTableRest[r] += -1;  
//     } else {
//       // i is already alone, so no need to split it. but we need fix the table
//       std::discrete_distribution<int> distribution(newdish_pr.begin(), newdish_pr.end());
//       int cluster_newtable = distribution(eng);
//       // if cluster_newtable = tables_dishes->K then it's a new cluster
//       // cout << "cluster_newtable " << cluster_newtable << endl;

//       if(tables_dishes->cluster_config[idish] > 1){
//         int *index1, *index2;
//         index1 = new int[ tables_dishes->cluster_config[idish] -1];
//         index2 = new int[1];
//         int count = 0;
//         for(int t = 0; t < tables_dishes->cluster_config[idish]; t++){
//           int tmp = tables_dishes->clusters[itable][t];
//           if(tmp != itableMainID){
//             index1[count] = tmp;
//             count++;
//           }
//         }
//         index2[0] = itableMainID;
//         if(cluster_newtable == tables_dishes->K){
//           // new dish
//           tables_dishes->Split(idish, index1, index2, tables_dishes->cluster_config[idish]-1, 1);
//         } else {
//           tables_dishes->Split_and_Merge(idish, index1, index2, tables_dishes->cluster_config[idish]-1, 1, idish, cluster_newtable);  
//         }
//         delete[] index1;
//         delete[] index2;
//       } else {
//         if(cluster_newtable != tables_dishes->K){
//           tables_dishes->Merge(idish, cluster_newtable);
//         } else {
//           // do nothing: you would be put in a new dish but you are already alone in your dish
//         }
//       }
//     }
//   } else {
//     P = zeros<vec>(restaurant->K+1);
//     l = zeros<vec>(restaurant->K+1);
//     for(int k = 0; k < restaurant->K;k++){
//       if(k == itable){
//         P[k] = restaurant->cluster_config[k]-1; // because we don't need to count i
//         // if i is the only one at that table, it should not be considered among the options
//       } else {
//         P[k] = restaurant->cluster_config[k];
//       }
//       tMainID = Table_RestId_MainId(r, k, nTableRest);
//       dish = tables_dishes->cluster_assignment[ tMainID ];
//       if(dish == idish){  
//         l[k] = dish_llikelihood(dish)-dish_llikelihood_minusCost(dish,iMainID);
//       } else {
//         l[k] = dish_llikelihood_plusCost(dish,iMainID)-dish_llikelihood(dish);
//       }
//     }
//     P[restaurant->K] = eta;
//     l[restaurant->K] = newdish_llikelihood(iMainID, idish, newdish_pr);
  
//     arma::vec logpr = log(P) + l;
//     logpr = logpr - max(logpr);
//     logpr = logpr - log(sum(exp(logpr)));
  
//     std::vector<double> pr(logpr.n_rows);
//     for(int k = 0; k < logpr.n_rows;k++){
//       pr[k] = exp(logpr[k]);
//       // cout << pr[k] << " ";
//     }
//     // cout << endl;
//     std::discrete_distribution<int> distribution(pr.begin(), pr.end());
//     int sampled = distribution(eng);
//     // cout << "Sampled " << sampled << endl; 
//     if(sampled < restaurant->K){
//       int *index1, *index2;
//       index1 = new int[ncl-1];
//       index2 = new int[1];
//       int count = 0;
//       for(int j = 0; j < ncl; j++){
//         int tmp = restaurant->clusters[itable][j];
//         if(tmp != i){
//           index1[count] = tmp;  
//           count++;
//         }
//       }
//       index2[0] = i;
//       // we can do a splitmerge for restaurant
//       restaurant->Split_and_Merge(itable, index1, index2, ncl-1, 1, itable, sampled);
//       delete[] index1;
//       delete[] index2;
//     } else {
//       std::discrete_distribution<int> distribution(newdish_pr.begin(), newdish_pr.end());
//       int *index1, *index2;
//       index1 = new int[ncl-1];
//       index2 = new int[1];
//       int count = 0;
//       for(int j = 0; j < ncl; j++){
//         int tmp = restaurant->clusters[itable][j];
//         if(tmp != i){
//           index1[count] = tmp;  
//           count++;
//         }
//       }
//       index2[0] = i;
//       restaurant->Split(itable, index1, index2, ncl-1, 1);
//       // we sample the new dish according to the probability of newdish_llikelihood
//       delete[] index1;
//       delete[] index2;
//       int cluster_newtable = distribution(eng);
//       nTableRest[r] += 1;
//       tables_dishes->AddTable_TableDish(r, restaurant->K-1, cluster_newtable, nTableRest);
//     }
//   }
//   return;
// }
// double Hdp::newdish_llikelihood(int i, int current_dish, std::vector<double>& pr){
//   // i is MainID, t is mainID
//   // when we remove i from its table, and (in case i is alone) its table from all the tables,
//   // if its table was the only one with that dish, we need to remove the dish
//   // int tMainID, tRestID, rest;
//   arma::vec P;
//   arma::vec l;
  
//   if(tables_dishes->cluster_config[current_dish] == 1){
//     P = zeros<vec>(tables_dishes->K);
//     l = zeros<vec>(tables_dishes->K);
//     for(int k=0; k < tables_dishes->K; k++){
//       if(k < current_dish){
//         P[k] = tables_dishes->cluster_config[k];
//         l[k] = dish_llikelihood_plusCost(k,i)-dish_llikelihood(k);
//       } else if(k > current_dish){
//         P[k-1] = tables_dishes->cluster_config[k];
//         l[k-1] = dish_llikelihood_plusCost(k,i)-dish_llikelihood(k);
//       }
//     }
//     P[tables_dishes->K -1] = gamma0;
//     l[tables_dishes->K -1] = cost_llikelihood(i);
//   } else {
//     P = zeros<vec>(tables_dishes->K+1);
//     l = zeros<vec>(tables_dishes->K+1);
//     for(int k=0; k < tables_dishes->K; k++){
//       P[k] = tables_dishes->cluster_config[k];
//       if(k == current_dish){
//         l[k] = dish_llikelihood(k)-dish_llikelihood_minusCost(k,i);
//       } else {
//         l[k] = dish_llikelihood_plusCost(k,i)-dish_llikelihood(k);
//       }
//     }
//     P[tables_dishes->K] = gamma0;
//     l[tables_dishes->K] = cost_llikelihood(i);
//   }
//   P = P/arma::accu(P);
//   pr.resize(P.n_rows);
//   for(int k=0; k < P.n_rows; k++){
//     pr[k] = P[k]*exp(l[k]);
//   }
//   return(arma::accu(log(P) + l)); 
// }

// double Hdp::cost_llikelihood(int i){
//   // i is MainID
//   int n = 1;
//   arma::vec y_cl = zeros<vec>( 1 );
//   y_cl(0) = Y( i );
//   double ybar = arma::mean(y_cl);
//   double k_n = k_0 + n;
//   // double mu_n = (k_0*mu_0 + n*ybar)/k_n;
//   double alpha_n = alpha_0 + n/2;
//   // double beta_n = beta_0 + sum(pow(y_cl-ybar,2)) + k_0 * n * pow(ybar-mu_0,2)/(2*k_n);
//   double beta_n = beta_0 + k_0 * n * pow(ybar-mu_0,2)/(2*k_n);
//   double loglike = lgamma(alpha_n) - lgamma(alpha_0) + alpha_0 * log(beta_0) - alpha_n * log(beta_n);
//   loglike += 0.5*(log(k_0)-log(k_n)) -n/2*log(2*M_PI);
//   return(loglike); 
// }
// double Hdp::dish_llikelihood_plusCost(int dish_k, int i){
//   // from dish_k we need to find all the tables, and all the costumers
//   // for this, we need to keep updated only nTableRest
//   // i should be the MainID for that costumer
//   int table_mainID, rest, table_restID, cost_restID, cost_mainID;
//   std::vector<int> costumers_dish_k;
//   for(int t = 0; t < tables_dishes->cluster_config[dish_k]; t++){
//     table_mainID = tables_dishes->clusters[dish_k][t];
//     Table_MainId_RestId(table_mainID, nTableRest, rest, table_restID);
//     for(int c = 0; c < costumers_tables[rest]->cluster_config[table_restID]; c++){
//       cost_restID = costumers_tables[rest]->clusters[table_restID][c];
//       cost_mainID = costumers_tables[rest]->indexing[cost_restID];
//       // cost_mainID = Costumer_RestId_MainId(rest, cost_restID, nObsRest);
//       costumers_dish_k.push_back(cost_mainID);
//     }
//   }
//   costumers_dish_k.push_back(i);
//   int n = costumers_dish_k.size();
//   if(n == 0){
//     cout << "dish_llikelihood_plusCost(dish_k,i) called with empty cluster!" << endl;
//   }
//   arma::vec y_cl = zeros<vec>( n );
//   for(int c = 0; c < n; c++){
//     y_cl(c) = Y( costumers_dish_k[c] );
//   }
//   double ybar = arma::mean(y_cl);
//   double k_n = k_0 + n;
//   // double mu_n = (k_0*mu_0 + n*ybar)/k_n;
//   double alpha_n = alpha_0 + n/2;
//   double beta_n = beta_0 + sum(pow(y_cl-ybar,2)) + k_0 * n * pow(ybar-mu_0,2)/(2*k_n);
//   double loglike = lgamma(alpha_n) - lgamma(alpha_0) + alpha_0 * log(beta_0) - alpha_n * log(beta_n);
//   loglike += 0.5*(log(k_0)-log(k_n)) -n/2*log(2*M_PI);
//   return(loglike); 
// }
// double Hdp::dish_llikelihood_minusCost(int dish_k, int i){
//   // from dish_k we need to find all the tables, and all the costumers
//   // for this, we need to keep updated only nTableRest
//   // i should be the MainID for that costumer
//   int table_mainID, rest, table_restID, cost_restID, cost_mainID;
//   std::vector<int> costumers_dish_k;
//   for(int t = 0; t < tables_dishes->cluster_config[dish_k]; t++){
//     table_mainID = tables_dishes->clusters[dish_k][t];
//     Table_MainId_RestId(table_mainID, nTableRest, rest, table_restID);
//     for(int c = 0; c < costumers_tables[rest]->cluster_config[table_restID]; c++){
//       cost_restID = costumers_tables[rest]->clusters[table_restID][c];
//       cost_mainID = costumers_tables[rest]->indexing[cost_restID];
//       // cost_mainID = Costumer_RestId_MainId(rest, cost_restID, nObsRest);
//       if(cost_mainID != i){
//         costumers_dish_k.push_back(cost_mainID);  
//       }
//     }
//   }
//   int n = costumers_dish_k.size();
//   if(n == 0){
//     cout << "dish_llikelihood_minusCost(dish_k,i) called with empty cluster!" << endl;
//   }
//   arma::vec y_cl = zeros<vec>( n );
//   for(int c = 0; c < n; c++){
//     y_cl(c) = Y( costumers_dish_k[c] );
//   }
//   double ybar = arma::mean(y_cl);
//   double k_n = k_0 + n;
//   // double mu_n = (k_0*mu_0 + n*ybar)/k_n;
//   double alpha_n = alpha_0 + n/2;
//   double beta_n = beta_0 + sum(pow(y_cl-ybar,2)) + k_0 * n * pow(ybar-mu_0,2)/(2*k_n);
//   double loglike = lgamma(alpha_n) - lgamma(alpha_0) + alpha_0 * log(beta_0) - alpha_n * log(beta_n);
//   loglike += 0.5*(log(k_0)-log(k_n)) -n/2*log(2*M_PI);
//   return(loglike); 
// }
// double Hdp::dish_llikelihood_plusminusCosts(int dish_k, std::list<int> index1, std::list<int> index2, int rest){
//   // index1, index2 contain local IDs
//   // I assume index1, index2 are contained in the same rest (it's fine)
//   // I need to assume index1 and index2 don't have elements in common
//   int table_mainID, trest, table_restID, cost_restID, cost_mainID;
//   bool found_index2;
//   bool found_index1 = false;
//   std::vector<int> costumers_dish_k;
//   if(dish_k < tables_dishes->K){
//     for(int t = 0; t < tables_dishes->cluster_config[dish_k]; t++){
//       table_mainID = tables_dishes->clusters[dish_k][t];
//       Table_MainId_RestId(table_mainID, nTableRest, trest, table_restID);
//       // cout << "dish_llikelihood_plusminusCosts " << trest << ": ";
//       if(trest != rest){
//         for(int c = 0; c < costumers_tables[trest]->cluster_config[table_restID]; c++){
//           cost_restID = costumers_tables[trest]->clusters[table_restID][c];
//           cost_mainID = costumers_tables[trest]->indexing[cost_restID];
//           // cost_mainID = Costumer_RestId_MainId(trest, cost_restID, nObsRest);
//           costumers_dish_k.push_back(cost_mainID);
//         }
//       } else {
//         for(int c = 0; c < costumers_tables[trest]->cluster_config[table_restID]; c++){
//           cost_restID = costumers_tables[trest]->clusters[table_restID][c];
//           found_index1 = (std::find(index1.begin(), index1.end(), cost_restID) != index1.end());
//           found_index2 = (std::find(index2.begin(), index2.end(), cost_restID) != index2.end());

//           if(found_index1 || found_index2){
//             // do not add it
//             // cout << cost_restID << " ";
//           } else {
//             cost_mainID = costumers_tables[trest]->indexing[cost_restID];
//             // cost_mainID = Costumer_RestId_MainId(trest, cost_restID, nObsRest);
//             costumers_dish_k.push_back(cost_mainID);
//           }
//         }
//       }
//       // cout << endl;
//     }
//   } 
//   for(std::list<int>::iterator it = index1.begin(); it != index1.end(); ++it){
//     cost_mainID = costumers_tables[rest]->indexing[*it];
//     // cost_mainID = Costumer_RestId_MainId(rest, *it, nObsRest);
//     costumers_dish_k.push_back(cost_mainID);
//   }
//   int n = costumers_dish_k.size();
//   arma::vec y_cl = zeros<vec>( n );
//   for(int c = 0; c < n; c++){
//     y_cl(c) = Y( costumers_dish_k[c] );
//   }
//   double ybar = arma::mean(y_cl);
//   // cout << ybar << " " << n<< endl;
//   double k_n = k_0 + n;
//   // double mu_n = (k_0*mu_0 + n*ybar)/k_n;
//   double alpha_n = alpha_0 + n/2;
//   double beta_n = beta_0 + sum(pow(y_cl-ybar,2)) + k_0 * n * pow(ybar-mu_0,2)/(2*k_n);
//   double loglike = lgamma(alpha_n) - lgamma(alpha_0) + alpha_0 * log(beta_0) - alpha_n * log(beta_n);
//   loglike += 0.5*(log(k_0)-log(k_n)) -n/2*log(2*M_PI);
//   return(loglike); 
// }

std::vector<double> Hdp::get_pr(int i, int rest){
  double Yi = Y( costumers_tables[rest]->indexing[i] );
  std::vector<double> pr(tables_dishes->K+1);
  std::vector<double> pr_norm(tables_dishes->K+1);
  double Yd;
  double C = 0;
  for(int d = 0; d < tables_dishes->K; d++){
    Yd = dish_avg(d);
    pr[d] = exp(- std::abs(Yi - Yd));
    C += pr[d];
  }
  pr[tables_dishes->K] = 1;
  pr_norm[tables_dishes->K] = 1/(C+1);
  for(int d = 0; d < tables_dishes->K; d++){
    pr_norm[d] = pr[d]/(C+1);
  }
  return(pr_norm);
}

double Hdp::dish_avg(int dish_k){
  // from dish_k we need to find all the tables, and all the costumers
  // for this, we need to keep updated only nTableRest
  int table_mainID, rest, table_restID, cost_restID, cost_mainID;
  std::vector<int> costumers_dish_k;
  for(int t = 0; t < tables_dishes->cluster_config[dish_k]; t++){
    table_mainID = tables_dishes->clusters[dish_k][t];
    Table_MainId_RestId(table_mainID, nTableRest, rest, table_restID);
    for(int c = 0; c < costumers_tables[rest]->cluster_config[table_restID]; c++){
      cost_restID = costumers_tables[rest]->clusters[table_restID][c];
      cost_mainID = costumers_tables[rest]->indexing[cost_restID];
      // cost_mainID = Costumer_RestId_MainId(rest, cost_restID, nObsRest);
      costumers_dish_k.push_back(cost_mainID);
      
      // cout << cost_mainID << " ";
    }
  }
  // cout << endl;
  
  int n = costumers_dish_k.size();
  if(n == 0){
    cout << "dish_llikelihood(dish_k) called with empty cluster!" << endl;
  }
  arma::vec y_cl = zeros<vec>( n );
  for(int c = 0; c < n; c++){
    y_cl(c) = Y( costumers_dish_k[c] );
  }
  double ybar = arma::mean(y_cl);
  return(ybar); 
}

double Hdp::dish_llikelihood(int dish_k, bool local_print){
  // from dish_k we need to find all the tables, and all the costumers
  // for this, we need to keep updated only nTableRest
  if(local_print)
    cout << "Hdp::dish_llikelihood, cost_mainID: ";
  int table_mainID, rest, table_restID, cost_restID, cost_mainID;
  std::vector<int> costumers_dish_k;
  for(int t = 0; t < tables_dishes->cluster_config[dish_k]; t++){
    table_mainID = tables_dishes->clusters[dish_k][t];
    Table_MainId_RestId(table_mainID, nTableRest, rest, table_restID);
    for(int c = 0; c < costumers_tables[rest]->cluster_config[table_restID]; c++){
      cost_restID = costumers_tables[rest]->clusters[table_restID][c];
      cost_mainID = costumers_tables[rest]->indexing[cost_restID];
      // cost_mainID = Costumer_RestId_MainId(rest, cost_restID, nObsRest);
      costumers_dish_k.push_back(cost_mainID);
      if(local_print)
        cout << cost_mainID << " ";
    }
  }
  if(local_print)
    cout << endl;
  
  int n = costumers_dish_k.size();
  if(n == 0){
    cout << "dish_llikelihood(dish_k) called with empty cluster!" << endl;
  }
  arma::vec y_cl = zeros<vec>( n );
  for(int c = 0; c < n; c++){
    y_cl(c) = Y( costumers_dish_k[c] );
  }
  if(local_print){
    cout << "Hdp::dish_llikelihood, Ys: ";
    cout << y_cl.t() << endl;
  }
  
  double ybar = arma::mean(y_cl);
  double k_n = k_H + n;
  // double mu_n = (k_H*mu_H + n*ybar)/k_n;
  double alpha_n = alpha_H + n/2;
  double beta_n = beta_H + sum(pow(y_cl-ybar,2)) + k_H * n * pow(ybar-mu_H,2)/(2*k_n);
  double loglike = lgamma(alpha_n) - lgamma(alpha_H) + alpha_H * log(beta_H) - alpha_n * log(beta_n);
  loglike += 0.5*(log(k_H)-log(k_n)) -n/2*log(2*M_PI);
  return(loglike); 
}

double Hdp::dish_llikelihood(std::list<int> index){
  // index contains tables
  // from dish_k we need to find all the tables, and all the costumers
  // for this, we need to keep updated only nTableRest
  int table_mainID, rest, table_restID, cost_restID, cost_mainID;
  std::vector<int> costumers_dish_k;
  for(std::list<int>::iterator it = index.begin(); it != index.end(); ++it){
    table_mainID = *it;
    Table_MainId_RestId(table_mainID, nTableRest, rest, table_restID);
    for(int c = 0; c < costumers_tables[rest]->cluster_config[table_restID]; c++){
      cost_restID = costumers_tables[rest]->clusters[table_restID][c];
      cost_mainID = costumers_tables[rest]->indexing[cost_restID];
      // cost_mainID = Costumer_RestId_MainId(rest, cost_restID, nObsRest);
      costumers_dish_k.push_back(cost_mainID);
    }
  }
  int n = costumers_dish_k.size();
  if(n == 0){
    cout << "dish_llikelihood(index) called with empty list!" << endl;
  }
  arma::vec y_cl = zeros<vec>( n );
  for(int c = 0; c < n; c++){
    y_cl(c) = Y( costumers_dish_k[c] );
  }
  
  double ybar = arma::mean(y_cl);
  double k_n = k_H + n; 
  // double mu_n = (k_H*mu_H + n*ybar)/k_n; 
  double alpha_n = alpha_H + n/2; 
  double beta_n = beta_H + sum(pow(y_cl-ybar,2)) + k_H * n * pow(ybar-mu_H,2)/(2*k_n);
  double loglike = lgamma(alpha_n) - lgamma(alpha_H) + alpha_H * log(beta_H) - alpha_n * log(beta_n);
  loglike += 0.5*(log(k_H)-log(k_n)) -n/2*log(2*M_PI);
  return(loglike);
}

double Hdp::table_llikelihood(std::list<int> index, int rest){
  // remember: index now contains costumers
  int cost_restID, cost_mainID;
  std::vector<int> costumers_dish_k;
  for(std::list<int>::iterator it = index.begin(); it != index.end(); ++it){
    cost_restID = *it;
    cost_mainID = costumers_tables[rest]->indexing[cost_restID];
    // cost_mainID = Costumer_RestId_MainId(rest, cost_restID, nObsRest);
    costumers_dish_k.push_back(cost_mainID);
  }
  int n = costumers_dish_k.size();
  if(n == 0){
    cout << "table_llikelihood called with empty list!" << endl;
  }
  arma::vec y_cl = zeros<vec>( n );
  for(int c = 0; c < n; c++){
    y_cl(c) = Y( costumers_dish_k[c] );
  }
  double ybar = arma::mean(y_cl);
  double k_n = k_H + n; 
  // double mu_n = (k_H*mu_H + n*ybar)/k_n; 
  double alpha_n = alpha_H + n/2; 
  double beta_n = beta_H + sum(pow(y_cl-ybar,2)) + k_H * n * pow(ybar-mu_H,2)/(2*k_n);
  double loglike = lgamma(alpha_n) - lgamma(alpha_H) + alpha_H * log(beta_H) - alpha_n * log(beta_n);
  loglike += 0.5*(log(k_H)-log(k_n)) -n/2*log(2*M_PI);
  return(loglike);
}

double Hdp::get_logprior(double eta_CT, double eta_TD){
  double res = 0.0;
  for(int r = 0; r < nRest; r++){
    res += costumers_tables[r]->get_logprior(eta_CT);
  }
  res += tables_dishes->get_logprior(eta_TD);

  return res;
}

// double Hdp::get_logprior2(){
//   double res = 0.0;
//   for(int r = 0; r < nRest; r++){
//     res += costumers_tables[r]->get_logprior2(eta_CT);
//   }
//   res += tables_dishes->get_logprior2(eta_TD);

//   return res;
// }

double Hdp::get_loglike(){
  double res = 0.0;
  for(int k = 0; k < tables_dishes->K; k++){
    res += dish_llikelihood(k);
  }

  return res;
}

void Hdp::get_CostDish(LPPartition costdish){
  costdish->Initialize_CostDish(nObsRest, nTableRest, costumers_tables, tables_dishes);
  return;
}

std::vector<double> Hdp::get_dishes_thetas(){
  std::vector<double> dishes_thetas(tables_dishes->K,0);
  int x, y; // x, y are used for main ID
  int rest, table, cost; // rest, table, cost are used for restaurant ID
  vector<vector<int> > clusters_v(tables_dishes->K); // using a vector to make things easier
  for(int k = 0; k < tables_dishes->K; k++){
    for(int i = 0; i < tables_dishes->cluster_config[k]; i++){
      y = tables_dishes->clusters[k][i];
      Table_MainId_RestId(y, nTableRest, rest, table); // get rest and table
      for(int ii = 0; ii < costumers_tables[rest]->cluster_config[table]; ii++){
        cost = costumers_tables[rest]->clusters[table][ii];
        x = costumers_tables[rest]->indexing[cost];
        clusters_v[k].push_back( x ); // which people at that table  
        dishes_thetas[k] += Y(x);
      }
    }
    dishes_thetas[k] = dishes_thetas[k]/clusters_v[k].size();
  }
  
  return dishes_thetas;
}

// void Hdp::Split_CostTable(int r, int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int cluster_newtable){
//   // Split a cluster of one restaurant (of the partition costumers_tables[r])
//   // new_cluster1 and new_cluster2 should have the indices with the names of the elements.
//   // So not the indices of where they are in their current clusters, but their actual "names" 
//   if(cluster_newtable > tables_dishes->K){
//     cout << "Hdp::Split_CostTable ERROR, cluster_newtable should not be larger than the current value of K" << endl;
//     return;
//   }
//   costumers_tables[r]->Split(split_k, new_cluster1, new_cluster2, size1, size2);
//   // before calling Initialize_costDish we need to updated nTableRest because it's used inside the function
//   nTableRest[r] += 1;
//   // if(cluster_newtable == tables_dishes->K){
//   //   K += 1;
//   // }
//   // the new table in restaurant r needs to be assigned to a cluster in tables_dishes
//   int newtable_id = nTableRest[r] -1; // because of 0-indexing
//   tables_dishes->AddTable_TableDish(r, newtable_id, cluster_newtable, nTableRest);
//   return;
// }
// void Hdp::Split_TableDish(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2){
//   tables_dishes->Split(split_k, new_cluster1, new_cluster2, size1, size2);
//   // K = K + 1;
//   return;
// }


void Hdp::Merge_CostTable(int r, int merge_k1, int merge_k2, double eta_CT, double eta_TD){
  /*
  Merges two tables within one restaurant. 
  Needs to change nTableRest[r], costumers_tables[r], tables_dishes
  merge_k1,2 are localID in r.
  */
  int k_max = max(merge_k1, merge_k2);
  // int k_min = min(merge_k1, merge_k2);
  // I need the main ID for the old cluster before changing things
  int mainID_k_max = Table_RestId_MainId(r, k_max, nTableRest);

  costumers_tables[r]->Merge(merge_k1, merge_k2, eta_CT);
  // before calling Initialize_costDish we need to updated nTableRest because it's used inside the function
  nTableRest[r] -= 1;
  
  tables_dishes->RemoveTable_TableDish(r, mainID_k_max, eta_TD);
  // K = tables_dishes->K;
  return;
}

// void Hdp::Merge_TableDish(int merge_k1, int merge_k2){
//   tables_dishes->Merge(merge_k1, merge_k2);
//   // K = K - 1;
//   return;
// }

// void Hdp::SplitMerge_CostTable(int r, int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int merge_k1, int merge_k2){
//   // the idea is that newcluster1 is merged with kstar1 and newcluster2 is merged with kstar2.
//   if(merge_k1 != merge_k2){
//     int y_split_k = Table_RestId_MainId(r, split_k, nTableRest);
//     int y_merge_k1 = Table_RestId_MainId(r, merge_k1, nTableRest);
//     int y_merge_k2 = Table_RestId_MainId(r, merge_k2, nTableRest);

//     if((split_k == merge_k1) & (split_k != merge_k2)){ // leave new_cluster1 alone and just attempt to merge new_cluster2 with k_star_2
//       costumers_tables[r]->Split_and_Merge(split_k, new_cluster1, new_cluster2, size1, size2, merge_k1, merge_k2);
//       return;
//     } else {
//       cout << "Hdp::SplitMerge_CostTable ERROR, merge_k1 should be equal to split_k (other cases not implemented)" << endl;
//       return;
//     }
//   }
// }

// void Hdp::SplitMerge_TableDish(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int k_star_1, int k_star_2){
//   tables_dishes->Split_and_Merge(split_k, new_cluster1, new_cluster2, size1, size2, k_star_1, k_star_2);
//   return;
// }

// void Hdp::KSplit_CostTable(int r, int split_k, int num_splits, std::vector<std::vector<int> > indices, std::vector<int> ns, std::vector<int> cluster_newtables){
//   // Split a cluster of one restaurant (of the partition costumers_tables[r]) into K parts
//   // indices should have the names of the elements, not the indices in the clusters.
//   // So not the indices of where they are in their current clusters, but their actual "names" 

//   // cluster_newtables has length = num_splits - 1
//   // ns has length num_splits
//   costumers_tables[r]->KSplit(split_k, num_splits, indices, ns);
//   // before calling Initialize_costDish we need to updated nTableRest because it's used inside the function
//   nTableRest[r] += num_splits-1;
  
//   // the new table in restaurant r needs to be assigned to a cluster in tables_dishes
//   std::vector<int> newtables_id(num_splits-1,0);
//   for(int i = 0; i < num_splits-1; i++){
//     newtables_id[i] = nTableRest[r] - num_splits + 1 + i; 
//   }
//   tables_dishes->AddTables_TableDish(r, newtables_id, cluster_newtables, nTableRest);
//   // K = tables_dishes->K;
//   return;
// }

// void Hdp::KSplit_TableDish(int split_k, int num_splits, std::vector<std::vector<int> > indices, std::vector<int> ns){
//   tables_dishes->KSplit(split_k, num_splits, indices, ns);
//   // K = K + num_splits - 1;
//   return;
// }

void Hdp::format_rest_table(double* vector_rest, double* vector_table, int starting_index){
  for(int r = 0; r < nRest; r++){
    for(int i = 0; i < costumers_tables[r]->nObs; i++){
      vector_rest[starting_index + costumers_tables[r]->indexing[i]] = r + 1;
      vector_table[starting_index + costumers_tables[r]->indexing[i]] = costumers_tables[r]->cluster_assignment[i] + 1;
    }
  }
  return;
}