#include <vector>
#include <list>
#include <iostream>
#include <armadillo>
#include <algorithm>
#include <math.h>
#include "multires.h"
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
extern bool print_bool;
extern bool yLR_bool;
extern double alpha1;
extern double alpha2;
 
Multires::Multires(){
  nObs_lowres = 0;
  nObs_highres = 0;
  rest_config = std::vector<int>();
  rests = std::vector< std::vector<int> >();
  lowresPart = NULL;
  highresPart = NULL;
  return;
}

Multires::Multires(LPMultires initial_multires){
  nObs_lowres = initial_multires->nObs_lowres;
  nObs_highres = initial_multires->nObs_highres;
  rest_config = std::vector<int>(initial_multires->rest_config);
  rests = std::vector<std::vector<int> >(initial_multires->rests);
  lowresPart = new Partition(initial_multires->lowresPart);
  highresPart = new Hdp(initial_multires->highresPart);
  return;
}

Multires::~Multires(){
  nObs_lowres = 0;
  nObs_highres = 0;
  delete lowresPart;
  delete highresPart;
  lowresPart = NULL;
  highresPart = NULL;
  return;
}

void Multires::Copy_Multires(LPMultires initial_multires){
  delete lowresPart;
  delete highresPart;
  
  nObs_lowres = initial_multires->nObs_lowres;
  nObs_highres = initial_multires->nObs_highres;

  rest_config = std::vector<int>(initial_multires->rest_config);
  rests = std::vector<std::vector<int> >(initial_multires->rests);
  
  lowresPart = new Partition(initial_multires->lowresPart);
  highresPart = new Hdp(initial_multires->highresPart);
  return;
}

void Multires::Initialize_Multires(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, double eta_LR, double eta_CT, double eta_TD){
  // !! If there is 1 cluster in lowres, there should be 1 restaurant in highres?
  // at the moment Rests corresponds to the stable (CT)
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  rests = std::vector<std::vector<int> >(Rests); // this remembers the "stable" division of BG in CT
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }
  lowresPart = new Partition();
  highresPart = new Hdp();
  lowresPart->Initialize_Partition(nObs_lowres, eta_LR); 
  // everyone is in one cluster, so just one restaurant
  std::vector<std::vector<int> > CTinRest, Rest;
  CTinRest.resize(1);
  Rest.resize(1);
  CTinRest[0].resize(nObs_lowres);
  Rest[0].reserve(nObs_highres); // does not change elements
  for(int i = 0; i < nObs_lowres; i++){
    CTinRest[0][i] = i;
    for(int j = 0; j < Rests[i].size(); j++){
      Rest[0].push_back(Rests[i][j]); 
    }
  }
  // in this way, the only restaurant in Hdp will have indexing in Partition that is ordered by CT and not by mainID
  highresPart->Initialize_Hdp(nObs_highres, CTinRest, Rest, eta_CT, eta_TD);
  return;
}

void Multires::Initialize_Multires(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, LPPartition lowres, LPHdp highres){
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  rests = std::vector<std::vector<int> >(Rests);
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }
  lowresPart = new Partition(lowres);
  highresPart = new Hdp(highres);
  return;
}

void Multires::Initialize_Multires_nclLR_oneHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, double eta_LR, double eta_CT, double eta_TD){
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  rests = std::vector<std::vector<int> >(Rests);
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }
  lowresPart = new Partition();
  highresPart = new Hdp();
  lowresPart->Initialize_Partition_nclusters(nObs_lowres, eta_LR); 
  // everyone is in a different cluster, so nObs_lowres restaurants
  std::vector<std::vector<int> > CTinRest;
  CTinRest.resize(nObsL);
  for(int i = 0; i < nObs_lowres; i++){
    CTinRest[i].resize(1);
    CTinRest[i][0] = i;
  }
  highresPart->Initialize_Hdp(nObs_highres, CTinRest, Rests, eta_CT, eta_TD);
  return;
}

void Multires::Initialize_Multires_oneLR_nclHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, double eta_LR, double eta_CT, double eta_TD){
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  rests = std::vector<std::vector<int> >(Rests);
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }
  lowresPart = new Partition();
  highresPart = new Hdp();
  lowresPart->Initialize_Partition(nObs_lowres, eta_LR); 
  // everyone is in one cluster, so just one restaurant
  std::vector<std::vector<int> > CTinRest, Rest;
  CTinRest.resize(1);
  Rest.resize(1);
  CTinRest[0].resize(nObs_lowres);
  Rest[0].reserve(nObs_highres); // does not change elements
  for(int i = 0; i < nObs_lowres; i++){
    CTinRest[0][i] = i;
    for(int j = 0; j < Rests[i].size(); j++){
      Rest[0].push_back(Rests[i][j]); 
    }
  }
  highresPart->Initialize_Hdp_nclusters(nObs_highres, CTinRest, Rest, eta_CT, eta_TD);
  // by doing this we have one restaurant in highres, with all the costumers in different
  // tables and all the tables in different dishes.
  return;
}

void Multires::Initialize_Multires_oneLR_partHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, int part_init_ptr[], int* K_init_ptr, double eta_LR, double eta_CT, double eta_TD){
  // Rests represents the highres areas in each lowres area
  // set lowres with one cluster
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  rests = std::vector<std::vector<int> >(Rests);
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }
  lowresPart = new Partition();
  highresPart = new Hdp();
  lowresPart->Initialize_Partition(nObs_lowres, eta_LR); 
  // everyone is in one cluster, so just one restaurant
  std::vector<std::vector<int> > CTinRest, Rest;
  CTinRest.resize(1);
  Rest.resize(1);
  CTinRest[0].resize(nObs_lowres);
  Rest[0].reserve(nObs_highres); // does not change elements
  for(int i = 0; i < nObs_lowres; i++){
    CTinRest[0][i] = i;
    for(int j = 0; j < Rests[i].size(); j++){
      Rest[0].push_back(Rests[i][j]); 
    }
  }
  // this only works if rests has length 1:
  highresPart->Initialize_Hdp_part(nObs_highres, CTinRest, Rest, part_init_ptr, K_init_ptr, eta_CT, eta_TD);
  // by doing this we have one restaurant in highres, with all the costumers in tables assigned by part
  // and all the tables in different dishes.
  return;
}
void Multires::Initialize_Multires_nclLR_partHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, int part_init_ptr[], double eta_LR, double eta_CT, double eta_TD){
  // Rests represents the highres areas in each lowres area
  // set lowres with one cluster
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  rests = std::vector<std::vector<int> >(Rests);
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }
  lowresPart = new Partition();
  highresPart = new Hdp();
  lowresPart->Initialize_Partition_nclusters(nObs_lowres, eta_LR); 
  // everyone is in a different cluster, so nObs_lowres restaurants
  std::vector<std::vector<int> > CTinRest;
  CTinRest.resize(nObsL);
  for(int i = 0; i < nObs_lowres; i++){
    CTinRest[i].resize(1);
    CTinRest[i][0] = i;
  }
  // basically CTinRest tells me which CT are in each LR cluster (rest)
  // instead Rests tells me which BG belong to which CT
  highresPart->Initialize_Hdp_part_nrest(nObs_highres, CTinRest, Rests, part_init_ptr, eta_CT, eta_TD);
  return;
}

void Multires::Initialize_Multires_partLR_oneHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, 
                                   int part_init_ptr[], int* K_init_ptr, double eta_LR, double eta_CT, double eta_TD){
  // part_init_ptr is the partition at the LR level
  // at the moment Rests corresponds to the stable (CT)
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  rests = std::vector<std::vector<int> >(Rests); // this remembers the "stable" division of BG in CT
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }
  lowresPart = new Partition();
  std::vector<int> ids;
  for(int i = 0; i < nObsL; i++) {
    ids.push_back(i);
  }
  lowresPart->Initialize_Partition(nObsL, ids, eta_LR, part_init_ptr, K_init_ptr);
  
  highresPart = new Hdp();
  std::vector<std::vector<int> > CTinRest, Rest;
  CTinRest.resize( lowresPart->K ); // this is an equivalent of clusters for the LR partition
  Rest.resize( lowresPart->K ); // this is an equivalent of zLR_HR (the HR units which fall into the LR clusters) (another sort of clusters)
  int ct;
  for(int k = 0; k < lowresPart->K; k++){
    CTinRest[k].resize( lowresPart->cluster_config[k] );
    for(int h = 0; h < lowresPart->cluster_config[k]; h++){
        ct = lowresPart->clusters[k][h];
        CTinRest[k][h] = ct;
        for(int j = 0; j < Rests[ct].size(); j++){
          Rest[k].push_back(Rests[ct][j]); 
        }
    }
  }
  highresPart->Initialize_Hdp(nObs_highres, CTinRest, Rest, eta_CT, eta_TD);
  return;
}

void Multires::Initialize_Multires_partLR_nclHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, 
                                   int part_init_ptr[], int* K_init_ptr, double eta_LR, double eta_CT, double eta_TD){
  // part_init_ptr is the partition at the LR level
  // at the moment Rests corresponds to the stable (CT)
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  // this remembers the "stable" division of BG in CT
  rests = std::vector<std::vector<int> >(Rests); 
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }

  lowresPart = new Partition();
  std::vector<int> ids;
  for(int i = 0; i < nObsL; i++) {
    ids.push_back(i);
  }
  lowresPart->Initialize_Partition(nObsL, ids, eta_LR, part_init_ptr, K_init_ptr);
  
  highresPart = new Hdp();
  std::vector<std::vector<int> > CTinRest, Rest;
  CTinRest.resize( lowresPart->K ); // this is an equivalent of clusters for the LR partition
  Rest.resize( lowresPart->K ); // this is an equivalent of zLR_HR (the HR units which fall into the LR clusters) (another sort of clusters)
  int ct;
  for(int k = 0; k < lowresPart->K; k++){
    CTinRest[k].resize( lowresPart->cluster_config[k] );
    for(int h = 0; h < lowresPart->cluster_config[k]; h++){
      ct = lowresPart->clusters[k][h];
      CTinRest[k][h] = ct;
      for(int j = 0; j < Rests[ct].size(); j++){
        Rest[k].push_back(Rests[ct][j]); 
      }
    }
  }
  highresPart->Initialize_Hdp_nclusters(nObs_highres, CTinRest, Rest, eta_CT, eta_TD);
  return;
}

void Multires::Initialize_Multires_partLR_partHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, 
                                                int partLR_init_ptr[], int* KLR_init_ptr, int partHR_init_ptr[], double eta_LR, double eta_CT, double eta_TD){
  // Rests represents the highres areas in each lowres area
  // set lowres with one cluster
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  // this remembers the "stable" division of BG in CT
  rests = std::vector<std::vector<int> >(Rests);
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }

  lowresPart = new Partition();
  std::vector<int> ids;
  for(int i = 0; i < nObsL; i++) {
    ids.push_back(i);
  }
  lowresPart->Initialize_Partition(nObsL, ids, eta_LR, partLR_init_ptr, KLR_init_ptr);

  highresPart = new Hdp();
  std::vector<std::vector<int> > CTinRest, Rest;
  CTinRest.resize( lowresPart->K ); // this is an equivalent of clusters for the LR partition
  Rest.resize( lowresPart->K ); // this is an equivalent of zLR_HR (the HR units which fall into the LR clusters) (another sort of clusters)
  int ct;
  for(int k = 0; k < lowresPart->K; k++){
    CTinRest[k].resize( lowresPart->cluster_config[k] );
    for(int h = 0; h < lowresPart->cluster_config[k]; h++){
      ct = lowresPart->clusters[k][h];
      CTinRest[k][h] = ct;
      for(int j = 0; j < Rests[ct].size(); j++){
        Rest[k].push_back(Rests[ct][j]); 
      }
    }
  }
  // basically CTinRest tells me which CT are in each LR cluster (rest)
  // instead Rests tells me which BG belong to which LR cluster (like zLR_HR)
  highresPart->Initialize_Hdp_part_nrest(nObs_highres, CTinRest, Rest, partHR_init_ptr, eta_CT, eta_TD);
  return;
}


void Multires::Initialize_Multires_nclLR_nclHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, double eta_LR, double eta_CT, double eta_TD){
  nObs_lowres = nObsL;
  nObs_highres = nObsH;
  rests = std::vector<std::vector<int> >(Rests);
  rest_config = std::vector<int>(nObs_lowres);
  for(int r = 0; r < nObs_lowres; r++){
    rest_config[r] = Rests[r].size();
  }
  lowresPart = new Partition();
  highresPart = new Hdp();
  lowresPart->Initialize_Partition_nclusters(nObs_lowres, eta_LR); 
  // everyone is in a different cluster, so nObs_lowres restaurants
  std::vector<std::vector<int> > CTinRest;
  CTinRest.resize(nObsL);
  for(int i = 0; i < nObs_lowres; i++){
    CTinRest[i].resize(1);
    CTinRest[i][0] = i;
  }
  // highresPart->Initialize_Hdp_nRestDish(nObs_highres, CTinRest, Rests);
  highresPart->Initialize_Hdp_nclusters(nObs_highres, CTinRest, Rests, eta_CT, eta_TD);
  // every restaurant has one table and each table is in a different cluster
  return;
}

void Multires::Print_Multires(){
  cout << "Number of census tracts: " << nObs_lowres << endl;
  cout << "Number of block groups: " << nObs_highres << endl;
  cout << "Number of block groups per census tract: " << endl;
  for(int r = 0; r < nObs_lowres; r++){
    cout << rest_config[r] << " ";
  }
  cout << endl;
  cout << "*** Partition of census tracts: " << endl;
  lowresPart->Print_Partition();
  cout << "*** Partition of block groups: " << endl;
  highresPart->Print_Hdp();
  return;
}

void Multires::Get_Highres_CostDish(LPPartition highres_costdish){
  highresPart->get_CostDish(highres_costdish);
  return;
}

std::vector<bool> Multires::SampleMix(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD, bool sampleLR, bool sampleHR, bool local_print){
  // the first of changed always refers to lowres,
  // the second is true if any of the cost_table got accepted,
  // the last of changed always refers to tables_dishes.
  std::vector<bool> changed(1, false);
  if(sampleHR==false && alpha1==0.0){
    changed[0] = SampleSM_LR_noHR(ngibbs, temp, eng, eta_LR, eta_CT, eta_TD);
  } else {
    if(sampleLR){
      changed[0] = SampleNew_LR(temp, eng, eta_LR, eta_CT, eta_TD, local_print);
    }
    if(sampleHR){
      std::vector<bool> vhigh = highresPart->SampleSM(ngibbs, temp, eng, eta_CT, eta_TD, local_print);
      for(int i = 0; i < vhigh.size(); i++){
        changed.push_back(vhigh[i]);
      }
    } else {
      changed.push_back(false);
      changed.push_back(false);
    }
  }
  return changed;
}

std::vector<bool> Multires::SampleSM(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD, bool sampleLR, bool sampleHR, bool local_print){
  // the first of changed always refers to lowres,
  // the second is true if any of the cost_table got accepted,
  // the last of changed always refers to tables_dishes.
  std::vector<bool> changed(1, false);
  if(sampleHR==false && alpha1==0.0){
    changed[0] = SampleSM_LR_noHR(ngibbs, temp, eng, eta_LR, eta_CT, eta_TD);
  } else {
    if(sampleLR){
      changed[0] = SampleSM_LR(ngibbs, temp, eng, eta_LR, eta_CT, eta_TD, local_print);
      // changed[0] = SampleSM_LR(ngibbs, eng);
    }
    if(sampleHR){
      std::vector<bool> vhigh = highresPart->SampleSM(ngibbs, temp, eng, eta_CT, eta_TD,local_print);
      for(int i = 0; i < vhigh.size(); i++){
        changed.push_back(vhigh[i]);
      }
    } else {
      changed.push_back(false);
      changed.push_back(false);
    }
  }
  return changed;
}
// Initialize_Hdp gives nj[i] costumers to the (i+1)th restaurants
// But we want to say which groups of costumers are in which restaurants
// we could pass a vector of vector 
// Also: local ID (in each restaurant partition) is converted to mainID by
// assuming the order is fixed 

// Rcpp::List Multires::format_multires(){
//   Rcpp::List output_list;
  
//   Rcpp::List tmp_list;
//   tmp_list = lowresPart->format_partition();
//   output_list["LowRes"] = tmp_list;

//   LPPartition CostDish = new Partition();
//   Get_Highres_CostDish(CostDish);
//   tmp_list = CostDish->format_partition();
//   output_list["HighRes_CostDish"] = tmp_list;
//   delete CostDish;

//   return output_list;  
// }

bool Multires::SampleSM_LR(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD, bool local_print){
  bool returned_value;
  LPMultires proposal = new Multires(this);
  double lqratio, lqratio2, lpratio, llikratio;
  proposal->ProposalSM_LR(ngibbs, lqratio, lqratio2, lpratio, llikratio, eng, eta_LR, eta_CT, eta_TD, local_print);
  double a = lqratio + lqratio2 + temp * (lpratio + llikratio);
  if(local_print){
    cout << lqratio << " " << lqratio2 << " " << lpratio << " " << llikratio << " " << a << endl;
  }

  a = exp(a);
  a = min(1.0,a);
  uniform_real_distribution<double> distribution(0.0,1.0);
  if(distribution(eng) < a){
    Copy_Multires(proposal);
    returned_value = true;
    // cout << a << " accepted" << endl;
  } else {
    returned_value = false;
    // cout << a << endl;
  }
  delete proposal;
  return returned_value;
} 

bool Multires::SampleSM_LR_noHR(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD){
  bool returned_value;
  LPMultires proposal = new Multires(this);
  double lqratio, lqratio2, lpratio, llikratio;
  proposal->ProposalSM_LR(ngibbs, lqratio, lqratio2, lpratio, llikratio, eng, eta_LR, eta_CT, eta_TD);
  double a = lqratio + temp * (lpratio + llikratio);  // does not use lqratio2

  a = exp(a);
  a = min(1.0,a);
  uniform_real_distribution<double> distribution(0.0,1.0);
  if(distribution(eng) < a){
    Copy_Multires(proposal);
    returned_value = true;
    // cout << a << " accepted" << endl;
  } else {
    returned_value = false;
    // cout << a << endl;
  }
  delete proposal;
  return returned_value;
} 

void Multires::ProposalSM_LR(int ngibbs, double& lqratio, double& lqratio2, double& lpratio, double& llikratio, std::default_random_engine& myRandomEngine, double eta_LR, double eta_CT, double eta_TD, bool local_print){
  /*
  ProposalSM_LR first picks two LR units, checks if they are in the same cluster or not, then either splits or merges.
  - Split: first updates lowresPart, then highresPart (where we need to split a restaurant and the tables within that overlap with the splits)
    - Split_Rest (SplitTables_TableDish which updates tables_dishes + Divide_Restaurants which updates costumers_tables)
    - get_tables_merge_reverse: computes the probability of merging the tables that have been split back into the merged configuration before splitting
  - Merge: first updates lowresPart, then highresPart (where we need to merge two restaurants and the tables within)
    - Merge_Rest: puts the tables in one unique restaurant but does not merge them
      (Combine_Restaurants changes costumers_tables + MoveTables_TableDish which changes tables_dishes)
    - get_tables_merge: figures out which pair of tables should be merged and the probability
      (Merge_CostTable operates on the pair of tables with Merge first and then RemoveTable_TableDish for tables_dishes)
  */
  LPPartition startingP = lowresPart;
  std::vector<int> i12 = sample_without_replacement(2, startingP->nObs, myRandomEngine);
  int i1 = i12[0];
  int i2 = i12[1];
  int c1 = startingP->cluster_assignment[i1];
  int c2 = startingP->cluster_assignment[i2];
  if(c2 < c1){
    // c1 should be smaller than c2
    // change the labels i1 and i2 so that c1 < c2
    int tmp;
    tmp = i1;
    i1 = i2;
    i2 = tmp;
    tmp = c1;
    c1 = c2;
    c2 = tmp;
  }
  // samecl_index contains the other people that belong to c1 or c2 (but not i1, i2) localID
  // if c1 == c2, then it's the other people in that cluster
  // if c1 != c2, then it's the other people in the two clusters
  std::vector<int> samecl_index = startingP->get_samecl_index(i1, i2, c1, c2);

  // launch will NOT modify the current partition, but update new_c1_ind, new_c2_ind, samecl_newcl
  std::list<int> new_c1_ind,new_c2_ind; // lists are optimized to remove elements in any position (use "remove")
  std::vector<int> samecl_newcl;
  int params[4] = {i1, i2, c1, c2};
  
  LaunchSM_LR(params, samecl_index, new_c1_ind, new_c2_ind, samecl_newcl, ngibbs, myRandomEngine, eta_TD);
  if(c1 == c2){
    lqratio = 0;
    // lqratio2 = 0;
    ProposalSM_LR_iter(true, params, samecl_index, new_c1_ind, new_c2_ind, samecl_newcl, lqratio, myRandomEngine, eta_CT);
    int nc1 = new_c1_ind.size();
    int nc2 = new_c2_ind.size();
    lpratio = log(eta_LR) + lbeta(nc1,nc2);

    llikratio = alpha1 * (- highresPart->costumers_tables[c1]->get_logprior(eta_CT)); 
    if(local_print){
      cout << (- highresPart->tables_dishes->get_logprior(eta_TD)) << endl;
    }
    llikratio += alpha1 * (- highresPart->tables_dishes->get_logprior(eta_TD)); 
    if(yLR_bool) llikratio += alpha2 * (- lowresPart->cluster_llikelihood(c1, true));
    // NOW WE UPDATE THE CURRENT PARTITION, BY SPLITTING
    // index1, index2, index1v, index2v contain the indices of the groups that get split in one or the other
    new_c1_ind.sort();
    new_c2_ind.sort();
    int *index1, *index2;
    std::vector<int> index1v(nc1);
    std::vector<int> index2v(nc2);
    index1 = new int[nc1];
    index2 = new int[nc2];
    if(local_print){
      cout << "LR Split: index1: ";
    }
    for(int i = 0; i < nc1; i++){
      index1[i] = new_c1_ind.front();
      index1v[i] = index1[i];
      new_c1_ind.pop_front();
      if(local_print){
        cout << index1[i] << " ";
      }
    }
    if(local_print){
      cout << " and index2: ";
    }
    for(int i = 0; i < nc2; i++){
      index2[i] = new_c2_ind.front();
      index2v[i] = index2[i];
      new_c2_ind.pop_front();
      if(local_print){
        cout << index2[i] << " ";
      }
    }
    if(local_print){
      cout << endl;
    }
    lowresPart->Split(c1, index1, index2, nc1, nc2, eta_LR);
    delete[] index1;
    delete[] index2;

    LPPartition oldmerged = new Partition(highresPart->costumers_tables[c1]);
    // now we need to split the restaurants, and the tables are divided like we devised.
    if(local_print){
      cout << "before split_rest" << endl;
    }
    highresPart->Split_Rest(c1, lowresPart->K-1, index1v, index2v, eta_CT, eta_TD, local_print);
    if(local_print){
      cout << "after split_rest" << endl;
    }
    // we also need to compute the probability of the opposite move (merging the tables)
    lqratio2 = get_tables_merge_reverse(c1, lowresPart->K-1, oldmerged, myRandomEngine);
    delete oldmerged;

    llikratio += alpha1 * (highresPart->costumers_tables[c1]->get_logprior(eta_CT) + highresPart->costumers_tables[lowresPart->K-1]->get_logprior(eta_CT)); 
    if(local_print){
      cout << (highresPart->tables_dishes->get_logprior(eta_TD)) << endl;
    }
    llikratio += alpha1 * (highresPart->tables_dishes->get_logprior(eta_TD));  
    if(yLR_bool) llikratio += alpha2 * (lowresPart->cluster_llikelihood(c1, true) + lowresPart->cluster_llikelihood(lowresPart->K-1, true));
  } else {
    lqratio = 0;
    lpratio = - log(eta_LR) - lbeta(lowresPart->cluster_config[c1], lowresPart->cluster_config[c2]);
    ProposalSM_LR_iter(false, params, samecl_index, new_c1_ind, new_c2_ind, samecl_newcl, lqratio, myRandomEngine, eta_CT);

    llikratio = alpha1 * (- highresPart->costumers_tables[c1]->get_logprior(eta_CT) - highresPart->costumers_tables[c2]->get_logprior(eta_CT)); 
    llikratio += alpha1 * (- highresPart->tables_dishes->get_logprior(eta_TD)); 
    if(yLR_bool) llikratio += alpha2 * (- lowresPart->cluster_llikelihood(c1, true) - lowresPart->cluster_llikelihood(c2, true));

    int threshold_K1old = highresPart->costumers_tables[c1]->K; // This will be useful in get_tables_merge
    lowresPart->Merge(c1,c2, eta_LR);
    highresPart->Merge_Rest(c1,c2); // we need to merge the restaurants 
    lqratio2 = -get_tables_merge(c1, threshold_K1old, myRandomEngine, eta_CT, eta_TD); // and the tables
    // get_tables_merge computes the probability of merging (not all of them are merged) BUT also merges! (calls Hdp::Merge_CostTable)

    // cout << "-- new partition --" << endl;
    // highresPart->costumers_tables[c1]->Print_Partition_Short();

    llikratio += alpha1 * (+ highresPart->costumers_tables[c1]->get_logprior(eta_CT)); 
    llikratio += alpha1 * (+ highresPart->tables_dishes->get_logprior(eta_TD)); 
    if(yLR_bool) llikratio += alpha2 * (+ lowresPart->cluster_llikelihood(c1, true));
  }
  // if(print_bool){
  //   if((i1 == 2 && i2 == 3) || (i1 == 3 && i2 == 2)){
  //     cout << lqratio << " " << lqratio2 << " " << lpratio << " " << llikratio << endl;
  //   }
  // }
  return;
}

void Multires::LaunchSM_LR(int params[], std::vector<int> samecl_index, std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, 
  std::vector<int>& samecl_newcl, int ngibbs, std::default_random_engine& gen, double eta_CT)
{
  // will NOT modify the current partition
  int i1, i2, c1, c2, new_c1, new_c2;
  i1 = params[0];
  i2 = params[1];
  c1 = params[2];
  c2 = params[3];
  if(c1 == c2){
    new_c1 = c1;
    new_c2 = lowresPart->K;
  } else {
    new_c1 = c1;
    new_c2 = c2;
  }
  // samecl_index contain the indices (from 1 to n) of people belonging to c1 or c2
  // new_c1_ind, new_c2_ind will contain the indices (from 1 to n - localID) divided into the two new clusters
  // samecl_newcl contains either new_c1 or new_c2 for the elements of samecl_index, with the same order of samecl_index

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
      } else {
        new_c2_ind.push_back(samecl_index[i]); // new_c2
        samecl_newcl[i] = new_c2;
      }
    }
    double lqratio_ignored = 0;
    for(int t = 0; t < ngibbs; t++){
      ProposalSM_LR_iter(true, params, samecl_index, new_c1_ind, new_c2_ind, samecl_newcl, lqratio_ignored, gen, eta_CT);
    }
  } 
  return;
}

void Multires::ProposalSM_LR_iter(bool SplitMerge, int params[], std::vector<int> samecl_index,
  std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, 
  double& lqratio, std::default_random_engine& myRandomEngine, double eta_CT)
{
  // SplitMerge = true: means we are splitting; false: we are merging
  int i1, i2, c1, c2, new_c1, new_c2;
  i1 = params[0];
  i2 = params[1];
  c1 = params[2];
  c2 = params[3];
  if(c1 == c2){
    new_c1 = c1;
    new_c2 = lowresPart->K;
  } else {
    new_c1 = c1;
    new_c2 = c2;
  }

  arma::vec P(2);
  arma::vec l(2);
  arma::vec logpr, pr;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  if(SplitMerge){
    for(int i = 0; i < samecl_index.size(); i++){
      if(samecl_newcl[i] == new_c1){
        P[0] = new_c1_ind.size()-1;
        P[1] = new_c2_ind.size();
        if(yLR_bool) {
          l[0] = alpha2 * (lowresPart->cluster_llikelihood(new_c1_ind, true) + lowresPart->cluster_llikelihood(new_c2_ind, true));
          l[0] += alpha1 * splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT);
        } else {
          l[0] = splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT); 
        }
        new_c1_ind.remove(samecl_index[i]);
        new_c2_ind.push_back(samecl_index[i]);
        if(yLR_bool) {
          l[1] = alpha2 * (lowresPart->cluster_llikelihood(new_c1_ind, true) + lowresPart->cluster_llikelihood(new_c2_ind, true));
          l[1] += alpha1 * splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT);
        } else {
          l[1] = splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT); 
        }
        new_c2_ind.pop_back();
        new_c1_ind.push_back(samecl_index[i]);
      } else {
        P[0] = new_c1_ind.size();
        P[1] = new_c2_ind.size()-1;
        if(yLR_bool) {
          l[1] = alpha2 * (lowresPart->cluster_llikelihood(new_c1_ind, true) + lowresPart->cluster_llikelihood(new_c2_ind, true));
          l[1] += alpha1 * splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT);
        } else {
          l[1] = splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT); 
        }
        new_c2_ind.remove(samecl_index[i]);
        new_c1_ind.push_back(samecl_index[i]);
        if(yLR_bool) {
          l[0] = alpha2 * (lowresPart->cluster_llikelihood(new_c1_ind, true) + lowresPart->cluster_llikelihood(new_c2_ind, true));
          l[0] += alpha1 * splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT);
        } else {
          l[0] = splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT); 
        }
        new_c1_ind.pop_back();
        new_c2_ind.push_back(samecl_index[i]);
      }
      logpr = log(P)+l;
      logpr = logpr - max(logpr);
      logpr = logpr - log(sum(exp(logpr)));
      pr = exp(logpr);
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
      }
    }
  } else {
    for(int i = 0; i < samecl_index.size(); i++){
      if(samecl_newcl[i] == new_c1){
        P[0] = new_c1_ind.size()-1;
        P[1] = new_c2_ind.size();
        if(yLR_bool) {
          l[0] = alpha2 * (lowresPart->cluster_llikelihood(new_c1_ind, true) + lowresPart->cluster_llikelihood(new_c2_ind, true));
          l[0] += alpha1 * splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT);
        } else {
          l[0] = splitted_tables(c1, c2, new_c1_ind, new_c2_ind, eta_CT);
        }
        new_c1_ind.remove(samecl_index[i]);
        new_c2_ind.push_back(samecl_index[i]);
        if(yLR_bool) {
          l[1] = alpha2 * (lowresPart->cluster_llikelihood(new_c1_ind, true) + lowresPart->cluster_llikelihood(new_c2_ind, true));
          l[1] += alpha1 * splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT);
        } else {
          l[1] = splitted_tables(c1, c2, new_c1_ind, new_c2_ind, eta_CT);
        }
        new_c2_ind.pop_back();
        new_c1_ind.push_back(samecl_index[i]);
      } else {
        P[0] = new_c1_ind.size();
        P[1] = new_c2_ind.size()-1;
        if(yLR_bool) {
          l[1] = alpha2 * (lowresPart->cluster_llikelihood(new_c1_ind, true) + lowresPart->cluster_llikelihood(new_c2_ind, true));
          l[1] += alpha1 * splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT);
        } else {
          l[1] = splitted_tables(c1, c2, new_c1_ind, new_c2_ind, eta_CT);
        }
        new_c2_ind.remove(samecl_index[i]);
        new_c1_ind.push_back(samecl_index[i]);
        if(yLR_bool) {
          l[0] = alpha2 * (lowresPart->cluster_llikelihood(new_c1_ind, true) + lowresPart->cluster_llikelihood(new_c2_ind, true));
          l[0] += alpha1 * splitted_tables(c1, new_c1_ind, new_c2_ind, eta_CT);
        } else {
          l[0] = splitted_tables(c1, c2, new_c1_ind, new_c2_ind, eta_CT);
        }
        new_c1_ind.pop_back();
        new_c2_ind.push_back(samecl_index[i]);
      }
      logpr = log(P)+l;
      logpr = logpr - max(logpr);
      logpr = logpr - log(sum(exp(logpr)));
      pr = exp(logpr);
      if(lowresPart->cluster_assignment[ samecl_index[i] ] == new_c1){
        lqratio += logpr[0];
      } else if(lowresPart->cluster_assignment[ samecl_index[i] ] == new_c2) {
        lqratio += logpr[1];
      } else {
        cout << "ProposalSM_LR_iter Problem!, eta_CT";
      }
    }
  }
  return;
}

bool Multires::SampleNew_LR(double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD, bool local_print){
  bool returned_value;
  
  double lqratio, lpratio, llikratio;
  proposal_info info;
  ProposalNew_LR(lqratio, lpratio, llikratio, &info, eng, eta_LR, eta_CT, eta_TD, local_print);
  double a = lqratio + temp * (lpratio + llikratio);
  if(local_print){
    cout << lqratio << " " << lpratio << " " << llikratio << " " << a << endl;
  }

  a = exp(a);
  a = min(1.0,a);
  uniform_real_distribution<double> distribution(0.0,1.0);
  if(distribution(eng) < a){
    // here we modify the partition using info.
    save_proposal(&info, eta_LR, eta_CT, eta_TD, local_print);
    // if(local_print) Print_Multires();
    // if(local_print) highresPart->tables_dishes->Print_Partition();
    returned_value = true;
    if(local_print) cout << "ACCEPTED" << endl;
  } else {
    returned_value = false;
  }
  return returned_value;
} 

void Multires::ProposalNew_LR(double& lqratio, double& lpratio, double& llikratio, proposal_info* info, std::default_random_engine& myRandomEngine, double eta_LR, double eta_CT, double eta_TD, bool local_print){
  /*
  ProposalNew_LR first picks one LR unit, samples the new cluster assignment.
  */
  double lpr_q = 0.0;         // log probability of the proposed change
  double lpr_q_inv = 0.0;     // log probability of the inverse move
  
  // Now we sample the LR unit (j) to change
  std::uniform_int_distribution<int> distribution_unif_int(0, lowresPart->nObs-1);
  int j = distribution_unif_int(myRandomEngine);
  info->jLR = j;
  info->jLR_index = lowresPart->indexing[j];

  // Now we sample the new proposed cluster assignment of j
  std::vector<double> pr(lowresPart->K + 1);
  for(int k = 0; k < lowresPart->K; k++){ 
    pr[k] = lowresPart->cluster_config[k]/(lowresPart->nObs-1 + eta_LR);
  }
  int kLR_old = lowresPart->cluster_assignment[j];
  pr[kLR_old] = pr[kLR_old] - 1/(lowresPart->nObs-1 + eta_LR); // we ignore the current assignment of j
  pr[lowresPart->K] = eta_LR/(lowresPart->nObs-1 + eta_LR);
  std::discrete_distribution<int> distribution(pr.begin(), pr.end());
  int kLR_new = distribution(myRandomEngine);

  info->kLR_old = kLR_old;
  info->kLR_new = kLR_new;

  if(local_print){
    cout << j << " " << info->jLR_index << " " << kLR_old << " " << kLR_new << endl;
  }

  lpr_q += log(pr[kLR_new]);
  if(pr[kLR_old] > 0){
    lpr_q_inv += log(pr[kLR_old]);
  } else {
    lpr_q_inv += log(eta_LR/(lowresPart->nObs-1 + eta_LR));
  }
  lpratio = lpr_q - lpr_q_inv;

  if(local_print){
    cout << lpr_q << " " << lpr_q_inv << endl;
  }
  
  // Info on the HR units (costumers) within the LR unit j:
  std::vector<int> cost_locID_j;
  std::vector<int> cost_locID_NOTj;
  std::vector<int> cost_indexing_j;
  std::vector<int> oldtables;
  std::vector<int> oldtables_mainID;
  std::vector<int> olddishes;
  int temp, temp2;
  for(int i = 0; i < highresPart->costumers_tables[kLR_old]->nObs; i++){
    if(InvMapping(highresPart->costumers_tables[kLR_old]->indexing[i],1) == info->jLR_index){
      cost_locID_j.push_back(i);
      cost_indexing_j.push_back(highresPart->costumers_tables[kLR_old]->indexing[i]);
      temp = highresPart->costumers_tables[kLR_old]->cluster_assignment[i];
      oldtables.push_back(temp);
      temp2 = Table_RestId_MainId(kLR_old, temp, highresPart->nTableRest);
      oldtables_mainID.push_back(temp2);
      olddishes.push_back(highresPart->tables_dishes->cluster_assignment[temp2]);
    } else {
      cost_locID_NOTj.push_back(i);
    }
  }
  int n_cost_j = cost_locID_j.size();

  if(local_print){
    cout << "cost_locID_j: ";
    for(int i = 0; i < n_cost_j; i++)
      cout << cost_locID_j[i] << " ";
    cout << endl;
    cout << "cost_indexing_j: ";
    for(int i = 0; i < n_cost_j; i++)
      cout << cost_indexing_j[i] << " ";
    cout << endl;
    cout << "oldtables: ";
    for(int i = 0; i < n_cost_j; i++)
      cout << oldtables[i] << " ";
    cout << endl;
    cout << "oldtables_mainID: ";
    for(int i = 0; i < n_cost_j; i++)
      cout << oldtables_mainID[i] << " ";
    cout << endl;
    cout << "olddishes: ";
    for(int i = 0; i < n_cost_j; i++)
      cout << olddishes[i] << " ";
    cout << endl;
  }
  
  // Info on the tables and dishes (minus the HR units in j):
  std::vector<int> remove_table(highresPart->costumers_tables[kLR_old]->K,0); // all initialized with 0
  std::vector<int> remove_dish(highresPart->tables_dishes->K,0);            // all initialized with 0
  int n_tables = highresPart->tables_dishes->nObs;  // number (over all restaurants) after removing the costumers of j
  int K_dishes = highresPart->tables_dishes->K;     // K_dishes is the size of nk_dishes (some elements of nk_dishes can be zero)
  std::vector<int> nk_dishes(K_dishes,0);           // number of tables after removing the costumers of j (can be 0)
  std::vector<double> pr_dishes(K_dishes+1);
  std::vector<double> dishes_thetas_copy(K_dishes); // copy (which will get modified) of highrestPart->dishes_thetas
  int temp_table, temp_rest;
  for(int k = 0; k < highresPart->costumers_tables[kLR_old]->K; k++){
    temp = 0; // how many from cost_locID_j were seated at table k
    for(int i = 0; i < oldtables.size(); i++){
      if(oldtables[i] == k){
        temp++;
      }
    }
    if(temp == highresPart->costumers_tables[kLR_old]->cluster_config[k]){
      remove_table[k] = 1;
      n_tables += -1;
      if(local_print) cout << "remove_table " << k << endl;
    }
  }
  for(int d = 0; d < highresPart->tables_dishes->K; d++){
    // if(local_print) cout << "dish " << d << endl;
    dishes_thetas_copy[d] = highresPart->dishes_thetas[d];

    for(int it = 0; it < highresPart->tables_dishes->cluster_config[d]; it++){
      temp = highresPart->tables_dishes->clusters[d][it];
      Table_MainId_RestId(temp, highresPart->nTableRest, temp_rest, temp_table);
      // if(local_print) cout << "table " << temp << " " << temp_table << " " << temp_rest << endl;
      if(temp_rest != kLR_old){
        nk_dishes[d]++;
      } else {
        if(remove_table[temp_table] != 1){
          nk_dishes[d]++;
        }
      }
    }
    if(nk_dishes[d] == 0){
      remove_dish[d] = 1;
      // K_dishes += -1;
      if(local_print) cout << "remove_dish " << d << endl;
    } 
  }

  // Here we will keep track of the prosal for table and dish assignment of the HR units within j
  std::vector<int> new_tables(n_cost_j);        // contains the new table assignment (localID) for the corresponding element of cost_locID_j
  std::vector<int> new_tables_mainID(n_cost_j); // contains the new table assignment (mainID) for the corresponding element of cost_locID_j
  std::vector<int> new_dishes(n_cost_j);        // contains the new dish assignment of the table of the corresponding element of cost_locID_j

  // let's prepare the probabilities for sampling the new proposal
  int n_cost;                 // number of HR units in the new restaurant after removing elements in j (gets updated within the loop)
  int K_tables;               // length of the vector nk_tables (including any 0 element)
  int K_tables_mr;            // length of the vector nk_tables "minus removed" (mr), ie we count if a table is empty after removing elements in j
  int K_tables_orig;          // number of tables before starting the resampling process. all the tables with new_tables[i]>=K_tables_orig are new tables
  std::vector<int> nk_tables; // vector with number of costumers seated at each table (after removing elements in j)
  std::vector<double> pr_tables; // probabilities used to sample the new table assignment
  std::vector<double> tables_dishes_thetas; // for each table in nk_tables, which theta is assigned to the dish serving there?
  
  // we'll use the struct to avoid having a thousand variables with a thousand names.
  struct info_sampling {
    int n_cost;                 // number of HR units in the new restaurant after removing elements in j (gets updated within the loop)
    int K_tables;               // length of the vector nk_tables (including any 0 element)
    int K_tables_mr;            // length of the vector nk_tables "minus removed" (mr), ie we count if a table is empty after removing elements in j
    int K_tables_orig;          // number of tables before starting the resampling process. all the tables with new_tables[i]>=K_tables_orig are new tables
    std::vector<int> nk_tables; // vector with number of costumers seated at each table (after removing elements in j)
    std::vector<double> pr_tables; // probabilities used to sample the new table assignment
    int n_tables;  // number (over all restaurants) after removing the costumers of j
    int K_dishes;     // K_dishes is the size of nk_dishes (some elements of nk_dishes can be zero)
    std::vector<int> nk_dishes;           // number of tables after removing the costumers of j (can be 0)
    std::vector<double> pr_dishes;
  } Pinfo, Pinfo_old; // Qinfo, Qinfo_old

  // let's compute the probabilities for computing the reverse move
  int n_cost_old = highresPart->costumers_tables[kLR_old]->nObs - n_cost_j;                 
  int K_tables_old = highresPart->costumers_tables[kLR_old]->K; 
  int K_tables_mr_old = K_tables_old;  
  for(int k = 0; k < highresPart->costumers_tables[kLR_old]->K; k++) // don't know if useful
    if(remove_table[k] == 1)
      K_tables_mr_old += -1;  
  int K_dishes_old = K_dishes; 
  int n_tables_old = n_tables;         
  std::vector<int> nk_tables_old(K_tables_old); 
  std::vector<double> pr_tables_old(K_tables_old+1);
  std::vector<int> nk_dishes_old(K_dishes_old); 
  std::vector<double> pr_dishes_old(K_dishes_old+1);
  std::vector<double> tables_dishes_thetas_old(K_tables_old);
  for(int k = 0; k < K_tables_old; k++)
    nk_tables_old[k] = highresPart->costumers_tables[kLR_old]->cluster_config[k];
  for(int i = 0; i < n_cost_j; i++)
    nk_tables_old[oldtables[i]] += -1;
  for(int k = 0; k < K_tables_old; k++){
    pr_tables_old[k] = nk_tables_old[k]/(n_cost_old+eta_CT);
    temp = Table_RestId_MainId(kLR_old, k, highresPart->nTableRest);
    tables_dishes_thetas_old[k] = highresPart->dishes_thetas[ highresPart->tables_dishes->cluster_assignment[temp] ];
  }
  pr_tables_old[K_tables_old] = eta_CT/(n_cost_old+eta_CT);

  nk_dishes_old.assign(nk_dishes.begin(), nk_dishes.end());
  
  Pinfo_old.n_cost = n_cost_old;
  Pinfo_old.K_tables = K_tables_old;
  Pinfo_old.K_tables_mr = K_tables_mr_old;
  Pinfo_old.n_tables = n_tables_old;
  Pinfo_old.K_dishes = K_dishes_old;
  Pinfo_old.nk_tables = nk_tables_old;
  Pinfo_old.pr_tables = pr_tables_old;
  Pinfo_old.nk_dishes = nk_dishes_old;
  Pinfo_old.pr_dishes = pr_dishes_old;
  // Qinfo_old = Pinfo_old;
  
  if(local_print){
    cout << "n_cost_old: " << n_cost_old << endl;
    cout << "nk_tables_old: ";
    for(int i = 0; i < K_tables_old; i++)
      cout << nk_tables_old[i] << " ";
    cout << endl;
  }

  if(kLR_new == lowresPart->K){
    // here nk_tables = 0
    n_cost = 0;
    K_tables = 0;
    K_tables_orig = 0;
    K_tables_mr = 0;
    pr_tables.resize(1);
    pr_tables[K_tables] = 1;    
  } else {
    if(kLR_new != kLR_old){
      n_cost = highresPart->costumers_tables[kLR_new]->nObs;
      K_tables = highresPart->costumers_tables[kLR_new]->K;
      K_tables_orig = K_tables;
      K_tables_mr = K_tables;
      nk_tables.resize(K_tables);
      pr_tables.resize(K_tables+1);
      tables_dishes_thetas.resize(K_tables);
      for(int k = 0; k < K_tables; k++){ // NOTE: k is a localID
        nk_tables[k] = highresPart->costumers_tables[kLR_new]->cluster_config[k];
        pr_tables[k] = nk_tables[k]/(n_cost+eta_CT);
        temp = Table_RestId_MainId(kLR_new, k, highresPart->nTableRest);
        tables_dishes_thetas[k] = highresPart->dishes_thetas[ highresPart->tables_dishes->cluster_assignment[temp] ];
      }
      pr_tables[K_tables] = eta_CT/(n_cost+eta_CT);      
    } else {
      n_cost = n_cost_old;
      K_tables = K_tables_old;
      K_tables_orig = K_tables;
      K_tables_mr = K_tables_mr_old;

      nk_tables.assign(nk_tables_old.begin(), nk_tables_old.end());
      // std::vector<int>::iterator it;
      // it=nk_tables_old.begin()+1;
      // nk_tables.assign (it,nk_tables_old.end()-1);
      pr_tables.assign(pr_tables_old.begin(), pr_tables_old.end());
      tables_dishes_thetas.assign(tables_dishes_thetas_old.begin(), tables_dishes_thetas_old.end());
    }
  }

  Pinfo.n_cost = n_cost;
  Pinfo.K_tables = K_tables;
  Pinfo.K_tables_mr = K_tables_mr;
  Pinfo.K_tables_orig = K_tables_orig;
  Pinfo.n_tables = n_tables;
  Pinfo.K_dishes = K_dishes;
  Pinfo.nk_tables = nk_tables;
  Pinfo.pr_tables = pr_tables;
  Pinfo.nk_dishes = nk_dishes;
  Pinfo.pr_dishes = pr_dishes;
  // Qinfo = Pinfo;

  if(local_print){
    cout << "n_cost: " << n_cost << endl;
    cout << "pr_tables: ";
    for(int i = 0; i < K_tables+1; i++)
      cout << pr_tables[i] << " ";
    cout << endl;
  }

  // Now we sample the table (and dish) assignment:
  int i, new_table_i, new_table_i_main, new_dish_i;
  double temp_loglik, beta_n, temp_sum;
  for(int it = 0; it < n_cost_j; it++){
    i = cost_locID_j[it];

    // update pr_tables with the "likelihood" part
    // copied and adapted from dish_llikelihood
    temp_sum = 0.0;
    for(int k = 0; k < K_tables; k++){
      temp_loglik = lgamma(alpha_H+0.5)-lgamma(alpha_H);
      temp_loglik += alpha_H*log(beta_H) - (alpha_H+0.5)*log(beta_H + pow(Y(cost_indexing_j[it])-tables_dishes_thetas[k],2));
      temp_loglik += -0.5*log(M_PI);
      pr_tables[k] *= exp(temp_loglik);
      temp_sum += pr_tables[k];
    }
    beta_n = beta_H + k_H * pow(Y(cost_indexing_j[it])-mu_H,2)/(2*(k_H + 1));
    temp_loglik = lgamma(alpha_H+0.5) - lgamma(alpha_H) + alpha_H * log(beta_H) - (alpha_H+0.5) * log(beta_n);
    temp_loglik += 0.5*(log(k_H)-log(k_H+1)) -1/2*log(2*M_PI);
    pr_tables[K_tables] *= exp(temp_loglik);
    temp_sum += pr_tables[K_tables];
    // normalize
    for(int k = 0; k < K_tables; k++){
      pr_tables[k] /= temp_sum;
    }
    pr_tables[K_tables] /= temp_sum;

    std::discrete_distribution<int> distribution(pr_tables.begin(), pr_tables.end());
    new_table_i = distribution(myRandomEngine);
    new_table_i_main = Table_RestId_MainId(kLR_new, new_table_i, highresPart->nTableRest);
    new_tables[it] = new_table_i;
    new_tables_mainID[it] = new_table_i_main;
    
    lpr_q += log(pr_tables[new_table_i]);

    if(new_table_i < K_tables){
      nk_tables[new_table_i]++;
      if((kLR_new < lowresPart->K) && (new_table_i < highresPart->costumers_tables[kLR_new]->K)){
        new_dishes[it] = highresPart->tables_dishes->cluster_assignment[new_table_i_main];
      } else { // table added in the new proposal
        for(int itt = 0; itt < it; itt++){
          if(new_tables[itt] == new_table_i){
            new_dishes[it] = new_dishes[itt];
            break;
          }
        }
      }
    } else {
      nk_tables.push_back(1);
      pr_tables.push_back(1);
      K_tables++;
      // Now we choose a new dish (and compute probabilities)
      temp_sum = 0.0;
      for(int d = 0; d < K_dishes; d++){
        pr_dishes[d] = nk_dishes[d]/(n_tables + eta_TD);
        // add the loglikelihood part
        temp_loglik = lgamma(alpha_H+0.5)-lgamma(alpha_H);
        temp_loglik += alpha_H*log(beta_H) - (alpha_H+0.5)*log(beta_H + pow(Y(cost_indexing_j[it])-dishes_thetas_copy[d],2));
        temp_loglik += -0.5*log(M_PI);
        pr_dishes[d] *= exp(temp_loglik);
        temp_sum += pr_dishes[d];
      }
      pr_dishes[K_dishes] = eta_TD/(n_tables + eta_TD);
      // add the loglikelihood part
      beta_n = beta_H + k_H * pow(Y(cost_indexing_j[it])-mu_H,2)/(2*(k_H + 1));
      temp_loglik = lgamma(alpha_H+0.5) - lgamma(alpha_H) + alpha_H * log(beta_H) - (alpha_H+0.5) * log(beta_n);
      temp_loglik += 0.5*(log(k_H)-log(k_H+1)) -1/2*log(2*M_PI);
      pr_dishes[K_dishes] *= exp(temp_loglik);
      temp_sum += pr_dishes[K_dishes];
      // normalize
      for(int d = 0; d < K_dishes; d++){
        pr_dishes[d] /= temp_sum;
      }
      pr_dishes[K_dishes] /= temp_sum;

      std::discrete_distribution<int> distribution(pr_dishes.begin(), pr_dishes.end());
      new_dish_i = distribution(myRandomEngine);
      new_dishes[it] = new_dish_i;
      lpr_q += log(pr_dishes[new_dish_i]);
      
      if(new_dish_i < K_dishes){
        nk_dishes[new_dish_i]++;
        tables_dishes_thetas.push_back( dishes_thetas_copy[new_dish_i] ); 
      } else {
        nk_dishes.push_back(1);
        pr_dishes.push_back(1);
        tables_dishes_thetas.push_back( Y(cost_indexing_j[it]) );
        dishes_thetas_copy.push_back( Y(cost_indexing_j[it]) );
        K_dishes++;
      }
      n_tables++;
    }
    // if(local_print) cout << new_dishes[it] << endl;
    // Now we compute the probabilities for sampling the next guy.
    n_cost++;
    for(int k = 0; k < K_tables; k++){
      pr_tables[k] = nk_tables[k]/(n_cost+eta_CT);
    }
    pr_tables[K_tables] = eta_CT/(n_cost+eta_CT);
  }
  // Now we compute the probability of the reverse move. 
  for(int it = 0; it < n_cost_j; it++){
    i = cost_locID_j[it];

    temp_sum = 0.0;
    for(int k = 0; k < K_tables_old; k++){
      temp_loglik = lgamma(alpha_H+0.5)-lgamma(alpha_H);
      temp_loglik += alpha_H*log(beta_H) - (alpha_H+0.5)*log(beta_H + pow(Y(cost_indexing_j[it])-tables_dishes_thetas_old[k],2));
      temp_loglik += -0.5*log(M_PI);
      pr_tables_old[k] *= exp(temp_loglik);
      temp_sum += pr_tables_old[k];
    }
    beta_n = beta_H + k_H * pow(Y(cost_indexing_j[it])-mu_H,2)/(2*(k_H + 1));
    temp_loglik = lgamma(alpha_H+0.5) - lgamma(alpha_H) + alpha_H * log(beta_H) - (alpha_H+0.5) * log(beta_n);
    temp_loglik += 0.5*(log(k_H)-log(k_H+1)) -1/2*log(2*M_PI);
    pr_tables_old[K_tables_old] *= exp(temp_loglik);
    temp_sum += pr_tables_old[K_tables_old];
    for(int k = 0; k < K_tables_old; k++){
      pr_tables_old[k] /= temp_sum;
    }
    pr_tables_old[K_tables_old] /= temp_sum;

    if(pr_tables_old[oldtables[it]] > 0){
      lpr_q_inv += log(pr_tables_old[oldtables[it]]);
    } else {
      // if i is the first one to seat at its table, we pretend it's a new one, 
      // and we compute the probability of sampling the right dish
      lpr_q_inv += log(pr_tables_old[K_tables_old]);
      for(int d = 0; d < K_dishes_old; d++){
        pr_dishes_old[d] = nk_dishes_old[d]/(n_tables_old + eta_TD);

        temp_loglik = lgamma(alpha_H+0.5)-lgamma(alpha_H);
        temp_loglik += alpha_H*log(beta_H) - (alpha_H+0.5)*log(beta_H + pow(Y(cost_indexing_j[it])-dishes_thetas_copy[d],2));
        temp_loglik += -0.5*log(M_PI);
        pr_dishes_old[d] *= exp(temp_loglik);
        temp_sum += pr_dishes_old[d];
      }
      pr_dishes_old[K_dishes_old] = eta_TD/(n_tables_old + eta_TD);
      beta_n = beta_H + k_H * pow(Y(cost_indexing_j[it])-mu_H,2)/(2*(k_H + 1));
      temp_loglik = lgamma(alpha_H+0.5) - lgamma(alpha_H) + alpha_H * log(beta_H) - (alpha_H+0.5) * log(beta_n);
      temp_loglik += 0.5*(log(k_H)-log(k_H+1)) -1/2*log(2*M_PI);
      pr_dishes_old[K_dishes_old] *= exp(temp_loglik);
      temp_sum += pr_dishes_old[K_dishes_old];
      for(int d = 0; d < K_dishes_old; d++){
        pr_dishes_old[d] /= temp_sum;
      }
      // pr_dishes_old[K_dishes_old] /= temp_sum;

      if(pr_dishes_old[olddishes[it]] > 0){
        lpr_q_inv += log(pr_dishes_old[olddishes[it]]);
      } else {
        lpr_q_inv += log( eta_TD/(n_tables_old + eta_TD) ); // == pr_dishes_old[K_dishes_old]
      }
      nk_dishes_old[olddishes[it]]++;
      n_tables_old++;
    }
    nk_tables_old[oldtables[it]]++;
    n_cost_old++;
    for(int k = 0; k < K_tables_old; k++)
      pr_tables_old[k] = nk_tables_old[k]/(n_cost_old+eta_CT);
    pr_tables_old[K_tables_old] = eta_CT/(n_cost_old+eta_CT);
  } 
  lqratio = lpr_q_inv - lpr_q;
  if(local_print){
    cout << "new_tables: ";
    for(int i = 0; i < n_cost_j; i++)
      cout << new_tables[i] << " ";
    cout << endl;
    cout << "new_tables_mainID: ";
    for(int i = 0; i < n_cost_j; i++)
      cout << new_tables_mainID[i] << " ";
    cout << endl;
    cout << "new_dishes: ";
    for(int i = 0; i < n_cost_j; i++)
      cout << new_dishes[i] << " ";
    cout << endl;
  } 

  // let's compute lpratio (part 1)
  for(int it = 0; it < n_cost_j; it++){
    i = cost_locID_j[it];
    new_table_i = new_tables[it];  
    lpratio += log(Pinfo.pr_tables[new_table_i]);

    if(new_table_i < Pinfo.K_tables){
      Pinfo.nk_tables[new_table_i]++;
    } else {
      Pinfo.nk_tables.push_back(1);
      Pinfo.pr_tables.push_back(1);
      Pinfo.K_tables++;
      // Now we choose a new dish (and compute probabilities)
      for(int d = 0; d < Pinfo.K_dishes; d++){
        Pinfo.pr_dishes[d] = Pinfo.nk_dishes[d]/(Pinfo.n_tables + eta_TD);
      }
      Pinfo.pr_dishes[Pinfo.K_dishes] = eta_TD/(Pinfo.n_tables + eta_TD);
      new_dish_i = new_dishes[it];
      lpratio += log(Pinfo.pr_dishes[new_dish_i]);
      
      if(new_dish_i < Pinfo.K_dishes){
        Pinfo.nk_dishes[new_dish_i]++;
      } else {
        Pinfo.nk_dishes.push_back(1);
        Pinfo.pr_dishes.push_back(1);
        Pinfo.K_dishes++;
      }
      Pinfo.n_tables++;
    }
    // if(local_print) cout << new_dishes[it] << endl;
    // Now we compute the probabilities for sampling the next guy.
    Pinfo.n_cost++;
    for(int k = 0; k < Pinfo.K_tables; k++){
      Pinfo.pr_tables[k] = Pinfo.nk_tables[k]/(Pinfo.n_cost+eta_CT);
    }
    Pinfo.pr_tables[Pinfo.K_tables] = eta_CT/(Pinfo.n_cost+eta_CT);
  }
  // let's compute lpratio (part 2)
  for(int it = 0; it < n_cost_j; it++){
    i = cost_locID_j[it];
    if(Pinfo_old.pr_tables[oldtables[it]] > 0){
      lpratio += -log(Pinfo_old.pr_tables[oldtables[it]]);
    } else {
      // if i is the first one to seat at its table, we pretend it's a new one, 
      // and we compute the probability of sampling the right dish
      lpratio += -log(Pinfo_old.pr_tables[Pinfo_old.K_tables]);
      for(int d = 0; d < Pinfo_old.K_dishes; d++){
        Pinfo_old.pr_dishes[d] = Pinfo_old.nk_dishes[d]/(Pinfo_old.n_tables + eta_TD);
      }
      if(Pinfo_old.pr_dishes[olddishes[it]] > 0){
        lpratio += -log(Pinfo_old.pr_dishes[olddishes[it]]);
      } else {
        lpratio += -log( eta_TD/(Pinfo_old.n_tables + eta_TD) ); // == pr_dishes_old[K_dishes_old]
      }
      Pinfo_old.nk_dishes[olddishes[it]]++;
      Pinfo_old.n_tables++;
    }
    Pinfo_old.nk_tables[oldtables[it]]++;
    Pinfo_old.n_cost++;
    for(int k = 0; k < Pinfo_old.K_tables; k++)
      Pinfo_old.pr_tables[k] = Pinfo_old.nk_tables[k]/(Pinfo_old.n_cost+eta_CT);
    Pinfo_old.pr_tables[Pinfo_old.K_tables] = eta_CT/(Pinfo_old.n_cost+eta_CT);
  }   

  // cout << "******" << endl;
  // cout << "lqratio " << lqratio << " lpratio " << lpratio << endl;
  // cout << "******" << endl;

  llikratio = 0.0;
  // we only need the dishes involved in the old or new configuration
  // llikratio for old configuration
  std::vector<int> unique_dishes;
  bool temp_bool;
  for(int it = 0; it < n_cost_j; it++){
    temp_bool = false; // is olddishes[it] contained in unique_dishes?
    for(int d = 0; d < unique_dishes.size(); d++)
      if(unique_dishes[d] == olddishes[it])
        temp_bool = true;
    if(temp_bool == false){ // then it's the first time we see it
      if(local_print) cout << "llikratio for old configuration, olddishes, dish " << olddishes[it] << endl;
      unique_dishes.push_back(olddishes[it]);
      llikratio += - highresPart->dish_llikelihood(olddishes[it], local_print);
    }
  }
  for(int it = 0; it < n_cost_j; it++){
    // if(local_print) cout << "new_dishes " << new_dishes[it] << endl;
    temp_bool = false; // is new_dishes[it] contained in unique_dishes?
    for(int d = 0; d < unique_dishes.size(); d++)
      if(unique_dishes[d] == new_dishes[it])
        temp_bool = true;
    if(temp_bool == false){ // then it's the first time we see it
      if(new_dishes[it] < highresPart->tables_dishes->K){ // only if it existed in the old configuration
        if(local_print) cout << "llikratio for old configuration, new_dishes, dish " << new_dishes[it] << endl;
        unique_dishes.push_back(new_dishes[it]);
        llikratio += - highresPart->dish_llikelihood(new_dishes[it], local_print);
      }
    }
  }
  // if(local_print) cout << "End llikratio for old configuration" << endl;
  // if(local_print) cout << "Size unique_dishes" << unique_dishes.size() << endl;

  // llikratio for new configuration
  int dish_k;
  int table_mainID, rest, table_restID, cost_restID, cost_mainID;
  for(int d = 0; d < unique_dishes.size(); d++){
    // code adapted from Hdp::dish_llikelihood(int dish_k)
    dish_k = unique_dishes[d];
    // if(local_print) cout << endl << endl << "dish_k " << dish_k << endl;
    std::vector<int> costumers_dish_k;
    if(dish_k < highresPart->tables_dishes->K){
      // let's check costumers in the old tables
      for(int t = 0; t < highresPart->tables_dishes->cluster_config[dish_k]; t++){
        table_mainID = highresPart->tables_dishes->clusters[dish_k][t];
        // if(local_print) cout << "table_mainID " << table_mainID << endl;
        Table_MainId_RestId(table_mainID, highresPart->nTableRest, rest, table_restID);

        if((rest != kLR_old) && (rest != kLR_new)){
          for(int c = 0; c < highresPart->costumers_tables[rest]->cluster_config[table_restID]; c++){
            cost_restID = highresPart->costumers_tables[rest]->clusters[table_restID][c];
            cost_mainID = highresPart->costumers_tables[rest]->indexing[cost_restID];
            costumers_dish_k.push_back(cost_mainID);
            // if(local_print)
            //   cout << cost_mainID << " ";
          }
        }
        if(rest == kLR_old){
          if(remove_table[table_restID] == 0){ // if we remove the table it means there are no elements left
            for(int c = 0; c < highresPart->costumers_tables[rest]->cluster_config[table_restID]; c++){
              cost_restID = highresPart->costumers_tables[rest]->clusters[table_restID][c];
              temp_bool = false; // is cost_restID in cost_locID_j?
              for(int it = 0; it < n_cost_j; it++)
                if(cost_locID_j[it] == cost_restID)
                  temp_bool = true;
              if(temp_bool == false){
                cost_mainID = highresPart->costumers_tables[rest]->indexing[cost_restID];
                costumers_dish_k.push_back(cost_mainID);
                // if(local_print)
                //   cout << cost_mainID << " ";
              } // else, we ignore it, since it's in j.
            }
          }
        }
        // here we consider only the costumers that were in kLR_new before the change
        if(rest == kLR_new){
          for(int c = 0; c < highresPart->costumers_tables[rest]->cluster_config[table_restID]; c++){
            cost_restID = highresPart->costumers_tables[rest]->clusters[table_restID][c];
            cost_mainID = highresPart->costumers_tables[rest]->indexing[cost_restID];
            costumers_dish_k.push_back(cost_mainID);
            // if(local_print)
            //   cout << cost_mainID << " ";
          }
          // // are there any elements in j that joined that table?
          // for(int it = 0; it < n_cost_j; it++){
          //   if(new_tables[it] == table_restID){
          //     costumers_dish_k.push_back(cost_indexing_j[it]);
          //   }
          // }
        }
        // if(local_print)
        //   cout << endl << endl;
      }
      // let's check the costumers in j if they are now assigned to dish_k
    } // if the dish was not in the old configuration, then it's only in new_dishes
    // if(local_print) cout << "among new_dishes: ";
    for(int it = 0; it < n_cost_j; it++){
      if(new_dishes[it] == dish_k){
        costumers_dish_k.push_back(cost_indexing_j[it]);
        // if(local_print)
        //   cout << cost_indexing_j[it] << " ";
      }
    }
    // if(local_print) cout << endl;
    
    int n = costumers_dish_k.size();
    if(n > 0){
      arma::vec y_cl = zeros<vec>( costumers_dish_k.size() );
      for(int c = 0; c < costumers_dish_k.size(); c++){
        y_cl(c) = Y( costumers_dish_k[c] );
      }

      double ybar = arma::mean(y_cl);
      double k_n = k_H + n;
      double alpha_n = alpha_H + n/2;
      double beta_n = beta_H + sum(pow(y_cl-ybar,2)) + k_H * n * pow(ybar-mu_H,2)/(2*k_n);
      double loglike = lgamma(alpha_n) - lgamma(alpha_H) + alpha_H * log(beta_H) - alpha_n * log(beta_n);
      loglike += 0.5*(log(k_H)-log(k_n)) -n/2*log(2*M_PI);
      
      llikratio += loglike;
    }
  }
  
  info->cost_locID_j = cost_locID_j;
  info->cost_locID_NOTj = cost_locID_NOTj;
  info->cost_indexing_j = cost_indexing_j;
  info->oldtables = oldtables;
  info->oldtables_mainID = oldtables_mainID;
  info->olddishes = olddishes;
  info->n_cost_j = n_cost_j;
  info->remove_table = remove_table;
  info->remove_dish = remove_dish;
  info->new_tables = new_tables; // *****
  // info->new_tables_mainID = new_tables_mainID;
  info->new_dishes = new_dishes;
  return;
}

void Multires::save_proposal(proposal_info* info, double eta_LR, double eta_CT, double eta_TD, bool local_print){
  // we first change LR
  int nc1 = lowresPart->cluster_config[ info->kLR_old ] - 1;
  int nc2 = 1;
  int *index1, *index2;
  index1 = new int[nc1];
  index2 = new int[nc2];
  int temp, count = 0;
  for(int i = 0; i < (nc1+1); i++){
    temp = lowresPart->clusters[info->kLR_old][i];
    if(temp != info->jLR){
      index1[count] = temp;
      count++;
    } else {
      index2[0] = info->jLR;
    }
  }
  if(nc1 == 0){
    if(info->kLR_new < lowresPart->K){
      // what happens here? that the largest is merged into the smallest
      lowresPart->Merge(info->kLR_old, info->kLR_new, eta_LR);
    } else {
      // do nothing!!
    }
  } else {
    if(info->kLR_new < lowresPart->K){
      lowresPart->Split_and_Merge(info->kLR_old, index1, index2, nc1, nc2, info->kLR_old, info->kLR_new, eta_LR);
    } else {
      lowresPart->Split(info->kLR_old, index1, index2, nc1, nc2, eta_LR);
    }
  }
  delete[] index1;
  delete[] index2;
  // we now change HR
  highresPart->Change_Rest(info, eta_CT, eta_TD, local_print);
}

double Multires::get_tables_merge(int rest, int threshold, std::default_random_engine& myRandomEngine, double eta_CT, double eta_TD){
  /*
  This function makes a merge of the tables (before they had just been collected together in the merged restaurant)
  but only if they are assigned to the same dish. 
  Thus, ntables_per_dish1/2 count how many tables are assigned to each dish
  and which_tables_per_dish is a sort of partTD->clusters but containing the tables localID
  Then, for each dish, we compute the possible number of pairs between the two sets of tables 
  (not choose(k,N) because it does not account for the order), where k<N
    this happens with probability 1/(n!/(n-k)!) and given the order with merge each pair with probability 0.5
  tobemerged1/2 contains the localID of the tables which will get merged 
    (they have the same length and contain each the index of one element of the pair)
    note: tobemerged1 contains someone from list1.
  For each pair, the merge is actually made by Merge_CostTable (using tobemerged1[i], tobemerged2[i])
  lqratio is computed and returned. 
  */
  // the tables (in the merged restaurant) that have localID < threshold, they were in list1
  double lqratio = 0.0;
  LPPartition part = highresPart->costumers_tables[rest];
  LPPartition partTD = highresPart->tables_dishes;
  // we need to figure out the dishes
  int tMainID, dish;
  std::vector<int> ntables_per_dish1(partTD->K, 0);
  std::vector<int> ntables_per_dish2(partTD->K, 0);
  std::vector<std::vector<int> > which_tables_per_dish(partTD->K); // should you use std::vector<int>() ?
  int shift = Table_RestId_MainId(rest, 0, highresPart->nTableRest);
  for(int t = 0; t < part->K; t++){
    tMainID = shift + t;
    // if(t == 0){
    //   tMainID = Table_RestId_MainId(rest, 0, highresPart->nTableRest);
    // } else {
    //   tMainID++;
    // }
    dish = partTD->cluster_assignment[tMainID]; 
    which_tables_per_dish[dish].push_back(t);
    if(t < threshold){
      ntables_per_dish1[dish]++;
    } else {
      ntables_per_dish2[dish]++;
    }
  }
  int t1, t2;
  arma::vec P(2);
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  std::vector<int> tobemerged1,tobemerged2;
  int N, k;
  for(int dish = 0; dish < partTD->K; dish++){
    if(ntables_per_dish1[dish]>0 && ntables_per_dish2[dish]>0){
      if(ntables_per_dish1[dish] <= ntables_per_dish2[dish]){
        k = ntables_per_dish1[dish]; 
        N = ntables_per_dish2[dish];
        // I changed it so I could remove the shuffling at the end (default: do shuffle. ie. here we shuffle so the order matters)
        std::vector<int> t2s = sample_without_replacement(k, N, myRandomEngine);
        // 1/Binomialcoeff(N,k): NO it should be n!/(n-k)! because I care about the order of the k tables
        // lqratio += -(lgamma(N+1) - lgamma(k+1) - lgamma(N-k+1));
        lqratio += -(lgamma(N+1) - lgamma(N-k+1));
        for(int i = 0; i < t2s.size(); i++){
          // since which_tables_per_dish[dish] is ordered, I can find the elements of the second list as ntables_per_dish1[dish]+t2
          t1 = which_tables_per_dish[dish][i];
          t2 = which_tables_per_dish[dish][ntables_per_dish1[dish]+t2s[i]];
          // given t1,t2, we compute the probability of merging vs the probability of keeping them separate.
          
          // P[0] = lgamma(part->cluster_config[t1]+part->cluster_config[t2]); // merged
          // P[1] = log(eta_CT) + lgamma(part->cluster_config[t1]) + lgamma(part->cluster_config[t2]); // separate
          // P = P - log(sum(exp(P)));
          P[0] = log(0.5); // merged
          P[1] = log(0.5); // separate
          if(distribution(myRandomEngine) < exp(P[0])){
            lqratio += P[0];
            // we merge
            tobemerged1.push_back(t1);
            tobemerged2.push_back(t2);
          } else {
            lqratio += P[1];
          }
        }
      } else {
        k = ntables_per_dish2[dish]; 
        N = ntables_per_dish1[dish];
        std::vector<int> t1s = sample_without_replacement(k, N, myRandomEngine);
        // lqratio += -(lgamma(N+1) - lgamma(k+1) - lgamma(N-k+1));
        lqratio += -(lgamma(N+1) - lgamma(N-k+1));
        for(int i = 0; i < t1s.size(); i++){
          t2 = which_tables_per_dish[dish][ntables_per_dish1[dish]+i];
          t1 = which_tables_per_dish[dish][t1s[i]];
          // given t1,t2, we compute the probability of merging vs the probability of keeping them separate.
          
          // P[0] = lgamma(part->cluster_config[t1]+part->cluster_config[t2]); // merged
          // P[1] = log(eta_CT) + lgamma(part->cluster_config[t1]) + lgamma(part->cluster_config[t2]); // separate
          // P = P - log(sum(exp(P)));
          P[0] = log(0.5); // merged
          P[1] = log(0.5); // separate
          if(distribution(myRandomEngine) < exp(P[0])){
            lqratio += P[0];
            // we merge
            tobemerged1.push_back(t1);
            tobemerged2.push_back(t2);
          } else {
            lqratio += P[1];
          }
        }
      }
    }
  } // end for
  // now we actually merge
  for(int i = 0; i < tobemerged1.size(); i++){
    highresPart->Merge_CostTable(rest, tobemerged1[i], tobemerged2[i], eta_CT, eta_TD);
    // after we merge one pair, one table from the ones which were in rmax "disappears" and 
    // we need to shift the names of the ones that came after. This is done for tobemerged2 only,
    // because the table from rmax gets merged into the one in rmin, and only for the ones that 
    // had names larger than the one which just got merged. Also, only for j > i because the 
    // ones before have already been merged so we don't care.
    for(int j = i+1; j < tobemerged1.size(); j++){
      if(tobemerged2[j] > tobemerged2[i]){
        tobemerged2[j] += -1;
      }
    }
  }
  return lqratio;
}

double Multires::get_tables_merge_reverse(int rest1, int rest2, LPPartition oldpart, std::default_random_engine& myRandomEngine){
  // used after the restaurant has been split
  /*
  This function computes a potential merge that would have happened to go back from rest1, rest2 to the original oldpart
  It does not actually make merges, but computes lqratio2.
  Thus, ntables_per_dish1/2 count how many tables in rest1/2 are assigned to each dish
  and which_tables_per_dish is a sort of partTD->clusters but containing the tables localID
  (since we use localID for both rest1/2, how do we find if t refers to rest1 or 2? Well, note that
  all the rest1 come before the rest2, so the first ntables_per_dish1[dish] correspond to rest1 and the ones after to rest2)

  Then, for each dish, we compute the possible number of pairs between the two sets of tables 
  (not choose(k,N) because it does not account for the order), where k<N
    this happens with probability 1/(n!/(n-k)!) and given the order with merge each pair with probability 0.5
    (we check if in the old partitions the two tables were in fact merged or not)
  lqratio is computed and returned. 
  */
  double lqratio = 0.0;
  LPPartition part1 = highresPart->costumers_tables[rest1];
  LPPartition part2 = highresPart->costumers_tables[rest2];
  LPPartition partTD = highresPart->tables_dishes;
  // we need to figure out the dishes
  int tMainID, dish;
  // how many tables that belong to list1/list2 are in each dish?
  std::vector<int> ntables_per_dish1(partTD->K, 0);
  std::vector<int> ntables_per_dish2(partTD->K, 0);
  std::vector<std::vector<int> > which_tables_per_dish(partTD->K); // should you use std::vector<int>() ?
  int shift = Table_RestId_MainId(rest1, 0, highresPart->nTableRest);
  for(int t = 0; t < part1->K; t++){
    tMainID = shift + t;
    // if(t == 0){
    //   tMainID = Table_RestId_MainId(rest1, 0, highresPart->nTableRest);
    // } else {
    //   tMainID++;
    // }
    dish = partTD->cluster_assignment[tMainID]; 
    which_tables_per_dish[dish].push_back(t);
    ntables_per_dish1[dish]++;
  }
  // tables in part2 are recognizable because their position in which_tables_per_dish[dish] is greater than ntables_per_dish1[dish]
  shift = Table_RestId_MainId(rest2, 0, highresPart->nTableRest);
  for(int t = 0; t < part2->K; t++){
    tMainID = shift + t;
    // if(t == 0){
    //   tMainID = Table_RestId_MainId(rest2, 0, highresPart->nTableRest);
    // } else {
    //   tMainID++;
    // }
    dish = partTD->cluster_assignment[tMainID]; 
    which_tables_per_dish[dish].push_back(t);
    ntables_per_dish2[dish]++;
  }
  int t1, t2;
  arma::vec P(2);
  int ind1, ind2, i1, i2;
  int k, N;
  for(int dish = 0; dish < partTD->K; dish++){
    if(ntables_per_dish1[dish]>0 && ntables_per_dish2[dish]>0){
      if(ntables_per_dish1[dish] <= ntables_per_dish2[dish]){
        k = ntables_per_dish1[dish]; 
        N = ntables_per_dish2[dish];
        std::vector<int> t2s = sample_without_replacement(k, N, myRandomEngine);
        // 1/(n!/(n-k)!) instead of 1/bin-coeff
        // lqratio += - ( lgamma(N+1) - lgamma(k+1) - lgamma(N-k+1) );
        lqratio += - ( lgamma(N+1) - lgamma(N-k+1) );
        for(int i = 0; i < t2s.size(); i++){
          // because which_tables_per_dish[dish] is ordered, I can find the elements of the second list as ntables_per_dish1[dish]+t2
          t1 = which_tables_per_dish[dish][i];
          t2 = which_tables_per_dish[dish][ntables_per_dish1[dish]+t2s[i]];
          // t1/2 is localID for a table in part1/2 
          // given t1,t2, we compute the probability of merging vs the probability of keeping them separate.
          
          // P[0] = lgamma(part1->cluster_config[t1]+part2->cluster_config[t2]); // merged
          // P[1] = log(eta_CT) + lgamma(part1->cluster_config[t1]) + lgamma(part2->cluster_config[t2]); // separate
          // P = P - log(sum(exp(P)));
          P[0] = log(0.5); // merged
          P[1] = log(0.5); // separate
          // checking if the two tables were together in the old partition
          ind1 = part1->indexing[ part1->clusters[t1][0] ];
          ind2 = part2->indexing[ part2->clusters[t2][0] ];
          for(int i = 0; i < oldpart->nObs; i++){
            if(oldpart->indexing[i] == ind1){
              i1 = i; // local index in oldpar that corresponds to a costumer seated at t1
            } else {
              if(oldpart->indexing[i] == ind2){
                i2 = i;
              }
            }
          }
          if(oldpart->cluster_assignment[i1] == oldpart->cluster_assignment[i2]){
            // if they were together, consider the probability of merging.
            lqratio += P[0];
          } else {
            lqratio += P[1];
          }
        }
      } else {
        k = ntables_per_dish2[dish]; 
        N = ntables_per_dish1[dish];
        std::vector<int> t1s = sample_without_replacement(k, N, myRandomEngine);
        // lqratio += -(lgamma(N+1) - lgamma(k+1) - lgamma(N-k+1));
        lqratio += - ( lgamma(N+1) - lgamma(N-k+1) );
        for(int i = 0; i < t1s.size(); i++){
          // because which_tables_per_dish[dish] is ordered, I can find the elements of the second list as ntables_per_dish1[dish]+t2
          t2 = which_tables_per_dish[dish][ntables_per_dish1[dish]+i];
          t1 = which_tables_per_dish[dish][t1s[i]];
          // given t1,t2, we compute the probability of merging vs the probability of keeping them separate.
          
          // P[0] = lgamma(part1->cluster_config[t1]+part2->cluster_config[t2]); // merged
          // P[1] = log(eta_CT) + lgamma(part1->cluster_config[t1]) + lgamma(part2->cluster_config[t2]); // separate
          // P = P - log(sum(exp(P)));
          P[0] = log(0.5); // merged
          P[1] = log(0.5); // separate
          // checking if the two tables were together in the old partition
          ind1 = part1->indexing[ part1->clusters[t1][0] ];
          ind2 = part2->indexing[ part2->clusters[t2][0] ];
          for(int i = 0; i < oldpart->nObs; i++){
            if(oldpart->indexing[i] == ind1){
              i1 = i;
            } else {
              if(oldpart->indexing[i] == ind2){
                i2 = i;
              }
            }
          }
          if(oldpart->cluster_assignment[i1] == oldpart->cluster_assignment[i2]){
            lqratio += P[0];
          } else {
            lqratio += P[1];
          }
        }
      }
    }
  } // end for
  // cout << "** " << endl;
  return lqratio;
}

double Multires::splitted_tables(int rest, std::list<int> new_c1_ind, std::list<int> new_c2_ind, double eta_CT){
  // This is used for split, when everyone in list1 and list2 belongs to rest
  // Assuming the cluster is covered by list1 and list2 
  LPPartition part = highresPart->costumers_tables[rest];
  LPPartition part_dish = highresPart->tables_dishes;
  // mainID correspond to one blockgroup (or costumer) and ctMainID is its corresponding census tract
  int mainID, ctMainID;
  std::list<int>::iterator it;
  double ll = 0;
  // each element contains how many costumers (of list1/2) are in each table
  // or: how many blockgroups in the census tracts of list (risp. list2) are in each table
  std::vector<int> subtables1(part->K, 0); 
  std::vector<int> subtables2(part->K, 0);
  std::vector<int> moretables(part_dish->K, 0); 
  bool nonzero1, nonzero2;
  int n1 = 0, n2 = 0;
  int tMainID, dish;
  for(int k = 0; k < part->K; k++){
    // consider table k in restaurant corresponding to the cluster of census tracts
    nonzero1 = false; // tells me if in table k there is at least one blockgroup whose census tract is in list1
    nonzero2 = false;
    for(int i = 0; i < part->cluster_config[k]; i++){
      // consider ith person in table k;
      mainID = part->indexing[ part->clusters[k][i] ];
      // ctMainID tells me the census tract in which this person (= block group) belongs
      ctMainID = InvMapping(mainID, 1);
      // Is the census tract in list1 or list2? (the lists contain the ctmainID)
      it = std::find(new_c1_ind.begin(), new_c1_ind.end(), ctMainID);
      if(it != new_c1_ind.end()){ // element exists in list1
        subtables1[k]++;
        nonzero1 = true;
      } else {
        subtables2[k]++;
        nonzero2 = true;
      }
    }
    if(nonzero1 & nonzero2){
      tMainID = Table_RestId_MainId(rest, k, highresPart->nTableRest);
      dish = part_dish->cluster_assignment[tMainID];
      moretables[dish]++;
    }
    if(nonzero1){
      n1 += subtables1[k];
      ll += lgamma(subtables1[k]) + log(eta_CT);
    }
    if(nonzero2){
      n2 += subtables2[k];
      ll += lgamma(subtables2[k]) + log(eta_CT);
    }
  }
  // this is the normalizing constant
  ll += -lgamma(eta_CT + n1) + lgamma(eta_CT); 
  ll += -lgamma(eta_CT + n2) + lgamma(eta_CT);
  // NOW we'll consider the prior of tables into dishes
  for(int d = 0; d < part_dish->K; d++){
    ll += lgamma(part_dish->cluster_config[d] + moretables[d]) - lgamma(part_dish->cluster_config[d]);
  }
  return ll;
}

double Multires::splitted_tables(int rest1, int rest2, std::list<int> new_c1_ind, std::list<int> new_c2_ind, double eta_CT){
  // This is used for merge, when everyone in list1 and list2 belongs to either rest1 or rest2
  LPPartition part1 = highresPart->costumers_tables[rest1];
  LPPartition part2 = highresPart->costumers_tables[rest2];
  
  // part1->Print_Partition_Short();
  // part2->Print_Partition_Short();
  // for(std::list<int>::iterator it=new_c1_ind.begin(); it != new_c1_ind.end(); ++it){
  //   cout << *it << " ";
  // }
  // cout << endl;
  // for(std::list<int>::iterator it=new_c2_ind.begin(); it != new_c2_ind.end(); ++it){
  //   cout << *it << " ";
  // }
  // cout << endl;
  // cout << endl;

  // mainID correspond to one blockgroup (or costumer) and ctMainID is its corresponding census tract
  int mainID, ctMainID;
  std::list<int>::iterator it;
  double ll = 0;
  // each element contains how many costumers (of list1/2) are in each table
  std::vector<int> subtables1(part1->K + part2->K, 0); 
  std::vector<int> subtables2(part1->K + part2->K, 0);
  bool nonzero1, nonzero2;
  int n1 = 0, n2 = 0;
  for(int k = 0; k < part1->K; k++){
    nonzero1 = false;
    nonzero2 = false;
    for(int i = 0; i < part1->cluster_config[k]; i++){
      mainID = part1->indexing[ part1->clusters[k][i] ];
      ctMainID = InvMapping(mainID, 1);
      // the lists contain the mainID
      it = std::find(new_c1_ind.begin(), new_c1_ind.end(), ctMainID);
      if(it != new_c1_ind.end()){ // element exists in list
        subtables1[k]++;
        nonzero1 = true;
      } else {
        subtables2[k]++;
        nonzero2 = true;
      }
    }
    if(nonzero1){
      n1 += subtables1[k];
      ll += lgamma(subtables1[k]) + log(eta_CT);
    }
    if(nonzero2){
      n2 += subtables2[k];
      ll += lgamma(subtables2[k]) + log(eta_CT);
    }
  }
  for(int k = 0; k < part2->K; k++){
    nonzero1 = false;
    nonzero2 = false;
    for(int i = 0; i < part2->cluster_config[k]; i++){
      mainID = part2->indexing[ part2->clusters[k][i] ];
      ctMainID = InvMapping(mainID, 1);
      // the lists contain the mainID
      it = std::find(new_c1_ind.begin(), new_c1_ind.end(), ctMainID);
      if(it != new_c1_ind.end()){ // element exists in list
        subtables1[part1->K + k]++;
        nonzero1 = true;
      } else {
        subtables2[part1->K + k]++;
        nonzero2 = true;
      }
    }
    if(nonzero1){
      n1 += subtables1[part1->K + k];
      ll += lgamma(subtables1[part1->K + k]) + log(eta_CT);
    }
    if(nonzero2){
      n2 += subtables2[part1->K + k];
      ll += lgamma(subtables2[part1->K + k]) + log(eta_CT);
    }
  }
  ll += -lgamma(eta_CT + n1) + lgamma(eta_CT);
  ll += -lgamma(eta_CT + n2) + lgamma(eta_CT);
  return ll;
}

double Multires::get_logposterior(double eta_LR, double eta_CT, double eta_TD){
  double res = 0.0;
  res += lowresPart->get_logprior(eta_LR);
  res += highresPart->get_logprior(eta_CT, eta_TD);
  res += alpha1 * highresPart->get_loglike();
  if(yLR_bool) res += alpha2 * lowresPart->get_loglike(true);
  return res;
}

// double Multires::tables_llikelihoood(int k, LPPartition i_alone){
//   double ll = 0;
//   double ll1 = 0;
//   double ll2 = 0;
//   int ncl;
//   int n = i_alone->nObs + highresPart->costumers_tables[k]->nObs;
//   for(int h = 0; h < highresPart->costumers_tables[k]->K){
//     ncl = highresPart->costumers_tables[k]->cluster_config[h];
//     ll += log(eta) + lgamma(ncl);
//   }
//   for(int h = 0; h < i_alone->K){
//     ncl = i_alone->cluster_config[h];
//     ll += log(eta) + lgamma(ncl);
//   }
//   ll += -lgamma(eta + nObs) + lgamma(eta);
//   return ll;
// }

// void Multires::ConditionalGibbs_lowres(int i, std::default_random_engine& eng){
//   int current_cl = lowresPart->cluster_assignment[i];
//   int ncl = lowresPart->cluster_config[current_cl];
//   int K = lowresPart->K;
//   arma::vec P, l, logpr;
//   if(ncl == 1){
//     P = zeros<vec>(K);
//     l = zeros<vec>(K);
//     for(int k = 0; k < K; k++){
//       if(k < current_cl){
//         double pp,ll;
//         pp = lowresPart->cluster_config[k]; 
//         ll = tables_llikelihoood(k, i_alone); //dish_llikelihood_plusCost(dish,iMainID)-dish_llikelihood(dish);
//         P[k] = pp;
//         l[k] = ll;
//       } else if(k > current_cl){
//         double pp,ll;
//         pp = cluster_config[k]; 
//         ll = something; //dish_llikelihood_plusCost(dish,iMainID)-dish_llikelihood(dish);
//         P[k-1] = pp;
//         l[k-1] = ll;
//       }
//     }
//     P[restaurant->K -1] = eta;
//     l[restaurant->K -1] = something; // newdish_llikelihood(iMainID, idish, newdish_pr);
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
//       Merge(itable, sampled); // merge will fix which one is larger
//       // now we really need to CHANGE the tables
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
