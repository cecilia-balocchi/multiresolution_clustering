#ifndef PARTITION_H_
#define PARTITION_H_
#include <armadillo>
#include <string.h>
#include <list>
#include <set>
#include <random>
using namespace std;

typedef class Partition* LPPartition;
class Partition
{
  //data members
public:
  int nObs; // number of indices
  int K; // number of clusters
  int* cluster_config; // sizes of each cluster
  int** clusters; // the actual clusters  // might be easier to use vectors.
  int* cluster_assignment; // size nObs, tells us which cluster each index belongs to
  int** pairwise_assignment; // array of pairwise assignments
  double* log_prior; // size K. holds the log-prior evaluated on each cluster
  // double* beta_hat; // will have size nObs. holds the point estimates of beta in each blockgroup

  std::vector<int> indexing; // from local ID to main ID (the index corresponds to local, the value to main)
  // double* log_det_Omegay;
  // double* quad_forms;
  // methods
public:
  Partition(); // constructor
  Partition(LPPartition initial_partition); // constructor. create a new partition from an old one
  ~Partition(); // destructor
public:
  void Initialize_Partition(int n, double eta);
  // void Initialize_Partition(int n, Rcpp::List gamma_init);
  void Initialize_Partition(int n, std::vector<int> ids, double eta);
  void Initialize_Partition(int n, std::vector<int> ids, double eta, int gamma_ptr[], int* K_init_ptr);
  std::vector<int> Initialize_Partition(int n, std::vector<int> ids, double eta, int gamma_ptr[]);
  void Initialize_Partition(int n, double eta, std::vector<int> gamma_ptr, std::vector<int> gamma_unique);
  void Initialize_Partition_nclusters(int n, double eta);
  void Initialize_Partition_nclusters(int n, std::vector<int> ids, double eta);
  void Initialize_CostDish(std::vector<int> nObsRest, std::vector<int> nTableRest, std::vector<LPPartition> costumers_tables, LPPartition tables_dishes);
  void Copy_Partition(LPPartition initial_partition);
  void get_pairwise();
  void log_pi_ep(int cluster_id, double eta);
  double get_logprior(double eta);
  // double get_logprior2(double eta);
  // void get_likelihood(int cluster_id);
  // void beta_postmean(int cluster_id);
  void Print_Partition();
  void Print_Partition_Short();
  void Print_Partition_ToFile(string file);
  void Read_Partition_FromFile(string file, int n);
  void Print_Y(bool LR = false);
  void AddTable_TableDish(int r, int newtable_id, int cluster_newtable, std::vector<int> nTableRest, double eta_TD);
  void MoveTables_TableDish(std::vector<int> tables, int starting_index);
  void SplitTables_TableDish(int r1, int r2, std::vector<int> tables1, std::vector<int> tables2, std::vector<int> nTableRest, double eta_TD, bool local_print = false);
  void AddTables_TableDish(int r, std::vector<int> newtables_id, std::vector<int> cluster_newtables, std::vector<int> nTableRest);
  void RemoveTable_TableDish(int r, int oldtable_mainid, double eta_TD);
  void Split(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, double eta); // split cluster split_k into two parts: new_cluster1, new_cluster2
  void KSplit(int split_k, int num_splits, std::vector<std::vector<int> > indices, std::vector<int> ns, double eta); // split cluster split_k into num_split parts
  void Merge(int k_1, int k_2, double eta); // merge cluster max(k_1, k_2) into min(k_1, k_2)
  void Split_and_Merge(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int k_star_1, int k_star_2, double eta);
  void Modify(int cl_ind);
  // void Find_Splits(int cluster_id, int **index1, int **index2, int &n1, int &n2);
  
  std::vector<int> SampleSM(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR);
  void ProposalSM(LPPartition startingP, int ngibbs, double& lqratio, double& lpratio, double& llikratio, std::default_random_engine& myRandomEngine, std::vector<int>& pars, double eta_LR);
  void LaunchSM(LPPartition startingP, int i1, int i2, std::vector<int> samecl_index, int c1, int c2, std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, 
    std::vector<int>& samecl_newcl, int ngibbs, std::default_random_engine& myRandomEngine);
  void ProposalSM_iter(bool SplitMerge, LPPartition startingP, std::vector<int>& samecl_index, int new_c1, int new_c2, std::list<int>& new_c1_ind, std::list<int>& new_c2_ind,
    std::vector<int>& samecl_newcl, double& lqratio, std::default_random_engine& myRandomEngine);
  
  double get_logposterior(double eta_LR);
  double cluster_llikelihood(std::list<int> index, bool LR = true);
  double cluster_llikelihood(int cl, bool LR = true);
  double get_loglike(bool LR = true);
  bool FindKMSplit(int split_k, arma::vec y, double eta);
  std::vector<int> get_samecl_index(int i1, int i2, int c1, int c2);
  void Combine_Restaurants(LPPartition part1, LPPartition part2);
  void Divide_Restaurants(LPPartition newpart, std::vector<int> index1, std::vector<int> index2, double eta_CT, bool local_print = false);
  void format_partition(double* vector, int starting_index = 0);
  void get_postmean(double* vector, int starting_index = 0, bool LR = true);
  // Rcpp::List format_partition();
};



#endif /* PARTITION_H_ */
