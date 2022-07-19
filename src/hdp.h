#ifndef HDP_H_
#define HDP_H_

#include <vector>
#include <random>
#include "various_functions.h"
using namespace std;

typedef class Hdp* LPHdp;
class Hdp
{
  //data members
public:
  int nObs; // number of costumer over all the restaurantes
  int nRest; // number of restaurants
  std::vector<int> nObsRest; // vector containing the number of people in each restaurant
  std::vector<int> nTableRest; // vector containing the number of tables in each restaurant
  std::vector<LPPartition> costumers_tables; // how the costumers are divided into tables, one for each restaurant
  LPPartition tables_dishes; // how the tables are divided into dishes
  // int K; // this could be the size of the overall partition costumers_dishes
  std::vector<std::vector<int> > CTinRest;
  std::vector<double> dishes_thetas;
  // methods
public:
  Hdp(); // constructor
  Hdp(LPHdp initial_hdp); // constructor. create a new hdp from an old one
  ~Hdp(); // destructor
public:
  void Copy_Hdp(LPHdp initial_hdp);
  void Initialize_Hdp(int n, vector<vector<int> > ctinrest, vector<vector<int> > rests, double eta_CT, double eta_TD);
  void Initialize_Hdp_nRestDish(int n, vector<vector<int> > ctinrest, vector<vector<int> > rests, double eta_CT, double eta_TD);
  void Initialize_Hdp_nclusters(int n, vector<vector<int> > ctinrest, vector<vector<int> > rests, double eta_CT, double eta_TD);
  void Initialize_Hdp_part(int n, vector<vector<int> > ctinrest, vector<vector<int> >  rests, int part_init_ptr[], int* K_init_ptr, double eta_CT, double eta_TD);
  void Initialize_Hdp_part_nrest(int n, vector<vector<int> > ctinrest, vector<vector<int> >  rests, int part_init_ptr[], double eta_CT, double eta_TD);
  
  void Print_Hdp_Short();
  void Print_Hdp();
  void Print_Y();

  void Split_CostTable(int r, int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int cluster_newtable);
  void Merge_CostTable(int r, int merge_k1, int merge_k2, double eta_CT, double eta_TD);
  void SplitMerge_CostTable(int r, int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int merge_k1, int merge_k2);
  void KSplit_CostTable(int r, int split_k, int num_splits, std::vector<std::vector<int> > indices, std::vector<int> ns, std::vector<int> cluster_newtables);
  
  void Split_TableDish(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2);
  void Merge_TableDish(int merge_k1, int merge_k2);
  void SplitMerge_TableDish(int split_k, int* new_cluster1, int* new_cluster2, int size1, int size2, int k_star_1, int k_star_2);
  void KSplit_TableDish(int split_k, int num_splits, std::vector<std::vector<int> > indices, std::vector<int> ns);
  
  std::vector<bool> SampleSM(int ngibbs, double temp, std::default_random_engine& eng, double eta_CT, double eta_TD, bool local_print=false);
  bool SampleSM_TableDishes(int gibbs, double temp, std::default_random_engine& eng, double eta_TD);
  bool SampleSM_CostTable(int r, int ngibbs, double temp, std::default_random_engine& eng, double eta_CT, double eta_TD, bool local_print=false);
  void ProposalSM_TableDishes(LPHdp startingH, int ngibbs, double& lqratio, double& lpratio, double& llikratio, std::default_random_engine& myRandomEngine, double eta_TD);
  void ProposalSM_CostTable(LPHdp startingH, int rest, int ngibbs, double& lqratio, double& lpratio, double& llikratio, std::default_random_engine& myRandomEngine, double eta_CT, double eta_TD, bool local_print=false);
  void LaunchSM_TableDishes(LPHdp startingH, int i1, int i2, std::vector<int> samecl_index, int c1, int c2, std::list<int>& new_c1_ind, 
    std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, int ngibbs, std::default_random_engine& gen);
  void LaunchSM_CostTable(LPHdp startingH, int rest, int params[], std::vector<int> samecl_index, 
    std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, int ngibbs, std::default_random_engine& gen);
  void ProposalSM_TableDishes_iter(bool SplitMerge, LPHdp startingH, std::vector<int>& samecl_index, int new_c1, int new_c2,
    std::list<int>& new_c1_ind, std::list<int>& new_c2_ind,std::vector<int>& samecl_newcl, double& lqratio, std::default_random_engine& myRandomEngine);
  void ProposalSM_CostTable_iter(bool SplitMerge, LPHdp startingH, int rest, int params[], std::vector<int> samecl_index,
    std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, double& lqratio, std::default_random_engine& myRandomEngine);
  
  void SampleGibbs_CostTable(int r, std::default_random_engine& eng);
  void ConditionalGibbs_CostTable(int r, int i, std::default_random_engine& eng);
  
  std::vector<double> get_pr(int i, int rest);
  double dish_avg(int dish_k);
  double dish_llikelihood(int dish_k, bool local_print = false);
  double dish_llikelihood(std::list<int> index);
  double table_llikelihood(std::list<int> index, int rest);
  double get_logprior(double eta_CT, double eta_TD);
  // double get_logprior2();
  double get_loglike();
  // double dish_llikelihood_plusCost(int dish_k, int i);
  // double dish_llikelihood_minusCost(int dish_k, int i);
  // double dish_llikelihood_plusminusCosts(int dish_k, std::list<int> index1, std::list<int> index2, int rest);
  // double newdish_llikelihood(int i, int current_dish, std::vector<double>& pr);
  // double cost_llikelihood(int i);
  
  void get_CostDish(LPPartition costdish);
  std::vector<double> get_dishes_thetas();
  void Merge_Rest(int r1, int r2);
  void Split_Rest(int r1, int r2, std::vector<int> index1, std::vector<int> index2, double eta_CT, double eta_TD, double local_print = false);
  void Change_Rest(proposal_info* info, double eta_CT, double eta_TD, double local_print = false);
  void format_rest_table(double* vector_rest, double* vector_table, int starting_index = 0);
};



#endif /* HPD_H_ */