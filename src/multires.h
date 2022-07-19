#ifndef MULTIRES_H_
#define MULTIRES_H_

#include <vector>
#include "partition.h"
#include "hdp.h"
#include "various_functions.h"
using namespace std;

typedef class Multires* LPMultires;
class Multires
{
public:
  // Problem: restaurant is used for both stable things (like CT) and unstable (like clusters)
  // Do I need a copy of Mapping for each object? no. 
  int nObs_lowres;
  int nObs_highres;
  std::vector<int> rest_config; // number of costumers in each restaurant CAREFUL: NEVER UPDATED OR USED
  std::vector< std::vector<int> > rests; // tells me which bg is in each cs (or which costumer is in which census tract)
  LPPartition lowresPart;
  LPHdp highresPart;
public:
  Multires(); // constructor
  Multires(LPMultires initial_multires); // constructor. create a new hdp from an old one
  ~Multires(); // destructor
public:
  void Copy_Multires(LPMultires initial_multires);
  void Initialize_Multires(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, double eta_LR, double eta_CT, double eta_TD);
  void Initialize_Multires(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, LPPartition lowres, LPHdp highres);
  void Initialize_Multires_nclLR_oneHR(int nObsL, int nObsH, std::vector<std::vector<int> > nj, double eta_LR, double eta_CT, double eta_TD);
  void Initialize_Multires_oneLR_nclHR(int nObsL, int nObsH, std::vector<std::vector<int> > nj, double eta_LR, double eta_CT, double eta_TD);
  void Initialize_Multires_oneLR_partHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, int part_init_ptr[], int* K_init_ptr, double eta_LR, double eta_CT, double eta_TD);
  void Initialize_Multires_nclLR_partHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, int part_init_ptr[], double eta_LR, double eta_CT, double eta_TD);
  void Initialize_Multires_partLR_oneHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, int part_init_ptr[], int* K_init_ptr, double eta_LR, double eta_CT, double eta_TD);
  void Initialize_Multires_partLR_nclHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, int part_init_ptr[], int* K_init_ptr, double eta_LR, double eta_CT, double eta_TD);
  void Initialize_Multires_partLR_partHR(int nObsL, int nObsH, std::vector<std::vector<int> > Rests, int partLR_init_ptr[], int* KLR_init_ptr, int partHR_init_ptr[], double eta_LR, double eta_CT, double eta_TD);
  void Initialize_Multires_nclLR_nclHR(int nObsL, int nObsH, std::vector<std::vector<int> > nj, double eta_LR, double eta_CT, double eta_TD);
  void Print_Multires();

  void Get_Highres_CostDish(LPPartition highres_costdish);
  // std::vector<bool> SampleSM(int ngibbs, std::default_random_engine& eng, bool sampleLR = true);
  std::vector<bool> SampleSM(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD, bool sampleLR = true, bool sampleHR = true, bool local_print = false);
  bool SampleSM_LR(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD, bool local_print = false);
  bool SampleSM_LR_noHR(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD);
  // bool SampleSM_LR(int ngibbs, std::default_random_engine& eng);
  void ProposalSM_LR(int ngibbs, double& lqratio, double& lqratio2, double& lpratio, double& llikratio, std::default_random_engine& myRandomEngine, double eta_LR, double eta_CT, double eta_TD, bool local_print = false);
  void LaunchSM_LR(int params[], std::vector<int> samecl_index, std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, int ngibbs, std::default_random_engine& gen, double eta_CT);
  void ProposalSM_LR_iter(bool SplitMerge, int params[], std::vector<int> samecl_index, std::list<int>& new_c1_ind, std::list<int>& new_c2_ind, std::vector<int>& samecl_newcl, double& lqratio, std::default_random_engine& myRandomEngine, double eta_CT);
  
  std::vector<bool> SampleMix(int ngibbs, double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD, bool sampleLR, bool sampleHR, bool local_print = false);
  bool SampleNew_LR(double temp, std::default_random_engine& eng, double eta_LR, double eta_CT, double eta_TD, bool local_print = false);
  void ProposalNew_LR(double& lqratio, double& lpratio, double& llikratio, proposal_info* info, std::default_random_engine& myRandomEngine, double eta_LR, double eta_CT, double eta_TD, bool local_print = false);
  void save_proposal(proposal_info* info, double eta_LR, double eta_CT, double eta_TD, bool local_print = false);

  // Rcpp::List format_multires();
  // double tables_llikelihoood(int k, LPPartition i_alone);
  double get_tables_merge(int rest, int threshold, std::default_random_engine& myRandomEngine, double eta_CT, double eta_TD);
  double get_tables_merge_reverse(int rest1, int rest2, LPPartition oldpart, std::default_random_engine& myRandomEngine);
  double splitted_tables(int rest, std::list<int> new_c1_ind, std::list<int> new_c2_ind, double eta_CT);
  double splitted_tables(int rest1, int rest2, std::list<int> new_c1_ind, std::list<int> new_c2_ind, double eta_CT);
  double get_logposterior(double eta_LR, double eta_CT, double eta_TD);
};



#endif /* MULTIRES_H_ */