/*
 * various_functions.h
 *
 *  Created on: January 23, 2018
 *      Author: Sameer
 */


#ifndef VARIOUS_FUNCTIONS_H_
#define VARIOUS_FUNCTIONS_H_

#include <stdio.h>
#include <armadillo>
#include <random>

using namespace std;
using namespace arma;

struct proposal_info {
  int jLR;
  int jLR_index;
  int kLR_old;
  int kLR_new;
  std::vector<int> cost_locID_j;    // used
  std::vector<int> cost_locID_NOTj; // used
  std::vector<int> cost_indexing_j;
  std::vector<int> oldtables;       // used
  std::vector<int> oldtables_mainID;
  std::vector<int> olddishes;
  int n_cost_j;
  std::vector<int> remove_table;    // used
  std::vector<int> remove_dish;     // used
  // int n_tables;
  // int K_dishes;
  // std::vector<int> nk_dishes;
  std::vector<int> new_tables;
  // std::vector<int> new_tables_mainID; // you should not use this.
  std::vector<int> new_dishes;
};

arma::mat Submatrix(arma::mat M, int n_rows, int n_cols, int* row_index, int* col_index);
arma::mat Submatrix(arma::mat M, int n_rows, int n_cols, std::vector<int> row_index, int* col_index);
arma::mat Submatrix(arma::mat M, int n_rows, int n_cols, std::vector<int> row_index, std::vector<int> col_index);
void Connected_Components(arma::mat M, int n, int* components, int* count);
// std::vector<std::vector<int> > Alternative_Connected_Components(std::vector<int> remain);
// void Alternative_Connected_Components(int element, std::vector<std::vector<int> >& current_conncomp);
void DFSUtil(arma::mat M, int n, int v, bool* visited, int* components, int* count);
// int* which_is_nearest(arma::mat centroids, arma::mat data);
int* which_is_nearest_k(arma::mat centroids, arma::mat data);
arma::mat Distance_matrix(double* beta_hat, int nObs);
double lbeta(double a, double b);
double lbinomial(int n, int k);
double lmgamma(int p, double x);
arma::mat block_inverse(arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22);
arma::mat block_inverse_ret(arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22, arma::mat* B11, arma::mat* B12, arma::mat* B21, arma::mat* B22);
double block_log_det(arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22);

int Costumer_RestId_MainId(int j, int i, std::vector<int> nj);
int Table_RestId_MainId(int j, int t, std::vector<int> mj);
void Costumer_MainId_RestId(int x, std::vector<int> nj, int& j, int& i);
void Table_MainId_RestId(int y, std::vector<int> mj, int& j, int& t);

std::vector<int> sample_without_replacement(int k, int N, std::default_random_engine& gen, bool shuffle = true);

#endif /* VARIOUS_FUNCTIONS_H_ */

