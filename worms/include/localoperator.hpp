#pragma once
#include <iostream>
#include <stdio.h>
#include <vector>
#include <array>
#include <string>
#include <numeric>
#include <random>
#include <math.h>
#include <bcl.hpp>
#include <assert.h> 
#include <fstream>
#include <tuple>
#include "outgoing_weight.hpp"
#include "funcs.hpp"


namespace model{
// template <class MC = bcl::heatbath>
// class local_operator;

/*
*params
-------
leg : number of sites bond operato acts on. typically 2.
size : number of hilbert space of bond operator.
sps : spin freedom per site.
*template argument
-------
MC : type of algorithm for generating transition matrix

*variables
-------
TPROB : type of transition matrix
*/
typedef std::mt19937 engine_type;
typedef bcl::markov<engine_type> markov_t;
typedef bcl::sparse_markov<engine_type> sparse_markov_t;


struct markov_v{
private:
  std::vector<markov_t> markov;
  std::vector<long long> s2i; 
public:
  markov_v(){}
  markov_v(std::vector<markov_t> markov, std::vector<long long> s2i)
  :markov(markov), s2i(s2i){}
  markov_t & operator[](size_t s) { 
    // auto i = s2i[s];
    auto i = s;
    if (i < 0) std::invalid_argument("probability that state s appears have to be 0"); 
    return markov[i]; 
  }
};

template <class MC=bcl::heatbath>
class local_operator{
public:
  using VECD = std::vector<double>;
  using TPROB = std::vector<VECD>; //type for transition probability. typically, this is 2D matrix with 4 x 4 elements( check notebook for detail definition of this type).
  std::vector<std::vector<double>> ham_;
  std::vector<std::vector<double>> ham; // virtual hamiltonian (or maybe absolute of original hamiltonian)
  std::vector<std::vector<double>> ham_rate; // original hamiltonian
private:
public:
  typedef MC MCT;
  outgoing_weight ogwt;

  const size_t sps;
  const int leg; // leg size.
  const int size; // size of operator (2**leg)
  double ene_shift = 0; //energy shift to ensure that diagonal elements of hamiltonian are non-negative
  double max_diagonal_weight_;
  double total_weights; //sum of diagonal elemtns of ham

  std::vector<int> signs; //list of sign defined via the sign of ham_;
  std::vector<TPROB> trans_prob; //num_configuration x 4 x 4 matrix.
  std::array<int, 2> num2index(int num);
  // markov_v markov;
  std::vector<sparse_markov_t> markov;
  std::vector<size_t> sps_base;

  local_operator(int leg, size_t sps = 2);

  void set_ham(double off_set = 0, double thres = 1E-8, bool dw = false);
  void set_trans_weights();
  void check_trans_prob();
  int index2num(std::array<int, 2> index);
};



template <class MC>
local_operator<MC>::local_operator(int leg, size_t sps)
  :leg(leg), size(pow(sps, leg)), ogwt(leg, sps), sps(sps)
  {
    std::cout << "malloc" << std::endl;
    ham_rate = std::vector<std::vector<double>>(size, std::vector<double>(size, 0));
    ham = std::vector<std::vector<double>>(size, std::vector<double>(size, 0));
  }

/*
setting various variable for local_operators 
this function should be called after manually define local hamiltonian.

*params
------
boolean zw : 1 = have a chance to delete a worm while updating.
*/
template <class MC>
void local_operator<MC>::set_ham(double off_set, double thres, bool zw){
  // std::cout << "Hi" << std::endl;
  int N = size*size;
  ene_shift=0;
  ham_ = ham;

  std::vector<double> ham_vector;
  std::vector<double> ham_rate_vector;

  ham_vector = std::vector<double>(size*size, 0);
  ham_rate_vector = std::vector<double>(size*size, 0);

  for (int i=0; i<ham_.size();i++){
    ene_shift = std::min(ene_shift, ham[i][i]);
    ene_shift = std::min(ene_shift, ham_rate[i][i]);
  }
  ene_shift *= -1;
  ene_shift += off_set;
  for (int i=0; i<ham_.size();i++){
    ham_[i][i] += ene_shift;
    ham_rate[i][i] += ene_shift;
  }

  for (int i=0; i<N; i++){
    auto index = num2index(i);
    ham_vector[i] = ham_[index[0]][index[1]];
    ham_rate_vector[i] = ham_rate[index[0]][index[1]];
  }


  total_weights = 0;
  double tmp=0;
  max_diagonal_weight_ = 0;
  for (int i=0; i<size; i++) {
    tmp += ham_[i][i];
    max_diagonal_weight_ = std::max(max_diagonal_weight_, ham_[i][i]);
  }

  for (int i=0; i<ham_vector.size(); i++){
    auto& x = ham_vector[i];
    auto& y = ham_rate_vector[i];
    signs.push_back(x >= 0 ? 1 : -1);
    x = std::abs(x);
    if (x < thres) x = 0;

    if (y!=0 && x==0){
      std::cerr << "cannot reweighting since support doesn't cover the original matrix" << std::endl;
      std::cerr << "y : " << y << "  x : " << x << std::endl;
      std::terminate();
    }
    if (x!= 0) y = y/x;
  }
  // set transition probability
  // std::cout << ham_vector << std::endl;
  // std::cout << "a" << std::endl;


  // std::vector<markov_t> markov_tmp;
  // std::vector<long long> state2index;

  // state2index.resize(ham_vector.size(), -1);
  for (size_t s=0; s < ham_vector.size(); s++){
    // if (ham_vector[s] == 0 ) continue;
    // state2index[s] = markov_tmp.size();
    std::tuple<std::vector<double>, std::vector<long long>, std::vector<size_t>> sparse_data = ogwt.init_table_sparse(ham_vector, s, zw);
    if (std::get<0>(sparse_data).size() != 0) markov.push_back(
      sparse_markov_t(MC(), std::get<0>(sparse_data), std::get<1>(sparse_data), std::get<2>(sparse_data))
      ); //* sparse markov
    else markov.push_back(sparse_markov_t());
  }

  // markov = markov_v(markov_tmp, state2index);
  

  //* free memories
  // ham.resize(0);
  // ham_rate.resize(0);
}


template <class MC>
std::array<int, 2> local_operator<MC>::num2index(int num){
  ASSERT(num < size*size, "num is invalid");
  std::array<int, 2> index;
  index[0] = num%size;
  index[1] = num/size;
  return index;
}

template <class MC>
int local_operator<MC>::index2num(std::array<int, 2> index){
  ASSERT(index[0] < size && index[1] < size, "index is invalid");
  int num = 0;
  num += index[0];
  num += index[1] * size;
  return num;
}
}