#ifndef __model__
#define __model__
#include <iostream>
#include <stdio.h>
#include <vector>
#include <array>
#include <string>
#include <numeric>
#include <random>
#include <math.h>
#include <bcl.hpp>
#include <lattice/graph.hpp>
#include <lattice/coloring.hpp>
#include <algorithm>
#include "outgoing_weight.hpp"


#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

#define TOR 1.0e-10


#ifdef TOR
  #define DGREATER(X1, X2) (X1 >= X2-TOR)
#else
  #define DGREATER(X1, X2) (x1 >= X2)
  #define TOR 0
#endif

namespace model {

  template <int N_op, size_t nls = 1, size_t max_L = 4, class MC = bcl::st2010>
  class base_spin_model;

  template <class MC = bcl::st2010>
  class local_operator;
  
  using SPIN = unsigned short;
  using STATE = std::vector<SPIN>;
  using BOND = std::vector<std::size_t>;
  inline std::vector<BOND> generate_bonds(lattice::graph lattice){
    std::vector<BOND> bonds;
    for (int b=0; b<lattice.num_bonds(); b++){
      std::vector<size_t> tmp(2);
      tmp[0] = lattice.source(b);
      tmp[1] = lattice.target(b);
      bonds.push_back(tmp);
    }
    return bonds;
  }

  inline std::vector<size_t> generate_bond_type(lattice::graph lattice){
    std::vector<size_t> bond_type;
    for (int b=0; b<lattice.num_bonds(); b++) bond_type.push_back(lattice.bond_type(b));
    return bond_type;
  }

  inline size_t num_type(std::vector<size_t> bond_type){
    std::sort(bond_type.begin(), bond_type.end());
    auto it = std::unique(bond_type.begin(), bond_type.end());
    return std::distance(bond_type.begin(), it);
  }
}



template <class MC>
class model::local_operator{
public:
  using VECD = std::vector<double>;
  using TPROB = std::vector<VECD>; //type for transition probability. typically, this is 2D matrix with 4 x 4 elements( check notebook for detail definition of this type).
  typedef MC MCT;
  int leg; // leg size.
  int size; // size of operator (2**leg)
  size_t nls;
  double ene_shift = 0; //energy shift to ensure that diagonal elements of hamiltonian are non-negative
  std::vector<std::vector<double>> ham;
  std::vector<std::vector<double>> ham_;
  std::vector<double> ham_vector;
  std::vector<int> signs; //list of sign defined via the sign of ham_;
  std::vector<TPROB> trans_prob; //num_configuration x 4 x 4 matrix.
  double max_diagonal_weight_;
  double total_weights; //sum of diagonal elemtns of ham
  std::array<int, 2> num2index(int num);
  // craete markov
  // table for worm update
  typedef std::mt19937 engine_type;
  typedef bcl::markov<engine_type> markov_t;
  outgoing_weight ogwt;
  std::vector<markov_t> markov;


  local_operator(int leg, size_t nls = 1);
  local_operator();

  void set_ham(double off_set = 0);
  void set_trans_weights();
  void check_trans_prob();
  int index2num(std::array<int, 2> index);



};

/*
//$\hat{H} = \sum_{<i,j>} [J \vec{S}_i \dot \vec{S}_j - h/Nb (S_i^z + S_j^z)]$ 
// map spin to binary number e.g. -1 \rightarrow 0, 1 \rightarrow 1
* S is local freedomness. 
*/
template <int N_op, size_t _nls, size_t _max_L, class MC>
class model::base_spin_model{
public:
  const int L;
  static const size_t max_L = _max_L;
  const int Nb; // number of bonds.
  static const int Nop = N_op; //number of local operator (1 for heisenberg model)
  static const size_t nls = _nls;
  double rho = 0;
  std::vector<double> shifts;
  typedef MC MCT;
  std::array<local_operator<MCT>, N_op> loperators; //in case where there are three or more body interactions.
  std::array<int, N_op> leg_size; //size of local operators;
  const std::vector<BOND> bonds;
  const std::vector<size_t> bond_type;
  std::array<size_t, N_op> bond_t_size;
  lattice::graph lattice;
  base_spin_model(int L_, int Nb_, std::vector<BOND> bonds)
  :L(L_), Nb(Nb_), bonds(bonds){}

  base_spin_model(lattice::graph lt)
  :L(lt.num_sites()), Nb(lt.num_bonds()), lattice(lt), 
    bonds(generate_bonds(lt)), bond_type(generate_bond_type(lt))
  {
    int sum = 0;
    for (int i=0; i<Nop; i++){
      bond_t_size[i] = 0;
      for (auto bt : bond_type){
        if (bt==i) bond_t_size[i]++;
      }
      sum += bond_t_size[i];
    }

    if (num_type(bond_type)!=Nop) {
      std::cerr << "Nop is not consistent with number of bond_type" << std::endl;
      std::terminate();
    }

    if (sum != bonds.size()) {
      std::cerr << "something wrong in bond_type" << std::endl;
      std::terminate();
    }
  }
  void initial_setting(std::vector<double>off_sets = std::vector<double>(N_op,0)){
    int i = 0;
    double tmp=0;
    for (auto& x : loperators){
      x.set_ham(off_sets[i]);
      shifts.push_back(x.ene_shift);
      i++;
    }
  }
};





// define functions for lolcal_operator class
template <class MC>
model::local_operator<MC>::local_operator()
  :local_operator(2){}

template <class MC>
model::local_operator<MC>::local_operator(int leg, size_t nls)
  :leg(leg), size(1<<nls * leg), ogwt(leg, nls), nls(nls){

  if (nls<=0) size = (1<<leg); // default size is 2**leg.
  ham = std::vector<std::vector<double>>(size, std::vector<double>(size, 0));
  ham_vector = std::vector<double>(size*size, 0);
}


/*
setting various variable for local_operators 
this function should be called after manually define 2D local hamiltonian.

- set 1D hamiltonian 
*/
template <class MC>
void model::local_operator<MC>::set_ham(double off_set){
  int N = ham_vector.size();
  ene_shift=0;
  ham_ = ham;

  for (int i=0; i<ham_.size();i++){
    ene_shift = std::min(ene_shift, ham[i][i]);
  }
  ene_shift *= -1;
  ene_shift += off_set;
  for (int i=0; i<ham_.size();i++){
    ham_[i][i] = ham_[i][i] + ene_shift;
  }

  for (int i=0; i<N; i++){
    auto index = num2index(i);
    ham_vector[i] = ham_[index[0]][index[1]];
  }


  total_weights = 0;
  // for (int i=0; i<size; i++) total_weights+= ham[i][i];

  double tmp=0;
  max_diagonal_weight_ = 0;
  for (int i=0; i<size; i++) {
    tmp += ham_[i][i];
    max_diagonal_weight_ = std::max(max_diagonal_weight_, ham_[i][i]);
  }



  // max_diagonal_weight_ = std::max(max_diagonal_weight_, weights_[p]);

  for (auto& x : ham_vector){
    signs.push_back(x >= 0 ? 1 : -1);
    x = std::abs(x);
  }

  // set transition probability
  ogwt.init_table(ham_vector);
  for (int c = 0; c < ogwt.size(); ++c) markov.push_back(markov_t(bcl::st2010(),ogwt[c]));

  // auto rand_src = engine_type(2021);
  // auto xxx = markov[0](0, rand_src);



  // check_trans_prob(); // check if transition probability is consistent with the definition of transition matrix

}




template <class MC>
std::array<int, 2> model::local_operator<MC>::num2index(int num){
  ASSERT(num < size*size, "num is invalid");
  std::array<int, 2> index;
  index[0] = num%size;
  index[1] = num/size;
  return index;
}

template <class MC>
int model::local_operator<MC>::index2num(std::array<int, 2> index){
  ASSERT(index[0] < size && index[1] < size, "index is invalid");
  int num = 0;
  num += index[0];
  num += index[1] * size;
  return num;
}




#endif