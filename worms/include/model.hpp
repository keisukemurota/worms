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

  template <int N_op>
  class base_spin_model;
  class local_operator;
  
  using SPIN = unsigned short;
  using STATE = std::vector<SPIN> ;
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

}



class model::local_operator{
public:
  using VECD = std::vector<double>;
  using TPROB = std::vector<VECD>; //type for transition probability. typically, this is 2D matrix with 4 x 4 elements( check notebook for detail definition of this type).


  int L; // number of site operator acts 
  int size; // size of operator (2**L)
  double ene_shift = 0; //energy shift to ensure that diagonal elements of hamiltonian are non-negative
  std::vector<std::vector<double>> ham;
  std::vector<std::vector<double>> ham_;
  std::vector<double> ham_vector;
  std::vector<std::vector<double>> trans_weights;
  std::vector<std::vector<double>> ori_trans_weights; //original weights, which can have negative elements;
  std::vector<int> signs; //list of sign defined via the sign of ham_;
  std::vector<TPROB> trans_prob; //num_configuration x 4 x 4 matrix.
  std::vector<double> diagonal_cum_weight; //normalized diagonal elements;
  std::vector<double> accept; //normalized diagonal elements;
  double max_diagonal_weight_;
  double total_weights; //sum of diagonal elemtns of ham
  std::array<int, 2> num2index(int num);

  // craete markov
  // table for worm update
  typedef std::mt19937 engine_type;
  typedef bcl::markov<engine_type> markov_t;
  outgoing_weight ogwt;
  std::vector<markov_t> markov;


  local_operator(int L);
  local_operator();

  void set_ham();
  void set_trans_weights();
  void set_trans_prob();
  void check_trans_status(VECD, TPROB);
  void check_trans_prob();
  int index2num(std::array<int, 2> index);



  void print_trans_weights(){
    for (int row=0; row<trans_weights.size(); row++){
        for(int column=0; column<trans_weights[0].size(); column++){
          printf("%.2f   ", trans_weights[row][column]);}
        printf("\n");
      }
  }
};

//$\hat{H} = \sum_{<i,j>} [J \vec{S}_i \dot \vec{S}_j - h/Nb (S_i^z + S_j^z)]$ 
// map spin to binary number e.g. -1 \rightarrow 0, 1 \rightarrow 1
template <int N_op = 1>
class model::base_spin_model{
public:
  const int L;
  const int Nb; // number of bonds.
  static const int Nop = N_op; //number of local operator (1 for heisenberg model)
  double rho = 0;
  std::vector<double> shifts;
  std::array<local_operator, N_op> loperators; //in case where there are three or more body interactions.
  std::array<int, N_op> leg_size; //size of local operators;
  std::array<double, N_op> operator_cum_weights;
  const std::vector<BOND> bonds;
  lattice::graph lattice;
  base_spin_model(int L_, int Nb_, std::vector<BOND> bonds)
  :L(L_), Nb(Nb_), bonds(bonds){}

  base_spin_model(lattice::graph lt)
  :L(lt.num_sites()), Nb(lt.num_bonds()), lattice(lt), bonds(generate_bonds(lt))
  {
    std::cerr << "lattice cannot be used for the base_spin_model with N_OP != 1" << std::endl;
  }
};










#endif