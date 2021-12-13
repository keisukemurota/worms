#ifndef __loop__
#define __loop__


#pragma once
#include <string.h>
#include <iostream>
#include <uftree.hpp>
#include <vector>
#include <random>
#include <fstream>
#include <ostream>
#include <strstream>
#include <sstream>
#include <algorithm>

#ifdef __APPLE__
#  include <mach-o/dyld.h>
#endif
#if defined(_WIN32)
#  include <windows.h>
#else
  #include <unistd.h>
#endif


#include <filesystem>
#include <unistd.h>

#include <ctime>
#include <math.h> 

#include "state.hpp"
#include "model.hpp"
#include "BC.hpp"
#define SEED 2030
/* inherit UnionFindTree and add find_and_flip function*/

// template <typename MODEL>

inline int positive_modulo(int i, int n) {
    return (i % n + n) % n;
}
using MODEL = model::heisenberg1D;
using OPS = std::vector<spin_state::OpStatePtr>;
using BSTATE = spin_state::BottomState;
using WORMS = spin_state::Worms;
using DOTS = std::vector<spin_state::Dot>;




class worm{
  public:
  MODEL model;
  double beta;
  int L;
  int W; //number of worms
  std::vector<int> front_dots;// initilize this by f[i] = -(i+1) so that we can check wether a operator is assigned for the ith site already
  std::vector<int> end_dots; //initialize in the same manner.
  OPS ops_tmp1; 
  OPS ops_tmp2; 


  OPS& ops_main = ops_tmp1;  //contains operators.
  OPS& ops_sub = ops_tmp2; // for sub.
  spin_state::BStatePtr pstate;
  spin_state::WormsPtr pworms;

  spin_state::Worms& worms = *pworms;
  spin_state::BottomState& state = *pstate;
  std::vector<int> worms_label;


  // BSTATE& state;
  // WORMS& worms;

  DOTS spacetime_dots; //contain dots in space-time.

  std::vector<double> worms_tau;
  std::vector<std::vector<int>> bonds;

  //declaration for random number generator

  #ifdef RANDOM_SEED
  std::mt19937 rand_src = std::mt19937(static_cast <unsigned> (time(0)));
  #else
  std::mt19937 rand_src = std::mt19937(SEED);
  #endif

  //for choosing bonds
  std::uniform_int_distribution<> dist;

  // for binary dice
  std::uniform_int_distribution<> binary_dice = std::uniform_int_distribution<>(0,1); 

  // random distribution from 0 to beta 
  std::uniform_real_distribution<> worm_dist;


  // reference of member variables from model class

  static const int N_op = MODEL::Nop;

  std::array<model::local_operator, N_op>& loperators; //holds multiple local operators
  std::array<int, N_op>& leg_sizes; //leg size of local operators;
  std::array<double, N_op>& operator_cum_weights;



  worm(double beta, MODEL model_, int W)
  :model(model_), L(model.L), beta(beta), W(W),
  dist(0, model.Nb-1), worm_dist(0.0, beta),
  bonds(model.bonds), 
  pstate(spin_state::BStatePtr(new BSTATE(L))), pworms(spin_state::WormsPtr(new WORMS(W))),
  loperators(model.loperators), leg_sizes(model.leg_size),
  operator_cum_weights(model.operator_cum_weights), worms_tau(W),
  front_dots(L,0), end_dots(L,0)
  {
    #ifdef RANDOM_SEED
    srand(static_cast <unsigned> (time(0)));
    #else
    srand(SEED);
    #endif
    worms_tau.resize(W);
    // pstate = new BSTATE(L);
    // pworms = new WORMS(W);
    
  }


  //* functions for initializing
  /* initialize worms */
  void init_worms_rand(){
    worms_label.resize(0);
    for (int i=0; i<W; i++){
      // pworms->worm_site[i] = dist(rand_src);
      // pworms->tau_list[i] =  worm_dist(rand_src);

      worms.worm_site[i] = dist(rand_src);
      worms.tau_list[i] =  worm_dist(rand_src);
    }
    std::sort(worms.tau_list.begin(), worms.tau_list.end(), std::less<double>());
  }

  /* 
  initialize front and end groups 
  */
  void init_front_n_end(){
    for (int i=0; i<L; i++){
      front_dots[i] = -(i+1);
      end_dots[i] = -(i+1);
    }
  }

  void init_states(){
  for (auto& x : state){
    x = binary_dice(rand_src);
    }
  }

  void init_dots(bool add_state = true){
    spacetime_dots.resize(0);
    if (add_state){
      for(int i=0; i<L; i++){
        set_dots(i, 0, 0, i);
      }
    }
  }

  //swapt main and sub
  void swap_oplist(){
    auto tmp_op = ops_sub;
    ops_sub = ops_main;
    ops_main = tmp_op;
  }


  //main functions

  void diagonal_update(){
    double tau = 0;
    int n_worm = 0;

    //initialization
    init_front_n_end();
    init_worms_rand();
    // init_ops_main()
    ops_main.resize(0);
    //init spacetime_dots
    init_dots();

    auto& worm_site = worms.worm_site;
    auto& worm_tau_list = worms.tau_list;

    double worm_tau = 0;
    if (worm_tau_list.size()) worm_tau = worm_tau_list[0];
    double op_sub_tau = 0;
    if (ops_sub.size()) op_sub_tau = ops_sub[0]->tau;
    int op_sub_label = 0;
    int N_op = model.Nop;


    int optau = 0;
    if (ops_sub.size()) optau = ops_sub[0]->tau;

    double sum;

    int s0, s1;
    int r_bond; // randomly choosen bond
    std::vector<int> cstate = state;





    //set worms
    while (true){
      // cout << "hi" << endl;
      double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
      double tau_prime = tau - log(r)/model.rho;

      // put worms on space.
      while(worm_tau<tau_prime && n_worm < W){
        int site = worm_site[n_worm];
        set_dots(site, worm_tau, 2 , n_worm);
        worms[n_worm] = cstate[site];
        n_worm++;
        worm_tau = worm_tau_list[n_worm];
      }

      checkODNFlip(op_sub_tau, tau_prime, op_sub_label, cstate);

      //choose and insert diagonal operator.

      if (tau_prime > beta) break;

      r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
      double max_ = *(operator_cum_weights.end()-1);
      int lop_label;
      for(lop_label=0; lop_label < N_op; lop_label++){
        if (operator_cum_weights[lop_label] >= r * max_) break;
      }

      int leg_size = leg_sizes[lop_label]; //size of choosen operator
      auto& lop = loperators[lop_label];
      auto& diag_cum_weight = lop.diagonal_cum_weight;
      max_ = lop.total_weights;

      r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
      int s_num; //choose local state 
      for (s_num=0; s_num < (1<<leg_size); s_num++){
        if (diag_cum_weight[s_num] >= r * max_) break;
      }

      // choose bond
      r_bond = dist(rand_src);
      int tuggle = 1;
      auto local_state = spin_state::num2state(s_num + (s_num<<leg_size ), 2*leg_size);
      std::vector<int> labels(leg_size);
      auto bond = bonds[r_bond];

      int n_dots = spacetime_dots.size();
      for (int i=0; i<leg_size; i++){
        labels[i] = n_dots;
        n_dots++;
        int s = bond[i];
        if (cstate[s] != local_state[i]) tuggle = 0;
      }



      if ( tuggle ){
        ops_main.emplace_back(
            new spin_state::OpState(
              local_state,
              &lop,
              bond,
              tau_prime
          )
        );

        int dot_label = spacetime_dots.size();
        int n = ops_main.size();
        for (int i=0; i<leg_size; i++){
          set_dots(bond[i], tau_prime, 1 , i);
          dot_label++;
        }
      }
      
      tau = tau_prime;
    } //end of while loop

  }

  void checkODNFlip(double& optau, double tau_prime, int& op_label,
                     std::vector<int>& cstate ){
    if (op_label < ops_sub.size()) optau = ops_sub[op_label]->tau;
    while(optau<tau_prime && op_label < ops_sub.size()){
      auto op_ptr = ops_sub[op_label];
      if (op_ptr->is_off_diagonal()){
        update_state(op_ptr, cstate);
      }

      ops_main.push_back(op_ptr);
      int dot_label = spacetime_dots.size();
      for (int i=0; i<op_ptr->L; i++){
        set_dots(op_ptr->bond[i], op_ptr->tau, 1 , i);
        dot_label++;
      }
      op_label++;
      if (op_label < ops_sub.size()) optau = ops_sub[op_label]->tau;
    }

    return;
  }

  void worm_update();

  /*
  this function will be called after assigining op_main
  */
  void set_dots(int site, double tau_prime, int dot_type, int index){

    int* sptr;
    spin_state::BaseStatePtr stateptr;

    // ASSERT(label == spacetime_dots.size()+1)
    int label = spacetime_dots.size();

    if (dot_type == 0) {
      sptr = state.data() + index;
      stateptr = pstate;
    }else if(dot_type == 1){
      int n = ops_main.size();
      sptr = ops_main[n-1]->data() + index;
      stateptr = ops_main[n-1];
    }else if(dot_type == 2){
      sptr = worms.data() + index;
      stateptr = pworms;
      worms_label.push_back(label);
    }

    if (end_dots[site] < 0){
      end_dots[site] = label;
    } 
    int prev = front_dots[site];
    int end = end_dots[site];


    spacetime_dots.emplace_back(
      site, tau_prime, prev, end, sptr,
      stateptr, dot_type
    );

    if (prev>=0) spacetime_dots[prev].set_next(label);
    
    if (end>=0) spacetime_dots[end].set_prev(label);

    front_dots[site] = label;

  }

  /*
  update given state by given operator ptr;
  */
  static void update_state(spin_state::OpStatePtr op_ptr, std::vector<int>& state){
    std::vector<int> local_state = *op_ptr;
    std::vector<int> state_(op_ptr->plop->L);
    int i=0;
    for (auto x : op_ptr->bond){
      state_[i] = state[x];
      i++;
    }
    ASSERT(is_same_state(local_state, state_), "the operator can not be applied to the state");
    if (op_ptr->is_off_diagonal()){
      int index = 0;
      for (auto x : op_ptr->bond){
        state[x] = local_state[state_.size() + index];
        index++;
      }
    }
  }


  //* functions for testing
  static void check_operators(spin_state::BottomState state, OPS ops){
    const auto state_ = state;
    for (auto op : ops){
      update_state(op, state);
    }
    ASSERT(is_same_state(state_, state), "operators are not consistent");
  }

  static bool is_same_state(int n, int m){
    return n==m;
  }

  static bool is_same_state(int n, std::vector<int> state){
    int m = spin_state::state2num(state, state.size());
    return n==m;
  }

  static bool is_same_state( std::vector<int> state_, std::vector<int> state){
    int m = spin_state::state2num(state, state.size());
    int n= spin_state::state2num(state_, state.size());
    return n==m;
  }


};

#endif 