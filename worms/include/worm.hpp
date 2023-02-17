#pragma once
#include <string.h>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <ostream>
#include <strstream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <unordered_set>
#include <bcl.hpp>
#include <stdlib.h>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <filesystem>
#include <unistd.h>

#include <ctime>
#include <math.h>

#include "state2.hpp"
#include "operator.hpp"
#include "automodel.hpp"
#include "funcs.hpp"
#include "autoobservable.hpp"
#define SEED 1662509963
/* inherit UnionFindTree and add find_and_flip function*/

// template <typename MODEL>

inline int positive_modulo(int i, int n)
{
  return (i % n + n) % n;
}

using spin_state::Dotv2;

// using MODEL = model::heisenberg1D;
using state_t = spin_state::VUS;
using SPIN = spin_state::US;
using BOND = model::VS;
using WORMS = spin_state::WORM_ARR;
using DOTS = std::vector<Dotv2>;
using size_t = std::size_t;

template <class MCT>
class Worm
{

public:
  using MODEL = model::base_model<MCT>;
  using LOPt = model::local_operator<MCT>;

private:
  model::MapWormObs _mp_worm_obs;
  alps::alea::batch_acc<double> _phys_cnt;

  // define observables
  // n*  number of physically meaningful configurations;
  double phys_cnt;
  // n*  sum of observables encountered while worm update. (observable must be non-diagonal operator)
  vector<double> obs_sum;

public:
  typedef spin_state::Operator OP_type;
  typedef std::vector<OP_type> OPS;
  typedef spin_state::StateFunc state_func;
  typedef std::mt19937 engine_type;
  typedef std::uniform_real_distribution<> uniform_t;
  typedef std::exponential_distribution<> expdist_t;
  typedef bcl::markov<engine_type> markov_t;

  uniform_t uniform;

  engine_type rand_src;
  engine_type test_src;

  MODEL spin_model;
  OPS ops_main; // n* contains operators.
  OPS ops_sub;  // n*  for sub.
  state_t state;
  state_t cstate;
  DOTS spacetime_dots; // n*  contain dots in space-time.
  WORMS worms_list;

  VVS pows_vec;
  size_t sps;
  std::vector<BOND> bonds;
  std::vector<size_t> bond_type;
  std::vector<state_func> state_funcs;
  std::unordered_set<size_t> can_warp_ops;
  std::vector<size_t> pres = std::vector<size_t>(0);
  std::vector<size_t> psop = std::vector<size_t>(0);
  std::vector<state_t> worm_states;
  std::vector<double> worm_taus;

  // n* reference of member variables from model class
  std::vector<model::local_operator<typename MODEL::MCT>> &loperators;
  std::vector<std::vector<double>> accepts; // n* normalized diagonal elements;

  // end of define observables

  int sign = 1;
  int cnt = 0;
  const int L; // n* number of sites
  const int N_op;
  size_t d_cnt = 0;
  size_t bocnt = 0;
  const size_t cutoff_length; // n* cut_off length
  size_t u_cnt = 0;
  double rho;
  const double beta;

  std::unordered_map<std::string, WormObs> &get_worm_obs() { return _mp_worm_obs(); }
  alps::alea::batch_acc<double> &get_phys_cnt() { return _phys_cnt; }

  Worm(double beta, MODEL model_, size_t cl = SIZE_MAX, int rank = 0)
      : Worm(beta, model_, model::WormObs(model_.sps_sites(0)), cl, rank) {}

  Worm(double beta, MODEL model_, model::MapWormObs mp_worm_obs_, size_t cl = SIZE_MAX, int rank = 0)
      : spin_model(model_), L(spin_model.L), beta(beta), rho(-1), N_op(spin_model.N_op),
        bonds(spin_model.bonds), bond_type(spin_model.bond_type), state(spin_model.L), cstate(spin_model.L), cutoff_length(cl), _mp_worm_obs(mp_worm_obs_), _phys_cnt(1),
        loperators(spin_model.loperators), sps(spin_model.sps_sites(0)), obs_sum(mp_worm_obs_().size())
  {

    srand(rank);
#ifdef NDEBUG
    unsigned rseed = static_cast<unsigned>(time(0)) + rand() * (rank + 1);
    rand_src = engine_type(SEED);
#else
    rand_src = engine_type(SEED);
    test_src = engine_type(SEED);
#endif
    double max_diagonal_weight = loperators[0].max_diagonal_weight_;
    for (auto const &lop : loperators)
    {
      max_diagonal_weight = std::max(max_diagonal_weight, lop.max_diagonal_weight_);
    }
    for (int i = 0; i < loperators.size(); i++)
    {
      LOPt const &lop = loperators[i];
      pows_vec.push_back(lop.ogwt.pows);
      state_funcs.push_back({lop.ogwt.pows, lop.ogwt.L});
      auto accept = std::vector<double>(lop.size, 0);

      auto const &ham = lop.ham_prime;
      for (int j = 0; j < lop.size; j++)
      {
        accept[j] = ham[j][j] / max_diagonal_weight;
      }
      accepts.push_back(accept);
      rho = max_diagonal_weight * spin_model.Nb;
    }
  }

  void initStates()
  { //* initialized to all up
    for (int i = 0; i < state.size(); i++)
    {
#ifdef RANDOM_SEED
      double r = uniform(rand_src);
      state[i] = static_cast<SPIN>(sps * r);
#else
      state[i] = 0;
#endif
    }
  }

  void initDots(bool add_state = true)
  {
    spacetime_dots.resize(0);
    if (add_state)
    {
      for (int i = 0; i < L; i++)
      {
        set_dots(i, -1, i);
      }
    }
  }

  // swapt main and sub
  void swapOps()
  {
    std::swap(ops_main, ops_sub);
  }

  // main functions

  void diagonalUpdate(double wdensity)
  {
    // dout << "random : " <<  uniform(rand_src) << endl;

    swapOps();
    // wdensity = 3;

    expdist_t expdist(rho * beta + wdensity);           // initialize exponential distribution
    double pstart = wdensity / (beta * rho + wdensity); // probability of choosing worms
    std::copy(state.begin(), state.end(), cstate.begin());
    size_t lop_label;
    lop_label = 0; // typically, lop_label is fixed to 0
    // int leg_size = leg_sizes[lop_label]; //size of choosen operator

    ops_main.resize(0);   //* init_ops_main()
    can_warp_ops.clear(); // forget previous record

    initDots(); //*init spacetime_dots

    // n*init worms
    worms_list.resize(0);
    worm_states.resize(0);
    worm_taus.resize(0);

    ops_sub.push_back(OP_type::sentinel(1)); //*sentinels
    double tau0 = uniform(rand_src);
    double tau = expdist(rand_src);
    bool set_atleast_one = false;
    for (typename OPS::iterator opi = ops_sub.begin(); opi != ops_sub.end();)
    {
      if (tau0 < tau && !set_atleast_one)
      {
        if (tau0 < opi->tau())
        {
          size_t s = static_cast<int>(L * uniform(rand_src));
          appendWorms(worms_list, s, spacetime_dots.size(), tau0);
          set_dots(s, -2, 0); //*index is always 0
          set_atleast_one = true;
          worm_states.push_back(cstate);
          worm_taus.push_back(tau0);
        }
      }
      // auto op_sub = *opi;
      if (tau < opi->tau())
      { //* if new point is behind the next operator is opsub.
        double r = uniform(rand_src);
        if (r < pstart)
        {
          size_t s = static_cast<int>(L * uniform(rand_src));
          appendWorms(worms_list, s, spacetime_dots.size(), tau);
          set_dots(s, -2, 0); //*index is always 0
          worm_states.push_back(cstate);
          worm_taus.push_back(tau);
        }
        else
        {
          size_t b = static_cast<size_t>(bonds.size() * uniform(rand_src));
          lop_label = bond_type[b];
          auto const &accept = accepts[lop_label];
          auto const &bond = bonds[b];
          size_t u = state_funcs[lop_label].state2num(cstate, bond);
          // printd("u = %lu\t", u);
          // printd("input = %lu\t", u);
          // printd("opi_type = %lu\n", opi->op_type());
          // dout << endl;
          // size_t u = spin_state::state2num(cstate, bond);

          r = uniform(rand_src);

          if (r < accept[u])
          {
            // dout << "append r " << r << endl;
            appendOps(ops_main, spacetime_dots, can_warp_ops,
                      &bond, &pows_vec[lop_label], u * pows_vec[lop_label][bond.size()] + u, lop_label, tau);
          }
        }
        tau += expdist(rand_src);
      }
      else
      { //*if tau went over the operator time.
        if (opi->is_off_diagonal())
        {

          update_state(opi, cstate);
          appendOps(ops_main, spacetime_dots, can_warp_ops,
                    opi->bond_ptr(), opi->pows_ptr(), opi->state(), opi->op_type(), opi->tau());
          // printStateAtTime(cstate, tau);
        }
        ++opi;
      }
    } // end of while loop

    //* comment out
    // #ifndef NDEBUG
    // for (typename OPS::iterator opi = ops_main.begin(); opi != ops_main.end();++opi){
    //   printf("[%d, %d]\n", opi->bond(0), opi->bond(1));
    // }

    // if (cstate != state){
    //   throw std::runtime_error("diagonalUpdate : state is not updated correctly");
    // }
    // #endif
  }

  /*
  *update Worm for W times.

  variables
  ---------
  dir : direction of worm head. 1 : upward, -1 : downward
  */
  void wormUpdate(double &wcount, double &wlength)
  {
    pres.resize(0);
    psop.resize(0);
    dout << "\nStart worm update" << endl;
    dout << "---------------------" << endl;
    std::copy(state.begin(), state.end(), cstate.begin());
    ops_main.push_back(OP_type::sentinel(1)); //*sentinels
    typename WORMS::iterator wsi = worms_list.begin();
    size_t w_index = 0;
    for (typename OPS::iterator opi = ops_main.begin(); opi != ops_main.end();)
    {
      if (d_cnt == 132)
      {
        cout << "debug" << endl;
      }
      if (opi->tau() < std::get<2>(*wsi))
      { // n* if operator is behind the worm tail.
        ++opi;
      }
      else
      {
        // t : tail, h : head. direction of tails is opposite to the direction of the initial head.
        // prime means the spin state in front of the worm.
        size_t wt_dot, site, wt_site, _t_spin, n_dot_label;
        double wt_tau, tau;                       // tau_prime is the time of the next operator.
        std::tie(wt_site, wt_dot, wt_tau) = *wsi; // contains site, dot label, tau
        tau = wt_tau;
        site = wt_site;
        Dotv2 *dot = &spacetime_dots[wt_dot];
        Dotv2 *_dot;
        double r = uniform(rand_src);
        size_t dir = (size_t)2 * r, t_dir = !dir;
        // size_t fl = 1;
        int fl = static_cast<int>((sps - 1) * uniform(rand_src)) + 1, t_fl = fl;
        int wl = wlength;
        int br = 0;
        size_t w_x = (worm_states[w_index][wt_site] + fl) % sps; // n* spin state at worm head
        size_t wt_x = w_x;
        bool wh = true; //* Worm head still exists.
        double wlength_prime = 0;
        wcount += 1;

        // if (cstate != worm_states[w_index])
        // {
        //   cout << "cstate != worm_states" << endl;
        //   exit(1);
        // }

        // n* initialize wobs variables
        phys_cnt = 0;
        obs_sum.assign(obs_sum.size(), 0);

        if (d_cnt == 385)
        {
          int x;
          cout << d_cnt << endl;
        }

        do
        {
          n_dot_label = dot->move_next(dir); // next label of dot.

          size_t status = wormOpUpdate(n_dot_label, dir, site, wlength_prime,
                                       fl, tau, wt_dot, wt_site, wt_tau, w_x, wt_x, t_fl, t_dir, w_index);
          if (status != 0)
          {
            if (status == 1)
            {
              std::cerr << "not implemented yet" << endl;
              exit(1);
              // wlength_prime = 0;
              // reset_ops();
              // std::copy(cstate.begin(), cstate.end(), state.begin());
              // br = 1;
            }
          }
          dot = &spacetime_dots[n_dot_label];
        } while ((n_dot_label != wt_dot || ((t_dir == dir ? 1 : -1) * t_fl + fl + sps) % sps != 0));
        _phys_cnt << phys_cnt;
        int obs_i = 0;
        for (auto &obs : _mp_worm_obs())
        {
          obs.second << obs_sum[obs_i];
          obs_i++;
        }

        if (br == 1)
        {
          bocnt++;
          break;
        }

        wlength += wlength_prime;
        checkOpsInUpdate(wt_dot, dir ? n_dot_label : dot->prev(), t_dir, t_fl, fl, dir);
        ++wsi;
        ++w_index;
      }
      if (wsi == worms_list.end())
        break;
    }
#ifndef NDEBUG
// if (cstate != state){
//   throw std::runtime_error("wormUpdate : state is not updated correctly");
// }
#endif
    ops_main.resize(ops_main.size() - 1);
  }

  /*
  This function will be called ever time the head of the worm cross the same propagation level.
  calculate $\langle x_{\tau} | \hat{o}_i \hat{o}_j |x^\prime_{\tau} \rangle$ and update the state of the worm. $x is lower states and x^\prime is upper states$

  note that this function is only called when the worm position of head and tail is apart.
  params
  ------
  tau : imaginary time of the worm head and tail.
  h_site : site of the worm head.
  t_site : site of the worm tail.
  h_x : lower state of the worm head.
  h_x_prime : upper state of the worm head.
  t_w : lower state of the worm tail.
  t_x_prime : upper state of the worm tail.
  */

  void calcHorizontalGreen(double tau, size_t h_site, size_t t_site,
                           size_t h_x, size_t h_x_prime, size_t t_x, size_t t_x_prime,
                           const state_t &_cstate)
  {
    // model::WormObs& _worm_obs = _mp_worm_obs().begin()->second;
    // double _add_phys = 0;
    // if (t_x == t_x_prime) cout << "h_x : " << h_x << " h_x_prime : " << h_x_prime << " t_x : " << t_x << " t_x_prime : " << t_x_prime << endl;
    int j = 0;
    for (auto &obs : _mp_worm_obs())
    {
      auto &_worm_obs = obs.second;
      if (h_site != t_site)
      {
        obs_sum[j] += _worm_obs.second()->operator()(t_x, h_x, t_x_prime, h_x_prime) * L * sign / 2.0;
      }
      else
      {
        if (t_x == t_x_prime)
        {
          // int i = uniform(rand_src) * L;
          for (size_t i = 0; i < L; i++)
          {
            size_t h_x = _cstate[i];
            if (i == t_site)
              obs_sum[j] += _worm_obs.first()->operator()(t_x, t_x) * L * sign / (double)(sps - 1);
            else
              obs_sum[j] += _worm_obs.second()->operator()(t_x, h_x, t_x, h_x) * L / 2.0 * sign / (double)(sps - 1);
          }
          phys_cnt = (double)sign / (sps - 1);
        }
        else
        { // n* This case could contribute to single flip operator but not implemented yet.
          ;
        }
      }
      j++;
    }
  }

  /*
  This function will be called ever time worm head warps.
  */
  void calcWarpGreen(double tau, size_t t_site, size_t t_x, size_t t_x_prime, const state_t &_cstate)
  {
    int j = 0;
    for (auto &obs : _mp_worm_obs())
    {
      auto &_worm_obs = obs.second;
      if (t_x == t_x_prime)
      {
        throw std::runtime_error("t_x == t_x_prime while wapr should never happen");
      }
      // int i = uniform(rand_src) * L;
      for (size_t i = 0; i < L; i++)
      {
        size_t h_x = _cstate[i];
        if (i == t_site)
          obs_sum[j] += _worm_obs.first()->operator()(t_x, t_x_prime) * L * sign * 2;
        else
        {
          obs_sum[j] += _worm_obs.second()->operator()(t_x, h_x, t_x_prime, h_x) * L * sign * 2;
        }
      }
      j++;
    }
  }

  // //*append to ops
  void appendOps(
      OPS &ops,
      DOTS &sp,
      std::unordered_set<size_t> &warp_sp,
      const BOND *const bp,
      const BOND *const pp,
      int state,
      int op_type,
      double tau)
  {

    int s = bp->size();
    ops.push_back(OP_type(bp, pp, state, op_type, tau));

    size_t n = ops.size();
    size_t label = sp.size();
    int site;
    if (loperators[op_type].has_warp(state))
      warp_sp.insert(sp.size()); // if the operator has warp, add the label of the leftmost dot to the set.
    for (int i = 0; i < s; i++)
    {
      // set_dots(bond[i], 0, i);
      site = bp->operator[](i);
      sp.push_back(Dotv2(sp[site].prev(), site, n - 1, i, site));
      sp[sp[site].prev()].set_next(label);
      sp[site].set_prev(label);
      label += 1;
    }
  }
  //* get dot state
  /*
  params
  ------
  ndot_label : label of the dot worm is directing to.
  dir : direction of the worm moves to dot. 1 : upwards, 0 : downwards. So if dir = 1, it means worm comes from below the dot.
  */
  inline size_t getDotState(size_t ndot_label, size_t dir)
  {
    Dotv2 *ndot = &spacetime_dots[ndot_label];
    if (ndot->at_origin())
    {
      return state[ndot->label()];
    }
    else if (ndot->at_operator())
    {
      OP_type &opstate = ops_main[ndot->label()];
      size_t cindex = ndot->leg(!dir, opstate.size()); // direction must be reversed here.
      return opstate.get_local_state(cindex);
    }
    else if (ndot->at_worm())
    {
      return getDotState(ndot->move_next(dir), dir);
    }
    else
    {
      throw std::invalid_argument("dots contains invalid dot type");
      return 0;
    }
  }

  //*append to worms
  inline void appendWorms(WORMS &wm, size_t site, size_t dot_label, double tau)
  {
    wm.push_back(std::make_tuple(site, dot_label, tau));
  }

  /*
  *perform one step from given Worm.
  If dot is operator then, Worm move to exit of the operator. otherwise just assigin spin to dot.
  params
  ------
  int next_dot : next dot.
  int dir : direction Worm is moving toward. 1 : move up, 0 : move down.
  int spin : current spin state of Worm.
  int site : site Worm is at.

  params(member variables)
  ------
  */
  int wormOpUpdate(size_t &next_dot, size_t &dir,
                   size_t &site, double &wlength, int &fl, double &tau,
                   const size_t wt_dot, const size_t wt_site, const double wt_tau,
                   size_t &w_x, size_t &wt_x, const int t_fl, const int t_dir,
                   const size_t w_index)
  {

    OP_type *opsp;
    double tau_prime;
    Dotv2 *dotp = &spacetime_dots[next_dot];

    // get imaginary temperature at the next dot.
    if (dotp->at_origin())
    {
      tau_prime = 0;
      dout << "at origin" << endl;
    }
    else if (dotp->at_operator())
    {
      opsp = &ops_main[dotp->label()];
      tau_prime = opsp->tau();
      dout << "at ops" << endl;
    }
    else if (dotp->at_worm())
    {
      tau_prime = std::get<2>(worms_list[dotp->label()]);
      dout << "at worm" << endl;
    }
    else
    {
      throw std::runtime_error("dot is not at any of the three places");
    }

    // dout << tau << " " << tau_prime << " " << wt_tau << " " << dir << " passed ? " << detectWormCross(tau, tau_prime, wt_tau, dir) << endl;
    dout << "fl :" << fl << " tau : " << tau << " tau_prime : " << tau_prime << " wt_tau : " << wt_tau << " dir : " << dir << " // passed ? " << detectWormCross(tau, tau_prime, wt_tau, dir);
    dout << "  site : [" << site << " " << wt_site << "] " << endl;

    size_t h_x, h_x_prime, t_x, t_x_prime;
    Dotv2 *wtdot = &spacetime_dots[wt_dot];
    // getSpinsDot(next_dot, dotp, dir, h_x, h_x_prime);
    for (size_t i = 0; i < worm_taus.size(); i++)
    {
      if (detectWormCross(tau, tau_prime, worm_taus[i], dir))
      {
        if (i == w_index)
        {
          if (abs(worm_taus[i] - wt_tau) > 1E-10)
            throw std::runtime_error("worms are crossing each other");
          if (wt_site == site)
          {
            int _fl = (dir == t_dir ? fl + t_fl : fl - t_fl) % sps;
            h_x = dir == 1 ? w_x : (w_x - _fl) % sps;
            h_x_prime = dir == 0 ? w_x : (w_x - _fl) % sps;
            t_x = h_x;
            t_x_prime = h_x_prime;
          }
          else
          {
            h_x = dir == 1 ? w_x : (w_x - fl) % sps;
            h_x_prime = dir == 0 ? w_x : (w_x - fl) % sps;
            t_x = t_dir == 1 ? wt_x : (wt_x - t_fl) % sps;
            t_x_prime = t_dir == 0 ? wt_x : (wt_x - t_fl) % sps;
          }

#ifndef NDEBUG
          size_t _h_x, _h_x_prime, _t_x, _t_x_prime;
          _t_x = getDotState(wtdot->move_next(0), 0);
          _t_x_prime = getDotState(wtdot->move_next(1), 1);
          if (dir == 1)
          {
            _h_x_prime = getDotState(next_dot, 1);
            _h_x = getDotState(dotp->move_next(0), 0);
          }
          else if (dir == 0)
          {
            _h_x_prime = getDotState(dotp->move_next(1), 1);
            _h_x = getDotState(next_dot, 0);
          }
          else
          {
            throw std::runtime_error("dir should be either 1 or 0");
          }

          if (h_x != _h_x || h_x_prime != _h_x_prime || t_x != _t_x || t_x_prime != _t_x_prime)
          {
            dout << "t_fl : " << t_fl << " fl : " << fl << " dir : " << dir << " t_dir : " << t_dir << endl;
            dout << "h_x : " << h_x << " " << _h_x << endl;
            dout << "h_x_prime : " << h_x_prime << " " << _h_x_prime << endl;
            throw std::runtime_error("spin is not consistent");
          }
#endif
          calcHorizontalGreen(tau, site, wt_site, h_x, h_x_prime, t_x, t_x_prime, worm_states[w_index]);
          if (wt_site == site)
            wt_x = (wt_x + fl) % sps; // n* fliped by worm head.
        }
        worm_states[i][site] = w_x; // n* assign new spin.
      }
    }

    if (detectWormCross(tau, tau_prime, wt_tau, dir))
    {
      // get spin over and under the worm head.
      // if head is on tail, the case is not regarded as 2point correlation.
    }
    size_t cur_dot = next_dot;
    // ASSERT(site == dotp->site(), "site is not consistent");
    if (dotp->at_origin())
    {
      state[dotp->label()] = (state[dotp->label()] + fl) % sps;
      if (dir || tau == 0)
        wlength += 1 - tau;
      else
        wlength += tau;
    }
    else if (dotp->at_operator())
    {
      dout << "update cnt : " << u_cnt << endl;
      u_cnt++;

      if (opsp->cnt() == 0)
      {
        psop.push_back(dotp->label());
        pres.push_back(opsp->state());
      }
      opsp->add_cnt();

      if (opsp->cnt() > cutoff_length)
      {
        return 1;
      }

      size_t state_num, tmp, nindex;
      size_t dir_in = !dir; // n* direction the Worm comes in from the view of operator.
      size_t leg_size = opsp->size();
      size_t cindex = dotp->leg(dir_in, leg_size);
      size_t index = dotp->leg(0, leg_size);
      size_t op_type = opsp->op_type();

      auto &lop = loperators[op_type];

      sign *= lop.signs[opsp->state()];
      state_num = opsp->update_state(cindex, fl);
      tmp = lop.markov[state_num](cindex * (sps - 1) + sps - fl, rand_src);

      int tmp_wlength = opsp->tau() - tau;
      if (!(dir ^ (tmp_wlength > 0)) & (tmp_wlength != 0))
        wlength += std::abs(tmp_wlength);
      else
        wlength += (1 - std::abs(tmp_wlength));

      //* if Worm stop at this operator
      if (tmp == 0)
      { // if head will warp.
        if (!lop.has_warp(state_num))
          can_warp_ops.erase(cur_dot - dotp->index());
        else
          can_warp_ops.insert(cur_dot - dotp->index());
        sign *= lop.signs[state_num]; // include an effect by start point

        // n* head wapred
        size_t _cur_dot, _state_num, _optype;
        Dotv2 *_dotp;
        std::unordered_set<size_t>::iterator it;
        do
        {
          t_x_prime = getDotState(wtdot->move_next(1), 1);
          t_x = getDotState(wtdot->move_next(0), 0);
          size_t n = static_cast<size_t>(can_warp_ops.size() * uniform(rand_src));
          it = std::begin(can_warp_ops);
          std::advance(it, n);
          if (it == can_warp_ops.begin())
            calcWarpGreen(tau, wt_site, t_x, t_x_prime, worm_states[w_index]);
          _cur_dot = *it;
          _dotp = &spacetime_dots[_cur_dot];
          opsp = &ops_main[_dotp->label()];
          _state_num = opsp->state();
          _optype = opsp->op_type();
          tmp = loperators[_optype].markov[_state_num](0, rand_src);
          leg_size = opsp->size();
        } while (tmp == 0);
        if (tmp == 0)
        {
          std::cerr << "cannot handle warp reject" << endl;
          exit(1);
        }

        tmp--;
        auto &_lop = loperators[_optype];
        sign *= _lop.signs[_state_num]; // warped point
        nindex = tmp / (sps - 1);       // sps_sites are filled with same value.
        fl = tmp % (sps - 1) + 1;
        sign *= _lop.signs[opsp->update_state(nindex, fl)]; // after procceeding
        dir = nindex / (leg_size);
        site = opsp->bond(nindex % leg_size);
        next_dot = opsp->next_dot(0, nindex, _cur_dot);
        if (!_lop.has_warp(_state_num))
          can_warp_ops.erase(_cur_dot - _dotp->index());
        else
          can_warp_ops.insert(_cur_dot - _dotp->index());
        tau_prime = opsp->tau();

        // n* contribute to worm operator which flip single spin.
        //  calcHorizontalGreen(tau, last_site, wt_site, last_state, last_state, t_x, t_x_prime);
        //  calcHorizontalGreen(tau, last_site, wt_site, last_state, last_state,
        //          t_x, t_x_prime, (L - 1) / (double)can_warp_ops.size());
      }
      else
      {
        tmp--;
        nindex = tmp / (sps - 1);
        fl = tmp % (sps - 1) + 1;
        state_num = opsp->update_state(nindex, fl);
        sign *= lop.signs[state_num];
        // n* assigin for next step
        dir = nindex / (leg_size);
        site = opsp->bond(nindex % leg_size);
        next_dot = opsp->next_dot(cindex, nindex, cur_dot);
        if (!lop.has_warp(state_num))
          can_warp_ops.erase(cur_dot - dotp->index());
        else
          can_warp_ops.insert(cur_dot - dotp->index());
      }

      w_x = opsp->get_local_state(nindex);

#ifndef NDEBUG
      int dsign = 1;
      for (auto &op : ops_main)
      {
        dsign *= loperators[op.op_type()].signs[op.state()];
      }
      if (dsign != sign)
      {
        std::cerr << "sign is wrong" << endl;
        exit(1);
      }
// int niter = 0;
// for (int i=0; i<niter; i++){
//   int tmp_ = loperators[opsp->op_type()].markov[state_num](cindex*(sps - 1) + sps-fl-1, test_src);
//   int nindex_ = tmp_/sps - 1;
//   int fl_ = tmp_ % sps - 1 + 1;
//   // printf("test tmp : %d, state : %d\n", tmp_, num ^ (fl_ << (nls*nindex_)));
// }
#endif
    }

    tau = tau_prime;
    return 0;
  }

  /*
   *this function will be called after assigining op_main
   */
  void set_dots(size_t site, size_t dot_type, size_t index)
  {

    size_t label = spacetime_dots.size();

    if (dot_type == -1)
    {
      ASSERT(label == site, "label must be equal to site");
      spacetime_dots.push_back(
          Dotv2::state(site));
    }
    else if (dot_type == -2)
    {
      size_t n = worms_list.size();
      spacetime_dots.push_back(
          Dotv2::worm(spacetime_dots[site].prev(), site, n - 1, site));
      spacetime_dots[spacetime_dots[site].prev()].set_next(label);
      spacetime_dots[site].set_prev(label);
    }
    else
    {
      std::cout << "dot_type must be either -1 or -2" << std::endl;
    }
  }

  void getSpinsDot(size_t next_dot, Dotv2 *dotp, int dir, size_t &h_x, size_t &h_x_prime)
  {
    if (dir == 1)
    {
      h_x_prime = getDotState(next_dot, 1);
      h_x = getDotState(dotp->move_next(0), 0);
    }
    else if (dir == 0)
    {
      h_x_prime = getDotState(dotp->move_next(1), 1);
      h_x = getDotState(next_dot, 0);
    }
  }

  // /*
  // *this function will be called after assigining op_main
  // */
  // void set_op_dots(size_t site, size_t index){
  //   size_t label = spacetime_dots.size();
  //   size_t n = ops_main.size();
  //   spacetime_dots.push_back(
  //     Dotv2(spacetime_dots[site].prev(), site, n-1,index, site)
  //   );
  //   spacetime_dots[spacetime_dots[site].prev()].set_next(label);
  //   spacetime_dots[site].set_prev(label);
  // }

  /*
   *update given state by given operator ptr;
   */
  void update_state(typename OPS::iterator opi, state_t &state)
  {
#ifndef NDEBUG
    state_t local_state = opi->get_state_vec();
    state_t state_(opi->size());
    int i = 0;
    for (auto x : *(opi->bond_ptr()))
    {
      state_[i] = state[x];
      i++;
    }
    ASSERT(is_same_state(local_state, state_, 0), "the operator can not be applied to the state");
#endif
    if (opi->is_off_diagonal())
      update_state_OD(opi, state);
  }

  /*
   *update given state by given offdiagonal operator ptr;
   */
  void update_state_OD(typename OPS::iterator opi, state_t &state)
  {
    int index = 0;
    auto const &bond = *(opi->bond_ptr());
    for (auto x : bond)
    {
      state[x] = opi->get_local_state(bond.size() + index);
      index++;
    }
  }

  /*
  * check the operator and state is consistent during the worm_updateg
  params
  ------
  worm_label : label of dot (Worm) we are aiming at.
  p_label : label of dot before reaching at the current position of Worm.

  */
  void checkOpsInUpdate(int worm_label, int p_label, int t_dir, int t_fl, int fl, int dir)
  {

#ifndef NDEBUG
    auto state_ = state;

    int label = 0;
    d_cnt++;
    std::cout << "debug cnt = " << d_cnt << std::endl;
    for (const auto &dot : spacetime_dots)
    {
      if (dot.at_operator() && (dot.leg(0, 0) == 0))
      { // if dot_type is operator
        auto opi = ops_main.begin() + dot.label();
        update_state(opi, state_);
      }
      else if (dot.at_origin())
      {
        int dot_spin = state[dot.label()];
        ASSERT(state_[dot.site()] == dot_spin, "spin is not consistent");
      }
      label++;
    }
    ASSERT(is_same_state(state_, state, 0), "operators are not consistent while update worms");
#endif
    return;
  }

  bool detectWormCross(double tau, double tau_prime, double wt_tau, int dir)
  {
    if (dir == 1)
    {
      if ((tau_prime == 0 ? 1 : tau_prime) >= wt_tau && wt_tau > tau)
        return true;
      else
        return false;
    }
    else
    { // if dir == 0
      if (tau_prime <= wt_tau && wt_tau < (tau == 0 ? 1 : tau))
        return true;
      else
        return false;
    }
  }

  void reset_ops()
  {
    for (size_t i = 0; i < psop.size(); i++)
    {
      ops_main[psop[i]].set_state(pres[i]);
    }
  }

  bool is_same_state(int n, int m)
  {
    return n == m;
  }

  bool is_same_state(int n, state_t state, size_t lopt)
  {
    int m = state_funcs[lopt].state2num(state, state.size());
    return n == m;
  }

  bool is_same_state(state_t state_, state_t state, size_t lopt)
  {
    int m = state_funcs[lopt].state2num(state, state.size());
    int n = state_funcs[lopt].state2num(state_, state.size());
    return n == m;
  }
  static void printStateAtTime(const state_t &state, double time)
  {
#ifndef NDEBUG
    std::cout << "current spin at time :" << time << " is : ";
    for (auto spin : state)
    {
      std::cout << spin << " ";
    }
    std::cout << std::endl;
#else
    return;
#endif
  }
};

extern template class Worm<bcl::heatbath>;
extern template class Worm<bcl::st2010>;
extern template class Worm<bcl::st2013>;
