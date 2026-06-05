/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(oxdna/awsem/excv/dh,PairOxdnaAwsemExcvDh);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA_AWSEM_EXCV_DH_H
#define LMP_PAIR_OXDNA_AWSEM_EXCV_DH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdnaAwsemExcvDh : public Pair {
 public:
  PairOxdnaAwsemExcvDh(class LAMMPS *);
  ~PairOxdnaAwsemExcvDh() override;
  void init_style() override;
  void coeff(int, char **) override;
  void settings(int, char **) override;
  double init_one(int, int) override;
  virtual void compute_interaction_sites(double *, double *, double *, double *, double *);
  void compute(int, int) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void *extract(const char *, int &) override;

 protected:
  double cut_coul_global;
  double **nxyz_xtrct;    // per-atom arrays for local unit vectors
  double **epsilon, **sigma, **cut_lj, **cut_ljsq; // LJ parameters
  double lambda, kappa_one, **kappa, *qeff_dh, **cut_coul, **cut_coulsq; // DH parameters
  double **lj1, **lj2, **lj3, **lj4, **offset;

  virtual void allocate();

  // potential file reading
  virtual void read_file(char *);
  virtual void setup_params();
  int N_VALUES = 80;
  int N_AA_TYPES = 20; // A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V
  int N_DNA_TYPES = 4; // A, C, T, G
  int N_TYPES = 24;    // total number of types

  class FixBackbone *fix_bb; // ptr to awsemmd backbone fix
  int se_map[26] = {0, 0, 4, 3, 6, 13, 7, 8, 9, 0, 11, 10, 12, 2, 0, 14, 5, 1, 15, 16, 0, 19, 17, 0, 18, 0};
  char one_letter_code[20] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};

  class FixOxdnaLRF *fix_lrf;    // ptr to oxdna/lrf fix
};

} // namespace LAMMPS_NS

#endif
#endif
