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
  double **nx_xtrct, **ny_xtrct, **nz_xtrct;    // per-atom arrays for local unit vectors
  double **epsilon, **sigma, **cut_lj, **cut_ljsq; // LJ parameters
  double lambda, kappa_one, **kappa, **qeff_dh, **cut_coul, **cut_coulsq; // DH parameters
  double **lj1, **lj2, **lj3, **lj4, **offset;

  virtual void allocate();

  // potential file reading
  virtual void read_file(char *);
  int N_VALUES = 20;
  struct File { // potential file data
   double *lj_epsilon;
   double *lj_sigma;
   double *qeff_dh;
  };
  File *file;

  class Fix *fix_lrf;    // ptr to oxdna/lrf fix
};

} // namespace LAMMPS_NS

#endif
#endif
