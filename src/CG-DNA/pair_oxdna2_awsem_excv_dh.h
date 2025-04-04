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
PairStyle(oxdna2/awsem/excv/dh,PairOxdna2AwsemExcvDh);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA2_AWSEM_EXCV_DH_H
#define LMP_PAIR_OXDNA2_AWSEM_EXCV_DH_H

#include "pair_oxdna_awsem_excv_dh.h"

namespace LAMMPS_NS {

class PairOxdna2AwsemExcvDh : public PairOxdnaAwsemExcvDh {
 public:
  PairOxdna2AwsemExcvDh(class LAMMPS *lmp) : PairOxdnaAwsemExcvDh(lmp) {}

  void compute_interaction_sites(double *, double *, double *, double *, double *);
};

} // namespace LAMMPS_NS

#endif
#endif
