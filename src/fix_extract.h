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

#ifdef FIX_CLASS
// clang-format off
FixStyle(extract,FixExtract);
// clang-format on
#else

#ifndef LMP_FIX_EXTRACT_H
#define LMP_FIX_EXTRACT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixExtract : public Fix {
 public:
  FixExtract(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void min_pre_force(int) override;
  double compute_scalar() override;

 private:
  bigint next_extract;
  char *var_extract, *fixname, *fieldname;
  int ivar_extract;
  class Fix *fstyle;
  double *fvoid;

  void extract_scalar();
};

}    // namespace LAMMPS_NS

#endif
#endif
