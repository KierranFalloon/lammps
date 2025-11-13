// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_extract.h"

#include "error.h"
#include "input.h"
#include "fix.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   this is a fix that allows for the extraction of a scalar value from
   another fix at regular intervals during a run through Fix::extract()
------------------------------------------------------------------------- */

FixExtract::FixExtract(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  scalar_flag = 1;

  if (narg < 4) utils::missing_cmd_args(FLERR, "fix extract", error);

  if (utils::strmatch(arg[3], "^v_")) {
    var_extract = utils::strdup(arg[3] + 2);
    nevery = 1;
  } else {
    nevery = utils::inumeric(FLERR,arg[3],false,lmp);
    if (nevery <= 0)
      error->all(FLERR,"Illegal fix extract nevery value: {}", nevery);
  }

  fixname = utils::strdup(arg[4]);
  fstyle = modify->get_fix_by_id(fixname);
  if (!fstyle || !fstyle->scalar_flag)
    error->all(FLERR,"Fix ID {} for fix extract does not exist or compute a scalar", fixname);
  extscalar = fstyle->extscalar;

  fieldname = utils::strdup(arg[5]);
}

FixExtract::~FixExtract()
{
  delete[] var_extract;
  delete[] fixname;
  delete[] fieldname;
}

/* ---------------------------------------------------------------------- */

int FixExtract::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExtract::init()
{
  // decide next timestep for extracting a value
  if (var_extract) {
    ivar_extract = input->variable->find(var_extract);
    if (ivar_extract < 0)
      error->all(FLERR, "Variable {} for fix extract timestep does not exist", var_extract);
    if (!input->variable->equalstyle(ivar_extract))
      error->all(FLERR, "Variable {} for fix extract timestep is invalid style", var_extract);
    next_extract = static_cast<bigint>(input->variable->compute_equal(ivar_extract));
    if (next_extract <= update->ntimestep)
      error->all(FLERR, "Fix extract timestep variable {} returned a bad timestep: {}", var_extract,
                 next_extract);
  } else {
    if (update->ntimestep == 0 || update->ntimestep % nevery)
      next_extract = (update->ntimestep / nevery) * nevery + nevery;
    else
      next_extract = update->ntimestep;
  }

  // ensure extractable field exists and is initially available
  extract_scalar();
  
}

/* ---------------------------------------------------------------------- */

void FixExtract::setup_pre_force(int /*vflag*/)
{
  extract_scalar();
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixExtract::pre_force(int /*vflag*/)
{
  if (update->ntimestep != next_extract) return;

  if (var_extract) {
    next_extract = static_cast<bigint>(input->variable->compute_equal(ivar_extract));
    if (next_extract <= update->ntimestep)
      error->all(FLERR, "Fix extract timestep variable returned a bad timestep: {}", next_extract);
  } else {
    next_extract = (update->ntimestep / nevery) * nevery + nevery;
  }

  int dim;
  fvoid = (double *) fstyle->extract(fieldname, dim);

}

/* ---------------------------------------------------------------------- */

void FixExtract::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixExtract::compute_scalar()
{
  if (!fvoid) return 0.0;
  return fvoid[0];
}

void FixExtract::extract_scalar()
{
  int dim;
  fvoid = (double *) fstyle->extract(fieldname, dim);
  if (!fvoid)
    error->all(FLERR,"Cannot extract field {} from fix {}", fieldname, fixname);
}
