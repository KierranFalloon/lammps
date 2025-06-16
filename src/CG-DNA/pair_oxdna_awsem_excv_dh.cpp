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
/* ----------------------------------------------------------------------
  Contributing authors: Kierran Falloon (University of Strathclyde, Glasgow)
                        Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "pair_oxdna_awsem_excv_dh.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "constants_oxdna.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOxdnaAwsemExcvDh::PairOxdnaAwsemExcvDh(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  reinitflag = 1;
  offset_flag = 1;
  trim_flag = 0;
  file = nullptr;
}

/* ---------------------------------------------------------------------- */

PairOxdnaAwsemExcvDh::~PairOxdnaAwsemExcvDh()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOxdnaAwsemExcvDh::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");

  memory->create(epsilon, np1, np1, "pair:epsilon");
  memory->create(sigma, np1, np1, "pair:sigma");
  memory->create(cut_lj, np1, np1, "pair:cut_lj");
  memory->create(cut_ljsq, np1, np1, "pair:cut_ljsq");
  memory->create(lj1, np1, np1, "pair:lj1");
  memory->create(lj2, np1, np1, "pair:lj2");
  memory->create(lj3, np1, np1, "pair:lj3");
  memory->create(lj4, np1, np1, "pair:lj4");
  memory->create(offset, np1, np1, "pair:offset"); // E_cut
  memory->create(cut_coul, np1, np1, "pair:cut_coul");
  memory->create(cut_coulsq, np1, np1, "pair:cut_coulsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdnaAwsemExcvDh::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_coul_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_coul[i][j] = cut_coul_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

// pair_coeff ilo*ihi jlo*jhi oxdna2/awsem/excv/dh temperature salt_conc param_file
void PairOxdnaAwsemExcvDh::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR, "Incorrect args for pair_coeff of pair style oxdna/awsem/excv/dh");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  // Simulation specific parameters
  double T = utils::numeric(FLERR, arg[2], false, lmp);
  double rhos = utils::numeric(FLERR, arg[3], false, lmp);

  // Coulombic interaction
  // Debye-Huckel Parameters (eps_r = 80.0)
  lambda = ConstantsOxdna::get_lambda_dh_one_prefactor() * sqrt(T/0.1/rhos);
  kappa = 1.0/lambda; // inverse Debye length

  // populate 'file' with data from the potential file
  file = new File();
  read_file(arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = file->lj_epsilon[count];
      sigma[i][j] = file->lj_sigma[count];
      cut_lj[i][j] = sigma[i][j] * 1.144714243; // 1.1447 = (3/2)^(1/3)
      cut_coul[i][j] = cut_coul_global;
      setflag[i][j] = 1;
      count++;
    }
  }

  memory->destroy(file->lj_epsilon);
  memory->destroy(file->lj_sigma);
  delete file;

  if (count == 0) error->all(FLERR, "Incorrect args for pair_coeff of pair style oxdna/awsem/excv/dh");
}

void PairOxdnaAwsemExcvDh::read_file(char *filename)
{
  memory->create(file->lj_epsilon, N_VALUES, "pair:lj_epsilon");
  memory->create(file->lj_sigma, N_VALUES, "pair:lj_sigma");

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, filename, "oxdna/awsem/excv/dh");
    try {
      reader.skip_line();
      reader.skip_line();

      // read in the lj epsilon, sigma and dh qeff values
      reader.next_dvector(file->lj_epsilon, N_VALUES);
      reader.next_dvector(file->lj_sigma, N_VALUES);

    } catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }

  MPI_Bcast(file->lj_epsilon, N_VALUES, MPI_DOUBLE, 0, world);
  MPI_Bcast(file->lj_sigma, N_VALUES, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdnaAwsemExcvDh::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  double cut = MAX(cut_lj[i][j], cut_coul[i][j]);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];

  lj1[i][j] = 36.0 * epsilon[i][j] * pow(sigma[i][j], 9.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 9.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 9.0) - pow(ratio, 6.0));
  } else
    offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  cut_coulsq[j][i] = cut_coulsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
    compute vector COM-excluded volume interaction sites in oxDNA
------------------------------------------------------------------------- */
void PairOxdnaAwsemExcvDh::compute_interaction_sites(double e1[3], double /*e2*/[3],
    double /*e3*/[3], double rs[3], double rb[3])
{
  double d_cs = ConstantsOxdna::get_d_cs();
  double d_cb = ConstantsOxdna::get_d_cb();

  rs[0] = d_cs*e1[0];
  rs[1] = d_cs*e1[1];
  rs[2] = d_cs*e1[2];

  rb[0] = d_cb*e1[0];
  rb[1] = d_cb*e1[1];
  rb[2] = d_cb*e1[2];

}

void PairOxdnaAwsemExcvDh::init_style()
{
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ----------------------------------------------------------------------
   compute function for AWSEM-oxDNA excluded volume & electrostatic interactions
   s=DNA sugar-phosphate backbone site, b=DNA base site
------------------------------------------------------------------------- */

void PairOxdnaAwsemExcvDh::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,iellipsoid,jellipsoid;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r2inv_s,r2inv_b,r6inv_s,r6inv_b,r3inv_s,r3inv_b,forcecoul,forcelj,factor_coul,factor_lj;
  double r,rinv,screening;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  double **torque = atom->torque;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double rtmp_s[3],rtmp_b[3];
  double delr_s[3],rsq_s;
  double delr_b[3],rsq_b;
  double delf[3],delti[3],deltj[3]; // force, torque increment;

  // vectors COM-backbone site, COM-base site in lab frame
  double ri_cs[3],ri_cb[3];
  double rj_cs[3],rj_cb[3];
  // Cartesian unit vectors in lab frame
  double ix[3],iy[3],iz[3];
  double jx[3],jy[3],jz[3];

  auto avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    iellipsoid = ellipsoid[i]; // < 0 if not an ellipsoid
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (iellipsoid >= 0) { // ellipsoid
      double *qn, nx_temp[3], ny_temp[3], nz_temp[3]; // quaternion and Cartesian unit vectors in lab frame

      qn = bonus[iellipsoid].quat;
      MathExtra::q_to_exyz(qn, nx_temp, ny_temp, nz_temp);

      ix[0] = nx_temp[0];
      ix[1] = nx_temp[1];
      ix[2] = nx_temp[2];
      iy[0] = ny_temp[0];
      iy[1] = ny_temp[1];
      iy[2] = ny_temp[2];
      iz[0] = nz_temp[0];
      iz[1] = nz_temp[1];
      iz[2] = nz_temp[2];

      compute_interaction_sites(ix, iy, iz, ri_cs, ri_cb);

      // vector COM-backbone site
      rtmp_s[0] = x[i][0] + ri_cs[0];
      rtmp_s[1] = x[i][1] + ri_cs[1];
      rtmp_s[2] = x[i][2] + ri_cs[2];

      // vector COM-base site
      rtmp_b[0] = x[i][0] + ri_cb[0];
      rtmp_b[1] = x[i][1] + ri_cb[1];
      rtmp_b[2] = x[i][2] + ri_cb[2];
    } else { // non ellipsoid
      // protein atom COM
      rtmp_s[0] = x[i][0];
      rtmp_s[1] = x[i][1];
      rtmp_s[2] = x[i][2];

      rtmp_b[0] = x[i][0];
      rtmp_b[1] = x[i][1];
      rtmp_b[2] = x[i][2];
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      jtype = type[j];
      jellipsoid = ellipsoid[j]; // < 0 if not an ellipsoid

      // sanity check - should always be ellipsoid <-> non-ellipsoid interaction
      if ((iellipsoid >= 0 && jellipsoid >= 0) || (iellipsoid < 0 && jellipsoid < 0)) {
        error->all(FLERR, "Detected a pair of ellipsoids or a pair of non-ellipsoids interacting via the oxdna/awsem/excv/dh potential");
      }

      if (jellipsoid >= 0) { // atom j is an ellipsoid, so atom i is not
        double *qn, nx_temp[3], ny_temp[3], nz_temp[3]; // quaternion and Cartesian unit vectors in lab frame

        qn = bonus[jellipsoid].quat;
        MathExtra::q_to_exyz(qn, nx_temp, ny_temp, nz_temp);

        jx[0] = nx_temp[0];
        jx[1] = nx_temp[1];
        jx[2] = nx_temp[2];
        jy[0] = ny_temp[0];
        jy[1] = ny_temp[1];
        jy[2] = ny_temp[2];
        jz[0] = nz_temp[0];
        jz[1] = nz_temp[1];
        jz[2] = nz_temp[2];

        compute_interaction_sites(jx, jy, jz, rj_cs, rj_cb);
        
        // as atom i is not an ellipsoid, rtmp_s and rtmp_b are the protein atom COM coordinates
        // protein atom COM - DNA backbone site
        delr_s[0] = rtmp_s[0] - (x[j][0] + rj_cs[0]);
        delr_s[1] = rtmp_s[1] - (x[j][1] + rj_cs[1]);
        delr_s[2] = rtmp_s[2] - (x[j][2] + rj_cs[2]);

        // protein atom COM - DNA base site
        delr_b[0] = rtmp_b[0] - (x[j][0] + rj_cb[0]);
        delr_b[1] = rtmp_b[1] - (x[j][1] + rj_cb[1]);
        delr_b[2] = rtmp_b[2] - (x[j][2] + rj_cb[2]);
      } else { // j is not an ellipsoid, so i is

        // as atom i is an ellipsoid, rtmp_s and rtmp_b are the DNA backbone and base site coordinates
        // DNA backbone site - protein atom COM
        delr_s[0] = rtmp_s[0] - x[j][0];
        delr_s[1] = rtmp_s[1] - x[j][1];
        delr_s[2] = rtmp_s[2] - x[j][2];

        // DNA base site - protein atom COM
        delr_b[0] = rtmp_b[0] - x[j][0];
        delr_b[1] = rtmp_b[1] - x[j][1];
        delr_b[2] = rtmp_b[2] - x[j][2];
      }
      
      // squared protein COM - DNA backbone and base site distances
      rsq_s = delr_s[0]*delr_s[0] + delr_s[1]*delr_s[1] + delr_s[2]*delr_s[2];
      rsq_b = delr_b[0]*delr_b[0] + delr_b[1]*delr_b[1] + delr_b[2]*delr_b[2];

      // backbone-COM interactions
      if (rsq_s < cutsq[itype][jtype]) {
        r2inv_s = 1.0/rsq_s;
        
        // Debye-Huckel
        if (rsq_s < cut_coulsq[itype][jtype]) {
          r = sqrt(rsq_s);
          rinv = 1.0/r;
          screening = exp(-kappa*r);
          // Note: ConstantsOxdna::qeff_dh_pf_one_prefactor === force->qqrd2e for eps_r = 80.0
          forcecoul = ConstantsOxdna::get_qeff_dh_pf_one_prefactor() * qtmp*q[j] * screening * (kappa + rinv);
        } else forcecoul = 0.0;

        // Lennard-Jones
        if (rsq_s < cut_ljsq[itype][jtype]) {
          r6inv_s = r2inv_s*r2inv_s*r2inv_s;
          r3inv_s = sqrt(r6inv_s);
          forcelj = r6inv_s * (lj1[itype][jtype]*r3inv_s - lj2[itype][jtype]);
        } else forcelj = 0.0;

        // Both Debye-Huckel and Lennard-Jones excluded volume between COM and backbone site
        fpair = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv_s;

        delf[0] = delr_s[0]*fpair;
        delf[1] = delr_s[1]*fpair;
        delf[2] = delr_s[2]*fpair;

        f[i][0] += delf[0];
        f[i][1] += delf[1];
        f[i][2] += delf[2];

        if (iellipsoid >= 0) {
          MathExtra::cross3(ri_cs,delf,delti);
          torque[i][0] += delti[0];
          torque[i][1] += delti[1];
          torque[i][2] += delti[2];
        }

        if (newton_pair || j < nlocal) {
          f[j][0] -= delf[0];
          f[j][1] -= delf[1];
          f[j][2] -= delf[2];

          if (jellipsoid >= 0) {
            MathExtra::cross3(rj_cs,delf,deltj);
            torque[j][0] -= deltj[0];
            torque[j][1] -= deltj[1];
            torque[j][2] -= deltj[2];
          }
        }

        if (eflag) {
          if (rsq_s < cut_coulsq[itype][jtype])
            ecoul = factor_coul * qqrd2e * qtmp*q[j] * rinv * screening;
          else ecoul = 0.0;
          if (rsq_s < cut_ljsq[itype][jtype]) {
            evdwl = r6inv_s*(lj3[itype][jtype]*r3inv_s-lj4[itype][jtype]) -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,ecoul,
            delf[0],delf[1],delf[2],delx,dely,delz);

      }

      // base-COM interactions
      if (rsq_b < cutsq[itype][jtype]) {
        r2inv_b = 1.0/rsq_b;

        // Lennard-Jones
        if (rsq_b < cut_ljsq[itype][jtype]) {
          r6inv_b = r2inv_b*r2inv_b*r2inv_b;
          r3inv_b = sqrt(r6inv_b);
          forcelj = r6inv_b * (lj1[itype][jtype]*r3inv_b - lj2[itype][jtype]);
        } else forcelj = 0.0;
        
        // Only Lennard-Jones excluded volume between COM and base site
        fpair = (factor_lj*forcelj) * r2inv_b;

        delf[0] = delr_b[0]*fpair;
        delf[1] = delr_b[1]*fpair;
        delf[2] = delr_b[2]*fpair;

        f[i][0] += delf[0];
        f[i][1] += delf[1];
        f[i][2] += delf[2];

        if (iellipsoid >= 0) {
          MathExtra::cross3(ri_cb,delf,delti);
          torque[i][0] += delti[0];
          torque[i][1] += delti[1];
          torque[i][2] += delti[2];
        }

        if (newton_pair || j < nlocal) {
          f[j][0] -= delf[0];
          f[j][1] -= delf[1];
          f[j][2] -= delf[2];

          if (jellipsoid >= 0) {
            MathExtra::cross3(rj_cb,delf,deltj);
            torque[j][0] -= deltj[0];
            torque[j][1] -= deltj[1];
            torque[j][2] -= deltj[2];
          }
        }

        if (eflag) {
          if (rsq_b < cut_ljsq[itype][jtype]) {
            evdwl = r6inv_b*(lj3[itype][jtype]*r3inv_b-lj4[itype][jtype]) -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,
            delf[0],delf[1],delf[2],delx,dely,delz);

      }

    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairOxdnaAwsemExcvDh::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) {
    fprintf(fp, "  %d %d %g %g %g %g %g\n",
      i, epsilon[i][i], sigma[i][i], cut_lj[i][i], cut_coul[i][i]);
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdnaAwsemExcvDh::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g\n",
        i, j, epsilon[i][j], sigma[i][j], cut_lj[i][j], cut_coul[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairOxdnaAwsemExcvDh::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return &epsilon[0][0];
  if (strcmp(str, "sigma") == 0) return &sigma[0][0];
  if (strcmp(str, "cut_lj") == 0) return &cut_lj[0][0];
  if (strcmp(str, "cut_coul") == 0) return &cut_coul[0][0];
  return nullptr;
}
