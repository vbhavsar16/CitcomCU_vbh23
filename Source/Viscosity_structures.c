/*
 * CitcomCU is a Finite Element Code that solves for thermochemical
 * convection within a three dimensional domain appropriate for convection
 * within the Earth's mantle. Cartesian and regional-spherical geometries
 * are implemented. See the file README contained with this distribution
 * for further details.
 * 
 * Copyright (C) 1994-2005 California Institute of Technology
 * Copyright (C) 2000-2005 The University of Colorado
 *
 * Authors: Louis Moresi, Shijie Zhong, and Michael Gurnis
 *
 * For questions or comments regarding this software, you may contact
 *
 *     Luis Armendariz <luis@geodynamics.org>
 *     http://geodynamics.org
 *     Computational Infrastructure for Geodynamics (CIG)
 *     California Institute of Technology
 *     2750 East Washington Blvd, Suite 210
 *     Pasadena, CA 91007
 *
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 2 of the License, or any
 * later version.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program; if not, write to the Free Software 
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/* Functions relating to the determination of viscosity field either
   as a function of the run, as an initial condition or as specified from
   a previous file */

#include <mpi.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"

/* ==================================================================================================== */
void viscosity_parameters(struct All_variables *E)
{
  int i, l;
  float temp;
  
  /* default values .... */
  E->viscosity.update_allowed = 0;
  E->viscosity.SDEPV = E->viscosity.TDEPV = 0;
  E->viscosity.EXPX = 0;
  E->viscosity.TSCDEPV = 0; /* MIB, JAN07 */
  E->control.sdepv_print_convergence = 0;  /* MIB, JAN07 */
  E->control.output_sdepvnewt_stepzero = 0; /* MIB, JAN07 */
  
  for(i = 0; i < 40; i++)
    {
      E->viscosity.N0[i] = 1.0;
      E->viscosity.T[i] = 0.0;
      E->viscosity.Z[i] = 0.0;
      E->viscosity.E[i] = 0.0;
      E->viscosity.T0[i] = 0.0;
      E->viscosity.V[i] = 0.0; /* MIB, JAN07 */
      E->viscosity.COH[i] = 0.0; /* MIB, JAN07 */
      E->viscosity.d[i] = 0.0; /* MIB, JAN07 */
      E->viscosity.sdepv_expt[i] = 0.0; /* MIB, JAN07 */
      E->viscosity.coh_expt[i] = 0.0; /* MIB, JAN07 */
      E->viscosity.grainsize_expt[i] = 0.0; /* MIB, JAN07 */
    }
  
  /* read in information */
  input_int("rheol", &(E->viscosity.RHEOL), "essential");
  input_int("num_mat", &(E->viscosity.num_mat), "1");
  input_int("weak_regions", &(E->viscosity.weak_region),"0"); 
 
  input_float_vector("viscT", E->viscosity.num_mat, (E->viscosity.T));	
  input_float_vector("viscT1", E->viscosity.num_mat, (E->viscosity.T)); /* redundant */
  input_float_vector("viscZ", E->viscosity.num_mat, (E->viscosity.Z));
  input_float_vector("viscE", E->viscosity.num_mat, (E->viscosity.E));
  input_float_vector("viscT0", E->viscosity.num_mat, (E->viscosity.T0)); /* not used */
  input_float_vector("visc0", E->viscosity.num_mat, (E->viscosity.N0));	/* redundant */
  input_float_vector("viscN0", E->viscosity.num_mat, (E->viscosity.N0));
  input_float_vector("viscV", E->viscosity.num_mat, (E->viscosity.V)); /* MIB, JAN07 */
  input_float_vector("viscCOH", E->viscosity.num_mat, (E->viscosity.COH)); /* MIB, JAN07 */
  input_float_vector("viscd", E->viscosity.num_mat, (E->viscosity.d)); /* MIB, JAN07 */

  input_boolean("TDEPV", &(E->viscosity.TDEPV), "on");
  input_boolean("SDEPV", &(E->viscosity.SDEPV), "off");
  input_boolean("TSCDEPV", &(E->viscosity.TSCDEPV), "off"); /* MIB, JAN07 */
  input_int("num_creepmech",&(E->viscosity.num_crpmech),"2"); /*MIB, JAN07 */
  input_int("sdepv_see_convergence",&(E->control.sdepv_print_convergence),"1"); /* MIB, JAN07 */
  input_int("output_sdepvnewt_stepzero", &(E->control.output_sdepvnewt_stepzero),"0"); /* MIB, JAN07 */

  input_boolean("YIELD", &(E->viscosity.YIELD), "off"); /* MIB, JAN07 */
  input_float_vector("coh_expt", E->viscosity.num_mat, (E->viscosity.coh_expt)); /* MIB, JAN07 */
  input_float_vector("grainsize_expt", E->viscosity.num_mat, (E->viscosity.grainsize_expt)); /* MIB, JAN07 */

  input_float("sdepv_misfit", &(E->viscosity.sdepv_misfit), "0.001");
  input_float_vector("sdepv_expt", E->viscosity.num_mat, (E->viscosity.sdepv_expt));
  input_float_vector("sdepv_trns", E->viscosity.num_mat, (E->viscosity.sdepv_trns));
 
  input_float("yield_surf", &(E->viscosity.yield_surf), "0.1e6"); /* MIB, JAN 07 */
  input_float("yield_mantle",  &(E->viscosity.yield_mantle), "1.0e9"); /* MIB, JAN 07 */
  input_float("yield_grad",  &(E->viscosity.yield_grad), "15e3"); /* MIB, JAN 07 */
  input_float("yield_temp",  &(E->viscosity.yield_temp), "1400"); /* MIB, JAN 07 */
  
  input_boolean("TDEPV_AVE", &(E->viscosity.TDEPV_AVE), "off");
  input_boolean("VFREEZE", &(E->viscosity.FREEZE), "off");
  input_boolean("VMAX", &(E->viscosity.MAX), "off");
  input_boolean("VMIN", &(E->viscosity.MIN), "off");
  input_boolean("VISC_UPDATE", &(E->viscosity.update_allowed), "on");
  
  input_float("freeze_thresh", &(E->viscosity.freeze_thresh), "0.0");
  input_float("freeze_value", &(E->viscosity.freeze_value), "1.0");
  input_float("visc_max", &(E->viscosity.max_value), "nodefault");
  input_float("visc_min", &(E->viscosity.min_value), "nodefault");
  
  input_boolean("VISC_GUESS", &(E->viscosity.guess), "off");
  input_string("visc_old_file", E->viscosity.old_file, " ");
  
  if (!E->viscosity.TSCDEPV) { /* MIB, JAN07; this is done inside composite_viscosity */

    if(E->viscosity.RHEOL == 1)
      {
	for(l = 1; l <= E->viscosity.num_mat; l++)
	  {
	    E->viscosity.E[l - 1] = E->viscosity.E[l - 1] / (E->data.gas_const * E->data.ref_temperature);
	    E->viscosity.T[l - 1] = E->viscosity.T[l - 1] / E->data.ref_temperature;
	    temp = exp(E->viscosity.E[l - 1] / (1.0 + E->viscosity.T[l - 1]));
	    E->viscosity.N0[l - 1] = E->viscosity.N0[l - 1] / temp;
	  }
      }
    else if(E->viscosity.RHEOL == 2)
      {
	if(E->parallel.me == 0)
	  fprintf(stderr, "this option for rheology is not supported\n");
	parallel_process_termination();
      }
    else if(E->viscosity.RHEOL == 3)
      {
	for(l = 1; l <= E->viscosity.num_mat; l++)
	  {
	    E->viscosity.E[l - 1] = E->viscosity.E[l - 1] / (E->data.gas_const * E->data.ref_temperature);
	    E->viscosity.Z[l - 1] = E->viscosity.Z[l - 1] * E->data.density * E->data.grav_acc * 
	      E->monitor.length_scale / (E->data.gas_const * E->data.ref_temperature);
	  }
      }
  }

  return;
}

/* ==================================================================================================== */
void get_viscosity_option(struct All_variables *E)
{
  int i;

	/* general, essential default */

	input_string("Viscosity", E->viscosity.STRUCTURE, NULL);	/* Which form of viscosity */

	input_boolean("VISC_EQUIVDD", &(E->viscosity.EQUIVDD), "off");	/* Whether to average it */
	input_int("equivdd_opt", &(E->viscosity.equivddopt), "1");
	input_int("equivdd_x", &(E->viscosity.proflocx), "1");
	input_int("equivdd_y", &(E->viscosity.proflocy), "1");

	input_int("update_every_steps", &(E->control.KERNEL), "1");

	input_boolean("VISC_SMOOTH", &(E->viscosity.SMOOTH), "off");	
	input_int("visc_smooth_loops", &(E->viscosity.smooth_loops),"0");	/*MAJ 8/8/07 number loops in apply_viscosity_smoother*/
	if(E->parallel.me == 0)
		fprintf(E->fp,"VISC_SMOOTH = %d and visc_smooth_loops = %d\n", E->viscosity.SMOOTH, E->viscosity.smooth_loops);		/*MAJ 8/15/07 Debugging */
	input_int("visc_smooth_cycles", &(E->viscosity.smooth_cycles), "0");

	if (E->viscosity.weak_region){
	get_weak_regions(E);  /* read from file */
	if(E->parallel.me == 0)
		fprintf(E->fp,"get_weak_regions is called. Flag E->viscosity.weak_region = %d\n",E->viscosity.weak_region);
	}

	for (i = 1; i<=E->lmesh.nel; i++)  /* Changed from nno to nel to match declaration. 9/12/12 MAJ */
	  E->visctyp[i] = 0.0;

	if(strcmp(E->viscosity.STRUCTURE, "system") == 0)	/* Interpret */
	{
		if(E->parallel.me == 0)
			fprintf(E->fp, "Viscosity derived from system state\n");
		E->viscosity.FROM_SYSTEM = 1;
		viscosity_for_system(E);
	}

	return;

}

/* ==================================================================================================== */

void viscosity_for_system(struct All_variables *E)
{
  if(!E->viscosity.update_allowed)
    {
      get_system_viscosity(E, 1, E->EVI[E->mesh.levmax], E->VI[E->mesh.levmax]);
    }
  //else /* MIB, JAN07: added else, before called get_system_viscosity twice is update_allowed is false */
  //  {
    //  get_system_viscosity(E, 1, E->EVI[E->mesh.levmax], E->VI[E->mesh.levmax]);
    //}
  return;
}

/* ==================================================================================================== */
void get_system_viscosity(struct All_variables *E, int propogate, float *evisc, float *visc)
{
  int i, j;
  //float *visc_old, *evisc_old;
  
  const int vpts = vpoints[E->mesh.nsd];
  
  if (E->viscosity.TSCDEPV) { /* MIB, JAN07, composite temp, stress, composition dep. viscosity */
    composite_viscosity(E, visc, evisc);
  }
  else {
    if(E->viscosity.TDEPV)
      visc_from_T(E, visc, evisc, propogate);
    else
      visc_from_mat(E, visc, evisc);
  
    if(E->viscosity.SDEPV)
      visc_from_S(E, visc, evisc, propogate);
  }

  if(E->viscosity.weak_region)
    weak_regions(E, evisc);
  
  if(E->viscosity.SMOOTH)
    apply_viscosity_smoother(E, visc, evisc);
  
  if(E->viscosity.MAX)
    {
      for(i = 1; i <= E->lmesh.nel; i++)
	for(j = 1; j <= vpts; j++)
	  {
	    if(evisc[(i - 1) * vpts + j] > E->viscosity.max_value)
	      evisc[(i - 1) * vpts + j] = E->viscosity.max_value;
	  }
    }
  
  if(E->viscosity.MIN)
    {
      for(i = 1; i <= E->lmesh.nel; i++)
	for(j = 1; j <= vpts; j++)
	  if(evisc[(i - 1) * vpts + j] < E->viscosity.min_value)
	    evisc[(i - 1) * vpts + j] = E->viscosity.min_value;
    }
  
  /*   v_to_nodes(E,evisc,visc,E->mesh.levmax);  */
  
  return;
}

/* ==================================================================================================== */
void apply_viscosity_smoother(struct All_variables *E, float *visc, float *evisc)
{
  double *ViscCentre;
  int i;
  
  ViscCentre = (double *)malloc((E->lmesh.nno + 10) * sizeof(double));
 
  /*fprintf(stderr,"You are in apply_viscosity_smoother DESPITE THAT THIS LOOP IS NOT BEING CALLED. F.O. TESTING- MAJ- 8/15/07/n");
  fprintf(stderr,"and flag E->viscosity.SMOOTH = %d, and E->viscosity.smooth_loops is %d\n", E->viscosity.SMOOTH, E->viscosity.smooth_loops);
  */

  /*for(i = 1; i <= E->viscosity.smooth_cycles; i++) Pulled to test 7/27/07 MAJ*/
  for(i = 1; i <= E->viscosity.smooth_loops; i++) /*3 in 7/27/07 MAJ. 1 in 7/30/2007. 3 in 7/31/07. 1 in 8/1/07. viscosity.smooth_loops 8/8/07 */
    {
      p_to_centres(E, visc, ViscCentre, E->mesh.levmax);
      p_to_nodes(E, ViscCentre, visc, E->mesh.levmax);
    }
  
  free((void *)ViscCentre);
  
  return;
}

/* ==================================================================================================== */
void visc_from_mat(struct All_variables *E, float *Eta, float *EEta)
{
  //int i, j, k, l, z, jj, kk;
  int i, jj;
  
  for(i = 1; i <= E->lmesh.nel; i++)
    for(jj = 1; jj <= vpoints[E->mesh.nsd]; jj++)
      {
	EEta[(i - 1) * vpoints[E->mesh.nsd] + jj] = E->viscosity.N0[E->mat[i] - 1];
      }
  
  return;
}
/* ==================================================================================================== */
void visc_from_T(struct All_variables *E, float *Eta, float *EEta, int propogate)
{
  //int i, j, k, l, z, e, jj, kk, imark;
  int i, l, e, jj, kk;
  //float c1, c2, c3, zero, e_6, one, eta0, Tave, depth, temp, tempa, TT[9];
  float zero, one, temp, tempa, TT[9];
  //double ztop, zbotm, zz, visc1, area1, temp1, temp2, temp3, temp4;
  double ztop, zbotm, zz, temp1, temp2, temp3, temp4;
  float *Xtmp[4];
  static int visits = 0;
  static float *Tadi;
  static double slope = 0;
  const int vpts = vpoints[E->mesh.nsd];
  const int ends = enodes[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  const int noz = E->lmesh.noz;
  
  one = 1.0;
  zero = 0.0;
  temp3 = temp4 = 0.0;
  
  if(visits == 0)
    {
      Tadi = (float *)malloc((noz + 1) * sizeof(float));
      
      if(E->control.Rsphere)
	slope = (1.0 - E->data.visc_factor) / (E->sphere.ro - E->sphere.ri);
      else if(E->control.CART3D)
	slope = (1.0 - E->data.visc_factor);
      
      for(i = noz; i >= 1; i--)
	E->Have.Tadi[i] = 0.0;
      
      E->data.T_adi0 = 0;
      E->data.T_adi1 = 0;
      
      return_horiz_ave(E, E->T, E->Have.T);

      
      fprintf(E->fp, "\tRheological option : %d\n", E->viscosity.RHEOL);
      
      for(l = 1; l <= E->viscosity.num_mat; l++)
	{
	  fprintf(E->fp, "\tlayer %d/%d: E=%g T1=%g N0=%g Z0=%g\n", l, E->viscosity.num_mat, 
		  E->viscosity.E[l - 1], E->viscosity.T[l - 1], E->viscosity.N0[l - 1], 
		  E->viscosity.Z[l - 1]);
	}
      fprintf(E->fp, "\tslope : %g %g\n", slope, E->data.visc_factor);
      fflush(E->fp);
      
    }
  
  if(E->control.Rsphere)
    {
      for(i = 1; i <= E->mesh.nsd; i++)
	Xtmp[i] = E->SX[i];
      ztop = E->sphere.ro;
      zbotm = E->sphere.ri;
    }
  else if(E->control.CART3D)
    {
      for(i = 1; i <= E->mesh.nsd; i++)
	Xtmp[i] = E->X[i];
      ztop = 1.0;
      zbotm = 0.0;
    }
  
  if(visits == 0 || E->monitor.solution_cycles % E->control.KERNEL == 0)
    {
      
      if(E->viscosity.RHEOL == 0)
		{
		  
		  temp2 = one;
		  
		  if(E->control.adi_heating)
		    {
		      temp = 0;
		      E->data.T_adi0 = 0;
		      
		      if(E->parallel.me_loc[3] == E->parallel.nprocz - 1)
			{
			  for(i = noz; i > 1; i--)
			    if(Xtmp[3][i] < (2 * E->viscosity.zlith - ztop))
			      {
				temp = E->Have.T[i];
				break;
			      }
			  E->data.T_adi0 = temp;
			}
		      
		      E->Have.Tadi[noz] = temp;
		      for(i = noz; i > 1; i--)
			{
			  if(Xtmp[3][i] < (2 * E->viscosity.zlith - ztop))
			    temp = temp + E->data.disptn_number * (E->expansivity[i] + E->expansivity[i - 1]) * 
			      (E->Have.T[i] + E->Have.T[i - 1] + 2 * E->data.surf_temp) / 
			      4.0 * (Xtmp[3][i] - Xtmp[3][i - 1]);
			  E->Have.Tadi[i - 1] = temp;
			}
		      
		      propogator_down_process(E, E->Have.Tadi);
		      
		      temp2 = one - E->data.T_adi1 + E->data.T_adi0;
		      
		    }					// end for adi_heating
		  
		  for(i = 1; i <= nel; i++)
		    {
		      l = E->mat[i];
		      e = (i - 1) % E->lmesh.elz + 1;
		      
		      tempa = E->viscosity.N0[l - 1];
		      temp1 = (E->Have.Tadi[e] + E->Have.Tadi[e + 1]) * 0.5 - E->data.T_adi0;
		      
		      for(kk = 1; kk <= ends; kk++)
			TT[kk] = E->T[E->ien[i].node[kk]];
		      
		      for(jj = 1; jj <= vpts; jj++)
			{
			  temp = 1.0e-32;
			  zz = 0.0;
			  for(kk = 1; kk <= ends; kk++)
			    {
			      temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
			      zz += Xtmp[3][E->ien[i].node[kk]] * E->N.vpt[GNVINDEX(kk, jj)];
			    }
			  EEta[(i - 1) * vpts + jj] = tempa * exp(E->viscosity.E[l - 1] * (temp2 - (temp - temp1)) / 
								  temp2) * (E->data.visc_factor + slope * (zz - zbotm));
			}
		    }
		}
      else if(E->viscosity.RHEOL == 1)
	{
	  for(i = 1; i <= nel; i++)
	    {
	      l = E->mat[i];
	      tempa = E->viscosity.N0[l - 1];
	      
	      for(kk = 1; kk <= ends; kk++)
		TT[kk] = E->T[E->ien[i].node[kk]];
	      
	      for(jj = 1; jj <= vpts; jj++)
		{
		  temp = 1.0e-32;
		  for(kk = 1; kk <= ends; kk++)
		    {
		      temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
		    }
		  EEta[(i - 1) * vpts + jj] = tempa * exp((E->viscosity.E[l - 1]) / (temp + E->viscosity.T[l - 1]));
		}
	    }
	}
      else if(E->viscosity.RHEOL == 2)
	{
	}
      else if(E->viscosity.RHEOL == 3)
	{
	  for(i = 1; i <= nel; i++)
	    {
	      l = E->mat[i];
	      tempa = E->viscosity.N0[l - 1];
	      
	      for(kk = 1; kk <= ends; kk++)
		TT[kk] = E->T[E->ien[i].node[kk]];
	      
	      for(jj = 1; jj <= vpts; jj++)
		{
		  temp = 1.0e-32;
		  zz = 0;
		  for(kk = 1; kk <= ends; kk++)
		    {
		      temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
		      zz += Xtmp[3][E->ien[i].node[kk]] * E->N.vpt[GNVINDEX(kk, jj)];
		    }
		  EEta[(i - 1) * vpts + jj] = tempa * exp((E->viscosity.E[l - 1] + (1 - zz) * E->viscosity.Z[l - 1]) / 
							  temp - (E->viscosity.E[l - 1] + E->viscosity.Z[l - 1]));
		  /*               EEta[(i-1)*vpts + jj] = tempa*exp(E->viscosity.E[l-1]*(E->viscosity.T[l-1]-temp));  */
		}
	    }
	}
      
      visits++;
      
    }
  
  
  /*
    fprintf(E->fp,"aaa\n"); 
    for(i=1;i<=nel;i++)   
    fprintf(E->fp,"%d %d %g\n",i,E->mat[i],EEta[(i-1)*vpts+1]); 
  */
  
  return;
}

/* ==================================================================================================== */
void visc_from_S(struct All_variables *E, float *Eta, float *EEta, int propogate)
{
  static int visits = 0;
  //float one, two, scale, stress_magnitude, depth, exponent1;
  float one, two, scale, exponent1;
  float *eedot;
  
  //int e, l, z, jj, kk;
  int e, jj;
  
  const int vpts = vpoints[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  
  eedot = (float *)malloc((2 + nel) * sizeof(float));
  one = 1.0;
  two = 2.0;
  
  if(visits == 0)
    {
      for(e = 1; e <= nel; e++)
	eedot[e] = one;
    }
  else
    strain_rate_2_inv(E, eedot, 1);
  
  for(e = 1; e <= nel; e++)
    {
      exponent1 = one - one / E->viscosity.sdepv_expt[E->mat[e] - 1];
      scale = pow(two * eedot[e] / E->viscosity.sdepv_trns[E->mat[e] - 1], exponent1);
      for(jj = 1; jj <= vpts; jj++)
	EEta[(e - 1) * vpts + jj] = two * EEta[(e - 1) * vpts + jj] / 
	  (one + scale * pow(EEta[(e - 1) * vpts + jj], exponent1));
    }
  
  
  visits++;
  
  free((void *)eedot);
  return;
}

/* ==================================================================================================== */
void strain_rate_2_inv(struct All_variables *E, float *EEDOT, int SQRT)
{
  double edot[4][4], dudx[4][4], rtf[4][9];
  float VV[4][9], Vxyz[9][9];
  
  //int e, i, j, p, q, n, nel, k;
  int e, i, j, p, q, n, nel;
  
  const int dims = E->mesh.nsd;
  const int ends = enodes[dims];
  const int lev = E->mesh.levmax;
  //const int nno = E->lmesh.nno;
  //const int vpts = vpoints[dims];
  const int ppts = ppoints[dims];
  //const int sphere_key = 1;
  
  nel = E->lmesh.nel;
  
  if(E->control.Rsphere)
    {
      
      for(e = 1; e <= nel; e++)
	{
	  
	  get_rtf(E, e, 2, rtf, lev);
	  
	  for(i = 1; i <= ends; i++)
	    {
	      n = E->ien[e].node[i];
	      VV[1][i] = E->V[1][n];
	      VV[2][i] = E->V[2][n];
	      VV[3][i] = E->V[3][n];
	    }
	  
	  for(j = 1; j <= ppts; j++)
	    {
	      Vxyz[1][j] = 0.0;
	      Vxyz[2][j] = 0.0;
	      Vxyz[3][j] = 0.0;
	      Vxyz[4][j] = 0.0;
	      Vxyz[5][j] = 0.0;
	      Vxyz[6][j] = 0.0;
	    }
	  
	  for(j = 1; j <= ppts; j++)
	    {
	      for(i = 1; i <= ends; i++)
		{
		  Vxyz[1][j] += (VV[1][i] * E->gNX[e].ppt[GNPXINDEX(0, i, j)] + 
				 VV[3][i] * E->N.ppt[GNPINDEX(i, j)]) * rtf[3][j];
		  Vxyz[2][j] += ((VV[2][i] * E->gNX[e].ppt[GNPXINDEX(1, i, j)] + 
				  VV[1][i] * E->N.ppt[GNPINDEX(i, j)] * cos(rtf[1][j])) / sin(rtf[1][j]) + 
				 VV[3][i] * E->N.ppt[GNPINDEX(i, j)]) * rtf[3][j];
		  Vxyz[3][j] += VV[3][i] * E->gNX[e].ppt[GNPXINDEX(2, i, j)];
		  
		  Vxyz[4][j] += ((VV[1][i] * E->gNX[e].ppt[GNPXINDEX(1, i, j)] - 
				  VV[2][i] * E->N.ppt[GNPINDEX(i, j)] * cos(rtf[1][j])) / sin(rtf[1][j]) + 
				 VV[2][i] * E->gNX[e].ppt[GNPXINDEX(0, i, j)]) * rtf[3][j];
		  Vxyz[5][j] += VV[1][i] * E->gNX[e].ppt[GNPXINDEX(2, i, j)] + rtf[3][j] * 
		    (VV[3][i] * E->gNX[e].ppt[GNPXINDEX(0, i, j)] - VV[1][i] * E->N.ppt[GNPINDEX(i, j)]);
		  Vxyz[6][j] += VV[2][i] * E->gNX[e].ppt[GNPXINDEX(2, i, j)] + rtf[3][j] * 
		    (VV[3][i] * E->gNX[e].ppt[GNPXINDEX(1, i, j)] / sin(rtf[1][j]) - VV[2][i] * E->N.ppt[GNPINDEX(i, j)]);
		}
	      edot[1][1] = 2.0 * Vxyz[1][j];
	      edot[2][2] = 2.0 * Vxyz[2][j];
	      edot[3][3] = 2.0 * Vxyz[3][j];
	      edot[1][2] = Vxyz[4][j];
	      edot[1][3] = Vxyz[5][j];
	      edot[2][3] = Vxyz[6][j];
	    }
	  
	  EEDOT[e] = edot[1][1] * edot[1][1] + edot[1][2] * edot[1][2] * 2.0 + edot[2][2] * edot[2][2] + 
	    edot[2][3] * edot[2][3] * 2.0 + edot[3][3] * edot[3][3] + edot[1][3] * edot[1][3] * 2.0;
	  
	}
      
    }
  
  else if(E->control.CART3D)
    {
      
      for(e = 1; e <= nel; e++)
	{
	  
	  for(i = 1; i <= ends; i++)
	    {
	      n = E->ien[e].node[i];
	      VV[1][i] = E->V[1][n];
	      VV[2][i] = E->V[2][n];
	      if(dims == 3)
		VV[3][i] = E->V[3][n];
	    }
	  
	  for(p = 1; p <= dims; p++)
	    for(q = 1; q <= dims; q++)
	      dudx[p][q] = 0.0;
	  
	  for(i = 1; i <= ends; i++)
	    for(p = 1; p <= dims; p++)
	      for(q = 1; q <= dims; q++)
		dudx[p][q] += VV[p][i] * E->gNX[e].ppt[GNPXINDEX(q - 1, i, 1)];
	  
	  for(p = 1; p <= dims; p++)
	    for(q = 1; q <= dims; q++)
	      edot[p][q] = dudx[p][q] + dudx[q][p];
	  
	  if(dims == 2)
	    EEDOT[e] = edot[1][1] * edot[1][1] + edot[2][2] * edot[2][2] + edot[1][2] * edot[1][2] * 2.0;
	  
	  else if(dims == 3)
	    EEDOT[e] = edot[1][1] * edot[1][1] + edot[1][2] * edot[1][2] * 2.0 + edot[2][2] * edot[2][2] + 
	      edot[2][3] * edot[2][3] * 2.0 + edot[3][3] * edot[3][3] + edot[1][3] * edot[1][3] * 2.0;
	  
	}
    }
  
  
  if(SQRT)
    for(e = 1; e <= nel; e++)
      EEDOT[e] = sqrt(0.5 * EEDOT[e]);
  else
    for(e = 1; e <= nel; e++)
      EEDOT[e] *= 0.5;
  
  return;
}


/* ==================================================================================================== */
int layers(struct All_variables *E, float x3)
{
  int llayers = 0;
  
  if(x3 >= E->viscosity.zlith)  /* note Z and R equal 1.0 at the top of the box. */
    llayers = 1; /* lithosphere */
  else if(x3 >= E->viscosity.z410)
    llayers = 2; /* asthenosphere */
  else if(x3 > E->viscosity.zlm)
    llayers = 3; /* transition zone */
  else
    llayers = 4; /* lower mantle */

  /* This is what use to be here 
     if(x3 >= E->viscosity.zlith)
     llayers = 1;
     else if(x3 < E->viscosity.zlith && x3 >= E->viscosity.zlm)
     llayers = 2;
     else if(x3 < E->viscosity.zlm)
     llayers = 3;
  */
  
  return (llayers);
}


/* ==================================================================================================== */
int weak_zones(struct All_variables *E, int node, float t_b) /* very tentative */
{
  int wweak_zones;
  
  wweak_zones = 0;
  return wweak_zones;
}

/* ==================================================================================================== */
float boundary_thickness(struct All_variables *E, float *H)
{
  float thickness;
  int i, j;
  
  thickness = 0.0;
  
  for(i = 1; i <= E->lmesh.noz; i++)
    H[i] = -H[i] / E->control.Atemp;
  
  if(E->parallel.me_loc[1] == 0)
    {
      for(i = 1; i <= E->lmesh.noz; i++)
	{
	  if(H[i] > 0.45)
	    {
	      j = i;
	      break;
	    }
	}
      thickness = E->X[3][j];
    }
  
  
  MPI_Bcast(&thickness, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  return (thickness);
}

/* =================================================================================
   Get weak_region data from file before parser is shutdown -- see weak_regions
   for details on what this is used for.
   MIB, JAN07
  ================================================================================== */

void get_weak_regions(struct All_variables *E)
{
  int i, j, jj, kk;
  int el;
  float *BW;
 	  
  const int ends=enodes[E->mesh.nsd];
  const int vpts = vpoints[E->mesh.nsd];  
    
  /* Allocate memory and initialize */
  BW = (float *)malloc((2+E->lmesh.nno)*sizeof(float));
  E->viscosity.blendweak = (float *) malloc(((E->lmesh.nel+2)*vpts)*sizeof(float));
  
  for (i = 1; i <= E->lmesh.nno; i++)
    BW[i] = 0.0;
  for (i = 1; i <= ((E->lmesh.nel+2)*vpts); i++)
    E->viscosity.blendweak[i] = 0;
  
  /* Read in blended weak region field from file */
  read_new_field(E,BW,"weakregion");
  
  /* Interpolate nodal weak region values to element integration points */
  for(el=1;el<=E->lmesh.nel;el++) {
    for(jj = 1; jj <= vpts; jj++) {
      for(kk = 1; kk <= ends; kk++) { 
	E->viscosity.blendweak[(el-1)*vpts +jj] += BW[E->ien[el].node[kk]]*E->N.vpt[GNVINDEX(kk,jj)];	  	
	//fprintf(stderr,"el %d jj %d BW %e\n",el,jj,BW[E->ien[el].node[kk]]);
      }
    }
  }   
  return;
}  
 
/* =================================================================================
   Use blendweak to smooth blend weak regions with viscosity field, by modifying the
     power of EEta by some fraction (1-blendweak).
   Note: if blendweak = 0, EEta is unchanged
         if blendweak = 1, EEta = 1 (background reference viscosity)
   Example: if EEta=10^5 and blendweak=0.8, then blended EEta = pow(10,5*0.2) = 10.

   Written by Margarete Jadamec 2006;  Modified for citcomcu by MIB, JAN07
   ================================================================================== */

void weak_regions(struct All_variables *E, float *EEta)
{
  int i, j, jj, kk;
  int el, num;
  float etatemp;
  const int vpts = vpoints[E->mesh.nsd];  

  /* Use blendweak to modify viscosity with smooth transition into weak regions */
  for(el = 1; el <= E->lmesh.nel; el ++) {
    for(jj = 1; jj <= vpts; jj++) {
      num = (el-1)*vpts + jj;    
      etatemp = pow(10,log10(EEta[num])*(1.0 - E->viscosity.blendweak[num]) );
      if (etatemp < EEta[num]) /* weak zone denotes maximum viscosity of regions */
	EEta[num] = etatemp;  /* but non-linear viscosity might already be lower */
    }
  }
 
  return;
}

/* =================================================================================
   Calculates an effective viscosity assuming linear combination strain-rate accommodated
   by two different creep mechanisms. For standard olivine model these would be  
   diffusion (newtonian) and dislocation (non-newtonian) creep.

   Also has option to include effects from chemistry/composition.

   In addition, plastic yielding reduces the viscosity where the effective viscosity 
   exceeds a yield stress.

   Other creep mechanisms (e.g., low-temp plasticity) can be added, but aren't included yet.

   Stress-dependence is written in terms of the second invariant of strain-rate tensor.

   Expects dimensional input parameters for flow laws. Uses reference viscosity from
   input file (E->data.ref_viscosity) to non-dimensionalize in last step.

   Note internal calculation is done using doubles to maintain precision of large
   viscosities (>10^32; these are usually over-written by yield stress near the end).

   Input parameters for each material group (4 layers + composition): 
   Creep Mechanism 1: 1 - 5 (e.g., diffusion creep).
   Creep Mechanism 2: 6 - 11 (e.g., dislocation creep).
   -- N,E,V,COH: for each mechanism and each layer.
      N - pre-exponential constant (Pascals; and for grain-size in meters, and COH in H/ppm-Si).
      E - activation energy (joules)
      V - activation volume (m^-3)
      COH - OH concentration (H/ppm-Si)
   -- n,r,p, : for each mechanism (d-grain-size could later on be an internal variable).
      n (sdepv_expt) - stress-dependence exponent
      r (coh_expt) - COH exponent.
      p (grainsize_expt) - grain size exponent
      d - grain size (meters).
  
   Makes no assumptions about where a particular mechanism is active... need to control
   this with input parameters.
   ================================================================================== */


void composite_viscosity(struct All_variables *E, float *Eta, float * EEta)
{
  static int visits = 0;

  /* counting */
  int i, j, k, l, el, jj, kk, node; 
  /* flow law variables */
  float n[10], n1[10], n2[10], r[10], p[10]; /* up to 9 creep mechanisms */
  double N[10], Es[10], V[10], COH[10], d[10];  /* up to 9 creep mechanisms */
  double eta[10];
  /* dimensionalization */
  double g, erad, rhor, beta, dT, R, eo;  
  /* Internal parameters */
  int nmat, emat[10], nmech, typ;
  double TT[9], ZZ[9], tempT, pres, Ttot, depthz;
  double A, sigmay, sigma;
  double etatop, etabot, eta_comp, eta_old;
  double zero, one, two;
  float *Xtmp[4], *eedot, vtyp;

  const int vpts = vpoints[E->mesh.nsd];
  const int ends = enodes[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  const int levmax = E->mesh.levmax;

  if (visits == 0) {
    fprintf(E->fp,"In composite_viscosity\n");
    for(l = 1; l <= E->viscosity.num_mat; l++)
      {
	fprintf(E->fp,"\tlayer %d/%d: N0=%g E=%g V=%g COH=%g d=%g n=%g p=%g r=%g\n", l, E->viscosity.num_mat, 
		E->viscosity.N0[l - 1], E->viscosity.E[l - 1], E->viscosity.V[l - 1], 
		E->viscosity.COH[l - 1],E->viscosity.d[l - 1],E->viscosity.sdepv_expt[l - 1],
		E->viscosity.grainsize_expt[l - 1],E->viscosity.coh_expt[l - 1]);
      }
    fprintf(E->fp,"YIELD=%d: yield_surf=%g yield_mantle=%g yield_grad=%g yield_temp=%g\n",E->viscosity.YIELD,
	    E->viscosity.yield_surf, E->viscosity.yield_mantle,E->viscosity.yield_grad,E->viscosity.yield_temp);
    fflush(E->fp);
  }

  zero = 0.0;
  one = 1.0;
  two = 2.0;

  nmech = E->viscosity.num_crpmech;    /* number of creep mechanisms */
  nmat = E->viscosity.num_mat/nmech;   /* 4 depth layers + composition */

  /* Allocate memory */
  eedot = (float *) malloc( (2+nel)*sizeof(float)); /* note constant over an element */

  /* Dimensionalization Variables */
  g = E->data.grav_acc;              /* m^2/s */
  rhor = E->data.density;            /* kg/m^3 */
  erad = E->monitor.length_scale;    /* default is in meters */
  R = E->data.gas_const;             /* J/mol*K */
  beta = E->data.adiabat_comp;       /* Pa^-1 */
  eo = E->data.therm_diff/(erad*erad);   /* per second */

  /* Internal variable pointer for coordinates */
  if(E->control.Rsphere)
    {
      for(i = 1; i <= E->mesh.nsd; i++)
	Xtmp[i] = E->SXX[levmax][i];      
    }
  else if(E->control.CART3D)
    {
      for(i = 1; i <= E->mesh.nsd; i++)
	Xtmp[i] = E->XX[levmax][i];    
    }

  /* Get second invariant strain-rate */
  if (visits==0)
    {
      for(el = 1; el <= nel; el++)
	eedot[el] = 1e-15/eo; /* gets redimensionalized further down */
    }
  else
    {
      strain_rate_2_inv(E,eedot,1); /* non-dimensional */
    }

  /* Get viscosity for each element */
  
  for(el = 1; el <= nel; el++)
    {

      for (i = 1; i <= nmech; i++ ) 
	{
	  emat[i] = E->mat[el] + (i - 1)*nmat - 1; 
	  /* Flow Law Variables  */
	  N[i] = E->viscosity.N0[emat[i]];
	  Es[i] = E->viscosity.E[emat[i]];
	  V[i] = E->viscosity.V[emat[i]];
	  COH[i] = E->viscosity.COH[emat[i]];
	  d[i] = E->viscosity.d[emat[i]];
	  p[i] = E->viscosity.grainsize_expt[emat[i]];
	  r[i] = E->viscosity.coh_expt[emat[i]];
	  n[i] = E->viscosity.sdepv_expt[emat[i]];
	  n1[i] = 1/n[i];
	  n2[i] = (1-n[i])/n[i];	 
	}
      
      for (kk = 1; kk <= ends; kk++ ) 
	{
	  TT[kk] = min(one,max(E->T[E->ien[el].node[kk]],zero));
	  ZZ[kk] = one - Xtmp[3][E->ien[el].node[kk]];	
	}
 

      /* loop for integration points */
      typ = 1;
      vtyp = zero;      
      for (jj = 1; jj <= vpts; jj++) {
	tempT = zero;
	pres = zero;
	depthz = zero;
	for (i = 1; i <= nmech; i++ ) 
	  eta[i] = zero;
	/* loop over element nodes */
	for (kk = 1; kk <= ends; kk++ )  
	{ 
	  tempT += E->data.ref_temperature*TT[kk]*E->N.vpt[GNVINDEX(kk,jj)];          /* K */
	  pres += (-1/beta)*log(1-rhor*g*beta*erad*ZZ[kk])*E->N.vpt[GNVINDEX(kk,jj)]; /* Pa */
	  depthz += erad*ZZ[kk]*E->N.vpt[GNVINDEX(kk,jj)];                            /* m */
	} /* ends */
	
	Ttot = E->data.Tsurf + E->data.dT_dz*depthz + tempT; /* K, surface + adiabatic gradient + dyn-temp */

	/* each creep mechanism */
	for (i = 1; i <= nmech; i++ ) 
	  {
	    A = pow(eo*eedot[el],n2[i])/pow( N[i]* pow(COH[i],r[i])/ pow(d[i],p[i]),n1[i]);
	    eta[i] = A*exp( (Es[i] + V[i]*pres)/(n[i]*R*Ttot) );	 
	    if (i>1) 
	      if (eta[i] < eta[i-1])
		typ = i;
	  } /* nmech */

	/* Get composite viscosity for this element-vpoint */
	etatop = 1;
	etabot = 0;
	eta_comp = 0;
	for (i = 1; i <= nmech; i++ ) 
	  etatop *= eta[i];  /* e.g. eta1*eta2*eta3 */
	for (i = 1; i <= nmech; i++ ) 
	  etabot += etatop/eta[i]; /* e.g. eta2*eta3 + eta1+eta3 + eta1*eta2 */

	eta_comp = etatop/etabot;

	if (E->viscosity.YIELD) 
	  {	   
	    if (Ttot < E->viscosity.yield_temp) 
	      {
		if (E->viscosity.yield_surf <= E->viscosity.yield_mantle) 
		  {
		    sigmay = min(E->viscosity.yield_surf + 
				 E->viscosity.yield_grad*depthz,E->viscosity.yield_mantle);
		  }
		else if (E->viscosity.yield_surf > E->viscosity.yield_mantle)
		  {
		    sigmay = max(E->viscosity.yield_surf + 
				 E->viscosity.yield_grad*depthz,E->viscosity.yield_mantle);
		  }
		sigma = eo*eedot[el]*eta_comp;
		if (sigma > sigmay )
		  {
		    eta_comp = sigmay/(eo*eedot[el]);		   
		    typ = nmech + 1;
		  }
	      } /* yield_temp */
	  } /* YIELD */

	/* ADD SPECIAL CONDITION FOR TOP ROW OF NODES WITH VELOCITY BC */
	/* need to change element loop to use elx, ely, elz, so know which elz you are */
	// if (E->control.vbcs_file || E->mech.topvbc ) {
	// if (E->parallel.me_loc[1] == 1) && (i==1)) {
	//    eta = eta[1];  
	// }
	// }

	/* Non-dimensionalize eta_comp */
	eta_comp /= E->data.ref_viscosity;
	
	/* check max value here because etacomp is double, but EEta is a float */
	if (E->viscosity.MAX) 
	  {
	    if (eta_comp > E->viscosity.max_value)
	      {
		eta_comp = E->viscosity.max_value;
	      }
	  }
	
	/* Assign eta_comp to EEta */
	if (visits == 0)
	  {
	    EEta[(el-1)*vpts + jj] = eta_comp;
	  }
	else
	  {
	    eta_old = EEta[(el-1)*vpts + jj];
	    EEta[(el-1)*vpts + jj] = 0.5*(eta_comp + eta_old);
	  }

	/* add vtyp for this element-vpt */
	vtyp += (float) typ;

      } /* vpts */

      E->visctyp[el] = vtyp/vpts; /* element average */

    } /* el */

  visits++;

  free ((void *) eedot);

  return;
}
