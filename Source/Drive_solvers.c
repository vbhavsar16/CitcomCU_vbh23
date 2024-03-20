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

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void general_stokes_solver(struct All_variables *E)
{
	//float vmag;
	double *delta_U;
	//double *force, Udot_mag, dUdot_mag;
	double Udot_mag, dUdot_mag;
	double time;
        double timestartLoop, timethisit; /* MAJ 10/28/008*/
	double timeLimit = 168.0*3600; /* 7/1/17 - Changed to 7 days (168 hrs) MAJ. 6/26/16 - Changed to be LARGE ** NEED to AUTOMATE THIS*MAJ.<> Now assumes 24 wallclock. MAJ 9/2012. Assumes 48 hour wallclock time MAJ 10/28/08*/
	//double timeLimit = 23.65*3600; /*Now assumes 24 wallclock. MAJ 9/2012. Assumes 48 hour wallclock time MAJ 10/28/08*/
	//int count, i, j, k;
	int count, i;

	static float *oldU;
	static int visits = 0;

	//const int nno = E->lmesh.nno;
	const int neq = E->lmesh.neq;
	//const int dims = E->mesh.nsd;
	//const int vpts = vpoints[E->mesh.nsd];

	if(visits == 0)
	{
		oldU = (float *)malloc((neq + 2) * sizeof(float));
		for(i = 1; i <= neq; i++)
			oldU[i] = 0.0;
		visits++;
	}

	dUdot_mag = 0.0;

	delta_U = (double *)malloc((neq + 2) * sizeof(double));

	/* FIRST store the old velocity field */

	E->monitor.elapsed_time_vsoln1 = E->monitor.elapsed_time_vsoln;
	E->monitor.elapsed_time_vsoln = E->monitor.elapsed_time;

	if(E->parallel.me == 0)
		time = CPU_time0();

	velocities_conform_bcs(E, E->U);

	assemble_forces(E, 0);

/*	
	if(E->parallel.me==0)
	{
		fprintf(stderr,"time1= %g seconds\n",CPU_time0()-time);
		time=CPU_time0();
	}
*/

	count = 1;

	do
	{
		timestartLoop = CPU_time0(); /*MAJ 10/28/08*/
		
		if (E->parallel.me == 0) { /* Added for timing on Lonestar 9/12/2012 MAJ.*/
			fprintf(stderr, "In Drive Solvers starting the do loop. cnt = %d\n", count);
		}

		if(E->viscosity.update_allowed)
			get_system_viscosity(E, 1, E->EVI[E->mesh.levmax], E->VI[E->mesh.levmax]);

		construct_stiffness_B_matrix(E);

		if (E->parallel.me == 0) { /* Added for timing on Lonestar 9/12/2012 MAJ.*/
			fprintf(stderr, "In Drive Solvers, Finished running construct_stiffness_B_matrix. cnt = %d\n", count);
		}

		solve_constrained_flow_iterative(E);

		if (E->parallel.me == 0) { /* Added for timing on Lonestar 9/12/2012 MAJ.*/                         
			fprintf(stderr, "In Drive Solvers, Finished running solve_constrained_flow_iterative. cnt = %d\n", count);
		}    

		if(E->viscosity.SDEPV || E->viscosity.TSCDEPV) /* MIB, JAN07 added TSCDEPV */
		{
			for(i = 1; i <= neq; i++)
			{
				delta_U[i] = E->U[i] - oldU[i];
				oldU[i] = E->U[i];
			}
			Udot_mag = sqrt(global_vdot(E, E->U, E->U, E->mesh.levmax));
			dUdot_mag = sqrt(global_vdot(E, delta_U, delta_U, E->mesh.levmax));

			if(Udot_mag != 0.0)
				dUdot_mag /= Udot_mag;

			//	if(E->control.sdepv_print_convergence < E->monitor.solution_cycles && E->parallel.me == 0)
			if(E->control.sdepv_print_convergence && E->parallel.me == 0)
			{
				fprintf(stderr, "Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d. And sdepv_misfit = %0.4e.\n", dUdot_mag, Udot_mag, count, E->viscosity.sdepv_misfit);
				fprintf(E->fp, "Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d. And sdepv_misfit = %0.4e.\n", dUdot_mag, Udot_mag, count, E->viscosity.sdepv_misfit);
				fflush(E->fp);
			}
			/* Output newtonian solution for zeroeth time-step only */
			if (E->control.output_sdepvnewt_stepzero && (E->monitor.solution_cycles == 0) && (count == 1)) {
				if (E->parallel.me == 0)
					fprintf(stderr,"Output newtonian solution for step zero\n");
				process_temp_field(E,E->monitor.solution_cycles-1);
				process_new_velocity(E,E->monitor.solution_cycles-1);
			}
			count++;
			timethisit = CPU_time0() - timestartLoop; /*MAJ 10/28/08*/
			if (E->parallel.me == 0) { /*MAJ 10/28/08. Appended for timing on Lonestar 9/10/2012 MAJ.*/
				fprintf(stderr, "count = %d. CPUtime = %0.4e. timestartLoop = %0.4e. begin_time = %0.4e. (CPUtime - begin time) = %0.4e. timethisit = %0.4e (%0.4e). timeLimit = %0.4e. timeLimit-timethisit = %0.4e. \n", count, (CPU_time0()/3600), (timestartLoop/3600), (E->monitor.begin_time/3600), ((CPU_time0() - E->monitor.begin_time)/3600), (timethisit/3600), ((CPU_time0() - timestartLoop)/3600), (timeLimit/3600), ((timeLimit-timethisit)/3600));
				fprintf(E->fp, "count = %d. CPUtime = %0.4e. timestartLoop = %0.4e. begin_time = %0.4e. (CPUtime - begin time) = %0.4e. timethisit = %0.4e (%0.4e). timeLimit = %0.4e. timeLimit-timethisit = %0.4e. \n", count, (CPU_time0()/3600), (timestartLoop/3600), (E->monitor.begin_time/3600), ((CPU_time0() - E->monitor.begin_time)/3600), (timethisit/3600), ((CPU_time0() - timestartLoop)/3600), (timeLimit/3600), ((timeLimit-timethisit)/3600));
				fflush(E->fp);
				}
		}						/* end for SDEPV   */

	} while((count < 50) && ((CPU_time0() - E->monitor.begin_time + (2.5*timethisit)) < timeLimit) && (dUdot_mag > E->viscosity.sdepv_misfit) && (E->viscosity.SDEPV || E->viscosity.TSCDEPV)); /*Added begin_time for timing on Lonestar MAJ 9/9/12. MAJ 10/28/08 timeLimit part. MIB, JAN07 added TSCDEPV */
	/*} while((count < 50) && ((CPU_time0() - E->monitor.begin_time) < (timeLimit - timethisit)) && (dUdot_mag > E->viscosity.sdepv_misfit) && (E->viscosity.SDEPV || E->viscosity.TSCDEPV));*/ /*Added begin_time for timing on Lonestar MAJ 9/9/12. MAJ 10/28/08 timeLimit part. MIB, JAN07 added TSCDEPV */

	free((void *)delta_U);

	return;
}
