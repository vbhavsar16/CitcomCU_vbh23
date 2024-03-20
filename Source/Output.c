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

/* Routine to process the output of the finite element cycles 
   and to turn them into a coherent suite  files  */


#include <fcntl.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>				/* for "system" command */
#ifndef __sunos__				/* string manipulations */
#include <strings.h>
#else
#include <string.h>
#endif

#include "element_definitions.h"
#include "global_defs.h"

/* ==========================================================================  
   Outputs model results at each node in local mesh, for record_every times.
   Everything, but temperature and coordinates.
   Modified from CitcomT, MIB JAN07.
   ========================================================================== */
void output_velo_related_binary(struct All_variables *E, int file_number)
{
	int i, j, k;
	int size2;
	char output_file[255];
	static float *SV, *visc, *pres, *vtyp, *edot;
	static int been_here = 0;
	
	if (been_here==0)
	{
		SV = (float *) malloc((3*E->lmesh.nnov+1)*sizeof(float));
		visc = (float *) malloc((E->lmesh.nno+1)*sizeof(float));
		pres = (float *) malloc((E->lmesh.nno+1)*sizeof(float));
		vtyp = (float *) malloc((E->lmesh.nno+1)*sizeof(float));
		edot = (float *) malloc((E->lmesh.nno+1)*sizeof(float));
		for (i = 1; i<=E->lmesh.nno; i++){
		  visc[i] = 0.0;
		  pres[i] = 0.0;
		  vtyp[i] = 0.0;
		  edot[i] = 0.0;
		}
		for (i = 1; i<=E->lmesh.nnov; i++){
		  SV[i] = 0.0;
		}
	}

	sprintf(output_file,"%s.velo.%d.%d",E->control.data_file,E->parallel.me,file_number);
	E->filed[8]=fopen(output_file,"w");
	sprintf(output_file,"%s.visc.%d.%d",E->control.data_file,E->parallel.me,file_number);
	E->filed[9]=fopen(output_file,"w");
	sprintf(output_file,"%s.pres.%d.%d",E->control.data_file,E->parallel.me,file_number);
	E->filed[10]=fopen(output_file,"w");
	sprintf(output_file,"%s.dstress.%d.%d",E->control.data_file,E->parallel.me,file_number);
	E->filed[11]=fopen(output_file,"w");		
	sprintf(output_file,"%s.edot.%d.%d",E->control.data_file,E->parallel.me,file_number);
	E->filed[12]=fopen(output_file,"w");
	sprintf(output_file,"%s.vtyp.%d.%d",E->control.data_file,E->parallel.me,file_number);
	E->filed[13]=fopen(output_file,"w");
	if(E->control.composition) {
		sprintf(output_file,"%s.comp.%d.%d",E->control.data_file,E->parallel.me,file_number);
		E->filed[14]=fopen(output_file,"w");
	}		
	sprintf(output_file,"%s.topo_t.%d.%d",E->control.data_file,E->parallel.me,file_number);
	E->filed[15]=fopen(output_file,"w");
	sprintf(output_file,"%s.topo_b.%d.%d",E->control.data_file,E->parallel.me,file_number);
	E->filed[16]=fopen(output_file,"w");
	//sprintf(output_file,"%s.fas410.%d.%d",E->control.data_file,E->parallel.me,file_number);
	//E->filed[17]=fopen(output_file,"w");
	//sprintf(output_file,"%s.fas670.%d.%d",E->control.data_file,E->parallel.me,file_number);
	//E->filed[18]=fopen(output_file,"w");
		
	size2 = (3*E->lmesh.nnov+1)*sizeof(float);
	get_volume_velo(E,SV);
	fwrite(SV,size2,1,E->filed[8]);
	
	size2 = (6*E->lmesh.nno+1)*sizeof(float);
	fwrite(E->dstress,size2,1,E->filed[11]);
		
	size2 = (E->lmesh.nno+1)*sizeof(float);
	get_visc_pres(E,visc,pres,edot,vtyp);
	fwrite(visc,size2,1,E->filed[9]);
	fwrite(pres,size2,1,E->filed[10]);
	fwrite(edot,size2,1,E->filed[12]);
	fwrite(vtyp,size2,1,E->filed[13]);

	if (E->control.composition) {
	   size2 = (E->lmesh.nno+1)*sizeof(float);
	   fwrite(E->C,size2,1,E->filed[14]);
	   fclose(E->filed[14]);
	}
	/* Observables on Surfaces */
	size2 = (( (E->parallel.me_loc[3]==0||E->parallel.me_loc[3]==
              E->parallel.nprocz-1)?E->lmesh.nsf:0 ) +1)*sizeof(float);
	fwrite(E->slice.tpg,size2,1,E->filed[15]);
	fwrite(E->slice.tpgb,size2,1,E->filed[16]);   
	//fwrite(E->filed[17],E->Fas410_b,size2);
	//fwrite(E->filed[18],E->Fas670_b,size2);

	fclose(E->filed[8]);
	fclose(E->filed[9]);
	fclose(E->filed[10]);
	fclose(E->filed[11]);
	fclose(E->filed[12]);
	fclose(E->filed[13]);
	fclose(E->filed[15]);
	fclose(E->filed[16]);
	//fclose(E->filed[17]);
	//fclose(E->filed[18]);
	

	return;
}

/* ==========================================================================  
   Reorders velocity as triplets at each node.
   Modified from CitcomT, MIB JAN07.
   ========================================================================== */
void get_volume_velo(struct All_variables *E, float *SV)
{
  int node;
  SV[0] = 0.0;

  for (node=1;node<=E->lmesh.nnov;node++)  {
      SV[(node-1)*3+1] = E->V[1][node];
      SV[(node-1)*3+2] = E->V[2][node];
      SV[(node-1)*3+3] = E->V[3][node];
  }

  return;
  }
  
/* ==========================================================================   
   Interpolates element values to nodes for pressure, viscosity and vtyp.  
   Modified from CitcomT, MIB JAN07 
   ========================================================================== */
void get_visc_pres(struct All_variables *E, float *visc, float *pres, float *edot, float *vtyp)
{
  
  int i,j,k,el,node,n;
  float ele_v;
  float *eedot;
 
  const int dims=E->mesh.nsd;
  const int nel=E->lmesh.nel;
  const int vpts=vpoints[E->mesh.nsd];
  const int lev = E->mesh.levmax;
  const int ends=enodes[dims];

  /* First need to get strain-rate invariant at elements */
  eedot = (float *)malloc((2 + nel) * sizeof(float));
  strain_rate_2_inv(E, eedot, 1);

  for (el=1;el<=E->lmesh.nel;el++)   {
    ele_v = 0.0;
    for (j=1;j<=vpts;j++) 
      ele_v += E->EVI[lev][(el-1)*vpts+j];

    ele_v /= vpts;
    for(j=1;j<=ends;j++) {
      node = E->ien[el].node[j];
      visc[node] += E->tw[node] * ele_v;
    }
  }
 
  
  for (el=1;el<=E->lmesh.nel;el++)
    for(j=1;j<=ends;j++)   {
      node = E->ien[el].node[j];
      pres[node] += E->tw[node] * E->P[el];
      edot[node] += E->tw[node] * eedot[el];
      vtyp[node] += E->tw[node] * E->visctyp[el];
    }


  exchange_node_f20(E,visc,lev);
  exchange_node_f20(E,pres,lev);
  exchange_node_f20(E,edot,lev);
  exchange_node_f20(E,vtyp,lev);
  
  return;
}

/* ==========================================================================  
   Output coordinates for first time step. Also outputs temperature at 
   record-every time-steps.
   Modified from CitcomT, MIB JAN07.
   ========================================================================== */

void output_temp_binary(struct All_variables *E, int file_number)
{
	int i, j;
	static int been_here = 0;
	static int size2, size1;
	char output_file[255];

	if (been_here==0)
	{
		
		/* First: write out mesh header information */
		sprintf(output_file,"%s.hdr",E->control.data_file);
		if ((E->filed[1]=fopen(output_file,"w")) == NULL){
			fprintf(stderr,"Unable to open mesh header file %s\n",output_file);
		}	
		else {
		  if (E->control.NMULTIGRID || E->control.EMULTIGRID ) {
		    fprintf(E->filed[1],"multigrid\n");
		    fprintf(E->filed[1],"%d %d %d %d\n",E->mesh.mgunitx,E->mesh.mgunity,
			    E->mesh.mgunitz,E->mesh.levels);
		    fprintf(E->filed[1],"%d %d %d\n",E->parallel.nprocx,
			    E->parallel.nprocy,E->parallel.nprocz);
		    fclose(E->filed[1]);
		  }
		  else {
		    fprintf(E->filed[1],"conj-grad\n");
		    fprintf(E->filed[1],"%d %d %d\n",E->mesh.nox,E->mesh.noy,E->mesh.noz);
		    fprintf(E->filed[1],"%d %d %d\n",E->parallel.nprocx,
			    E->parallel.nprocy,E->parallel.nprocz);
		    fclose(E->filed[1]);
		  }
		}
			
		/* Second: write out mesh once, does not depend on time-step */
		sprintf(output_file,"%s.x.%d",E->control.data_file,E->parallel.me);
		E->filed[2]=fopen(output_file,"w");
		sprintf(output_file,"%s.y.%d",E->control.data_file,E->parallel.me);
		E->filed[3]=fopen(output_file,"w");
		sprintf(output_file,"%s.z.%d",E->control.data_file,E->parallel.me);
		E->filed[4]=fopen(output_file,"w");

		size2 = (E->lmesh.nno+1)*sizeof(float);
		
		if (E->control.CART3D) {
			fwrite(E->X[1],size2,1,E->filed[2]);
			fwrite(E->X[2],size2,1,E->filed[3]);
			fwrite(E->X[3],size2,1,E->filed[4]);
		}
		else {
			fwrite(E->SX[1],size2,1,E->filed[2]);
			fwrite(E->SX[2],size2,1,E->filed[3]);
			fwrite(E->SX[3],size2,1,E->filed[4]);
		}
		fclose(E->filed[2]);
		fclose(E->filed[3]);
		fclose(E->filed[4]);

		been_here++;
	}
	
	/* Third: write out temperature field every time-step */	
		sprintf(output_file,"%s.temp.%d.%d",E->control.data_file,
			E->parallel.me,file_number);
		E->filed[5]=fopen(output_file,"w");
		fwrite(E->T,size2,1,E->filed[5]);
		fclose(E->filed[5]);

	return;
}


/* ====================================================================== */

void process_restart(struct All_variables *E)
{
	//int fileid[20];
	//int i, j, k, ii, size2;
	//char output_file[255], in_file[255];
/*
	sprintf(output_file,"%s.temp",E->control.data_file);
	E->filed[5]=open(output_file,O_RDWR | O_CREAT, 0644);
	sprintf(output_file,"%s.velo",E->control.data_file);
	E->filed[10]=open(output_file,O_RDWR | O_CREAT, 0644);
	sprintf(output_file,"%s.topo_t",E->control.data_file);
	E->filed[11]=open(output_file,O_RDWR | O_CREAT, 0644);
	sprintf(output_file,"%s.topo_b",E->control.data_file);
	E->filed[12]=open(output_file,O_RDWR | O_CREAT, 0644);
	sprintf(output_file,"%s.visc",E->control.data_file);
	E->filed[13]=open(output_file,O_RDWR | O_CREAT, 0644);
	sprintf(output_file,"%s.fas670",E->control.data_file);
	E->filed[14]=open(output_file,O_RDWR | O_CREAT, 0644);
	sprintf(output_file,"%s.stress",E->control.data_file);
	E->filed[9]=open(output_file,O_RDWR | O_CREAT, 0644);

	sprintf(in_file,"%s.temp",E->control.data_file1);
	fileid[5]=open(in_file,O_RDONLY,0);
	sprintf(in_file,"%s.velo",E->control.data_file1);
	fileid[10]=open(in_file,O_RDONLY,0);
	sprintf(in_file,"%s.topo_t",E->control.data_file1);
	fileid[11]=open(in_file,O_RDONLY,0);
	sprintf(in_file,"%s.topo_b",E->control.data_file1);
	fileid[12]=open(in_file,O_RDONLY,0);
	sprintf(in_file,"%s.visc",E->control.data_file1);
	fileid[13]=open(in_file,O_RDONLY,0);
	sprintf(in_file,"%s.fas670",E->control.data_file1);
	fileid[14]=open(in_file,O_RDONLY,0);
	sprintf(in_file,"%s.stress",E->control.data_file1);
	fileid[9]=open(in_file,O_RDONLY,0);


	for (i=1;i<=E->control.restart_frame;i++)
	{
		ii = E->mesh.nsf;
		size2 = (2*ii+2)*sizeof(float);
		read(fileid[10],E->T,size2);
		write(E->filed[10],E->T,size2);

		size2 = (E->mesh.nsf+1)*sizeof(float);
		read(fileid[11],E->T,size2);
		write(E->filed[11],E->T,size2);
		read(fileid[12],E->T,size2);
		write(E->filed[12],E->T,size2);

		size2 = (E->mesh.nel+1)*sizeof(float);
		read(fileid[13],E->T,size2);
		write(E->filed[13],E->T,size2);

		size2 = (E->mesh.nsf+1)*sizeof(float);
		read(fileid[14],E->T,size2);
		write(E->filed[14],E->T,size2);

		size2 = (12*E->mesh.nsf+1)*sizeof(float);
		read(fileid[9],E->T,size2);
		write(E->filed[9],E->T,size2);

		size2= (E->mesh.nno+1)*sizeof(float);
		read(fileid[5],E->T,size2);
		write(E->filed[5],E->T,size2);
	}

	close(fileid[5]);
	close(fileid[10]);
	close(fileid[11]);
	close(fileid[12]);
	close(fileid[13]);
	close(fileid[14]);
	close(fileid[9]);
*/
	return;
}

/*   ======================================================================
    ======================================================================  */

void print_field_spectral_regular(struct All_variables *E, float *TG, float *sphc, float *sphs, int proc_loc, char *filen)
{
	FILE *fp, *fp1;
	char output_file[255];
	int i, node, j, ll, mm;
	float minx, maxx, t, f, rad;
	rad = 180.0 / M_PI;

	maxx = -1.e26;
	minx = 1.e26;
	if(E->parallel.me == proc_loc)
	{
		sprintf(output_file, "%s.%s_intp", E->control.data_file, filen);
		fp = fopen(output_file, "w");
		for(i = E->sphere.nox; i >= 1; i--)
			for(j = 1; j <= E->sphere.noy; j++)
			{
				node = i + (j - 1) * E->sphere.nox;
				t = 90 - E->sphere.sx[1][node] * rad;
				f = E->sphere.sx[2][node] * rad;
				fprintf(fp, "%.3e %.3e %.4e\n", f, t, TG[node]);
				if(TG[node] > maxx)
					maxx = TG[node];
				if(TG[node] < minx)
					minx = TG[node];
			}
		fprintf(stderr, "lmaxx=%.4e lminx=%.4e for %s\n", maxx, minx, filen);
		fprintf(E->fp, "lmaxx=%.4e lminx=%.4e for %s\n", maxx, minx, filen);
		fclose(fp);

		sprintf(output_file, "%s.%s_sharm", E->control.data_file, filen);
		fp1 = fopen(output_file, "w");
		fprintf(fp1, "lmaxx=%.4e lminx=%.4e for %s\n", maxx, minx, filen);
		fprintf(fp1, " ll   mm     cos      sin \n");
		for(ll = 0; ll <= E->sphere.output_llmax; ll++)
			for(mm = 0; mm <= ll; mm++)
			{
				i = E->sphere.hindex[ll][mm];
				fprintf(fp1, "%3d %3d %.4e %.4e \n", ll, mm, sphc[i], sphs[i]);
			}

		fclose(fp1);
	}

	return;
}

/*   ======================================================================
    ======================================================================  */
void output_velo_related(struct All_variables *E, int file_number)
{
	
	int i;	
	char output_file[255];	
	static int been_here = 0;
	

	if(been_here == 0)
	{
		sprintf(output_file, "%s.coord.%d", E->control.data_file2, E->parallel.me);
		E->filed[13] = fopen(output_file, "w");
		been_here++;
		fprintf(E->filed[13], "%6d\n", E->lmesh.nno);
		if(E->control.CART3D)
			for(i = 1; i <= E->lmesh.nno; i++)
				fprintf(E->filed[13], "%.5e %.5e %.5e\n", E->X[1][i], E->X[2][i], E->X[3][i]);
		else
			for(i = 1; i <= E->lmesh.nno; i++)
				fprintf(E->filed[13], "%.5e %.5e %.5e\n", E->SX[1][i], E->SX[2][i], E->SX[3][i]);
		fclose(E->filed[13]);
	}

	if(E->parallel.me < E->parallel.nprocz)
	{
		sprintf(output_file, "%s.ave.%d.%d", E->control.data_file2, E->parallel.me, file_number);
		E->filed[10] = fopen(output_file, "w");

		fprintf(E->filed[10], "%6d %6d %.5e %.5e %.5e %.4e %.4e %.5e %.5e\n", E->lmesh.noz, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nut, E->slice.Nub, E->data.T_adi0, E->data.T_adi1, E->monitor.Sigma_interior, E->monitor.Sigma_max);
		for(i = 1; i <= E->lmesh.noz; i++)
			fprintf(E->filed[10], "%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", E->Have.T[i], E->Have.vrms[i], E->Have.Vi[i], E->Have.Rho[i], E->Have.F[i], E->Have.f[i], E->Have.C[i], E->Have.Tadi[i]);
		fclose(E->filed[10]);
	}


	if(file_number % (10 * E->control.record_every) == 0)
	{
		sprintf(output_file, "%s.temp.%d.%d", E->control.data_file2, E->parallel.me, file_number);
		E->filed[10] = fopen(output_file, "w");
		fprintf(E->filed[10], "%6d %6d %.5e\n", E->lmesh.nno, E->advection.timesteps, E->monitor.elapsed_time);
		if(file_number % (20 * E->control.record_every) == 0)
		{
			if(E->control.composition == 0)
				for(i = 1; i <= E->lmesh.nno; i++)
//        fprintf(E->filed[10],"%.5e %.4e %.4e %.5e %.5e %.5e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->V[1][i],E->V[2][i],E->V[3][i]);
//       fprintf(E->filed[10],"%.5e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
					fprintf(E->filed[10], "%.5e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
			else if(E->control.composition)
				for(i = 1; i <= E->lmesh.nno; i++)
//      fprintf(E->filed[10],"%.5e %.4e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
					fprintf(E->filed[10], "%.5e %.4e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
		}
		else
		{
			if(E->control.composition == 0)
				for(i = 1; i <= E->lmesh.nno; i++)
//      fprintf(E->filed[10],"%.5e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
					fprintf(E->filed[10], "%.5e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
			else if(E->control.composition)
				for(i = 1; i <= E->lmesh.nno; i++)
//      fprintf(E->filed[10],"%.5e %.4e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
					fprintf(E->filed[10], "%.5e %.4e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
		}
		fclose(E->filed[10]);

		if(E->parallel.me_loc[3] == E->parallel.nprocz - 1)
		{
			sprintf(output_file, "%s.th_t.%d.%d", E->control.data_file2, E->parallel.me, file_number);
			E->filed[11] = fopen(output_file, "w");
			fprintf(E->filed[11], "%6d %6d %.5e %.5e\n", E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nut);
			for(i = 1; i <= E->lmesh.nsf; i++)
			{
				fprintf(E->filed[11], "%.5e %.5e %.5e %.5e\n", E->slice.tpg[i], E->slice.shflux[i], E->Fas410_b[i], E->Fas670_b[i]);
			}
			fclose(E->filed[11]);
		}
		if(E->parallel.me_loc[3] == 0)
		{
			sprintf(output_file, "%s.th_b.%d.%d", E->control.data_file2, E->parallel.me, file_number);
			E->filed[11] = fopen(output_file, "w");
			fprintf(E->filed[11], "%6d %6d %.5e %.5e\n", E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nub);
			for(i = 1; i <= E->lmesh.nsf; i++)
			{
				fprintf(E->filed[11], "%.5e %.5e %.5e %.5e\n", E->slice.tpgb[i], E->slice.bhflux[i], E->Fas410_b[i], E->Fas670_b[i]);
			}
			fclose(E->filed[11]);
		}

	}

	if(E->control.composition && file_number % (10 * E->control.record_every) == 0)
	{
		sprintf(output_file, "%s.traces.%d", E->control.data_file2, E->parallel.me);
		E->filed[10] = fopen(output_file, "w");
		fprintf(E->filed[10], "%6d %6d %.5e\n", E->advection.markers, E->advection.timesteps, E->monitor.elapsed_time);
		for(i = 1; i <= E->advection.markers; i++)
			fprintf(E->filed[10], "%g %g %g %d %d\n", E->XMC[1][i], E->XMC[2][i], E->XMC[3][i], E->CElement[i], E->C12[i]);
		for(i = 1; i <= E->lmesh.nel; i++)
			fprintf(E->filed[10], "%g\n", E->CE[i]);
		fclose(E->filed[10]);
	}



	return;
}

