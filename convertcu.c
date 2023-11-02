#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>


#define max(A,B) (((A) > (B)) ? (A) : (B))


main(argc,argv)
     int argc;
     char **argv;
{
  char *prefix, *timestep, *dimen;
  char hdrfile[50], dimenfile[50], infile[50], outfile[50], solver[50];
  FILE *fp;
  int d, i, j, k, m, ig, jg, kg, me;
  int nodel, nodeg, nodesl, nodesg, nsd, nsf, nsfl;
  int mgunitx, mgunity, mgunitz, levels, nprocx, nprocy, nprocz;
  int nprocxz, nprocxy, nproczy, nproc, levmax;
  int nox, noy, noz, nno, elx, ely, elz, noxl, noyl, nozl, nnol;
  int *locx, *locy, *locz;
  int size1, size2, size3, size4;
  int size1g, size2g, size3g, size4g;
  int fcomp;
  float *X, *Y, *Z, *TEMP, *VISC, *PRES, *VEL, *VTYP, *COMP, *EDOT, *STRESS;
  float *TPGT, *TPGB;
  float *x, *y, *z, *temp, *visc, *pres, *vel, *vtyp, *comp, *edot, *stress;
  float *tpgt, *tpgb;

  float Ra, eta, kappa, rhoref, alpha, dT, To, g, erad, rhotop, rhobot, rhochem;
  float sec2yrs, r2d, Ra1, econ, pcon, vcon, drhotop, drhobot, drhochem;

  /* Dimensionalization Constants */
  sec2yrs = 3600*24*365.25;               /* seconds per year */
  fprintf(stderr,"sec2yrs %g\n",sec2yrs);
  r2d = 180.0/3.14157;
    
  nsd = 3;
  /* Read Parameter File */  
  if (argc < 3)   {
      fprintf(stderr,"Usage: executable prefix time-step [dimen-file] \n");
   exit(10);
   }
  
   prefix = argv[1];
   timestep = argv[2];
   dimen = argv[3];
  
   fprintf(stderr,"run: %s, time-step %s\n",prefix,timestep);

   sprintf(hdrfile,"%s.hdr",prefix);
   if ((fp = fopen(hdrfile,"r")) == NULL)  {
      fprintf(stderr,"File: %s is unreadable\n",argv[1]);
      exit(11);
   }
   fscanf(fp,"%s",&(solver));
   fprintf(stderr,"solver %s\n",solver);
   if (!(strcmp(solver,"multigrid")) ) {
     fscanf(fp,"%d %d %d %d",&(mgunitx),&(mgunity),&(mgunitz),&(levels));
     fscanf(fp,"%d %d %d",&(nprocx),&(nprocy),&(nprocz));
     fclose(fp);
     
     levmax = levels - 1;
     noz = mgunitz*pow(2,levmax) + 1;
     nox = mgunitx*pow(2,levmax) + 1;
     noy = mgunity*pow(2,levmax) + 1;
     fprintf(stderr,"%d %d %d %d\n",mgunitx, mgunity, mgunitz, levels);

   }
   else {
     fscanf(fp,"%d %d %d",&(nox),&(noy),&(noz));
     fscanf(fp,"%d %d %d",&(nprocx),&(nprocy),&(nprocz));
     fclose(fp);  
   }  
   fprintf(stderr,"mesh: %d %d %d\n",nox, noy, noz);
   fprintf(stderr,"proc: %d %d %d\n",nprocx, nprocy, nprocz);

   /* Read in dimensionalization values if given as input */
   if (argc > 3) {
     /* To use default values */
     sprintf(dimenfile,"%s",dimen);
     if (!(strcmp(dimenfile,"default"))) {
       fprintf(stderr,"Using default dimensionalization\n");
       g = 9.81;          /* m/s^2 */
       eta = 1e20;        /* Pa s */
       kappa = 1e-6;      /* m^2/s */
       rhoref = 3300;     /* kg/m^3 */
       alpha = 2e-5;      /* K^{-1} */
       dT = 1400;
       erad = 6371.137e3; /* m */
       Ra = g*rhoref*alpha*dT*erad*erad*erad/(eta*kappa);

       rhotop = 2300;     /* kg/m^3 */
       rhobot = 7700;     /* kg/m^3 */
       rhochem = 2800;    /* kg/m^3 */
       To = 0;            /* K */
     }
     else { /* read from file */
       if ((fp = fopen(dimenfile,"r")) == NULL)  {
	 fprintf(stderr,"File: %s is unreadable\n",dimenfile);
	 exit(11);
       }
       fprintf(stderr,"Using dimensionalization from input file %s\n",dimenfile);
       fscanf(fp,"%f",&(g));
       fscanf(fp,"%f",&(rhoref));
       fscanf(fp,"%f",&(alpha));
       fscanf(fp,"%f",&(dT));
       fscanf(fp,"%f",&(erad));
       fscanf(fp,"%f",&(eta));
       fscanf(fp,"%f",&(kappa));
       fscanf(fp,"%f",&(Ra1));      
       fscanf(fp,"%f",&(rhotop));
       fscanf(fp,"%f",&(rhobot));
       fscanf(fp,"%f",&(rhochem));
       fscanf(fp,"%f",&(To));

       Ra = g*rhoref*alpha*dT*erad*erad*erad/(eta*kappa);
       fprintf(stderr,"Using internally calculate Ra if input & internal do not agree\n");
       fprintf(stderr,"Check: Ra %g Ra1 %g \n", Ra, Ra1);
     }
     /* Dimensionalization Variables */
     econ = kappa/(erad*erad);             /* s^-1 */
     pcon = (eta*kappa)/(erad*erad);       /* kg/(m*s^2) which is Pa  */
     vcon = kappa*100*sec2yrs/erad;      /* CM/YR */
     drhotop = rhoref - rhotop;
     drhobot = rhobot-rhoref;
     drhochem = rhoref - rhochem;
     fprintf(stderr,"econ %g pcon %g vcon %g drhotop %g drhobot %g drhochem %g\n",
	     econ,pcon,vcon,drhotop,drhobot,drhochem); 
     
   }

   /* Calculate Mesh Parameters */
   nprocxz = nprocx*nprocz;
   nprocxy = nprocx*nprocy;
   nproczy = nprocz*nprocy;
   nproc = nprocx*nprocy*nprocz;

   /* Global Mesh */
   nno = nox*noy*noz;
   nsf = nox*noy;
   elz = noz - 1;
   elx = nox - 1;
   ely = noy - 1;

   /* Local Mesh */
   nozl = elz/nprocz + 1;
   noxl = elx/nprocx + 1;
   noyl = ely/nprocy + 1;
   nnol = nozl*noxl*noyl;
   nsfl = noxl*noyl;
	
   //   fprintf(stderr,"nox %d noy %d noz %d nno %d nsf %d \n",
   //	  nox,noy,noz,nno,nsf);
   // fprintf(stderr,"noxl %d noyl %d nozl %d nnol %d nsfl %d\n",
   //	  noxl,noyl,nozl,nnol,nsfl);
	
   /* Allocate memory */
   locx = (int *) malloc((nno+1)*sizeof(int));
   locy = (int *) malloc((nno+1)*sizeof(int));
   locz = (int *) malloc((nno+1)*sizeof(int));
   /* Local */
   x = (float *) malloc((nnol+1)*sizeof(float));
   y = (float *) malloc((nnol+1)*sizeof(float));
   z = (float *) malloc((nnol+1)*sizeof(float));
   temp = (float *) malloc((nnol+1)*sizeof(float));
   pres  = (float *) malloc((nnol+1)*sizeof(float));
   visc = (float *) malloc((nnol+1)*sizeof(float));
   vel = (float *) malloc((3*nnol+1)*sizeof(float));
   vtyp  = (float *) malloc((nnol+1)*sizeof(float));
   comp  = (float *) malloc((nnol+1)*sizeof(float));
   edot  = (float *) malloc((nnol+1)*sizeof(float));
   stress = (float *) malloc((6*nno+1)*sizeof(float));
   tpgt = (float *) malloc((nsfl+1)*sizeof(float));
   tpgb = (float *) malloc((nsfl+1)*sizeof(float));
   /* Global */  
   X = (float *) malloc((nno+1)*sizeof(float));
   Y = (float *) malloc((nno+1)*sizeof(float));
   Z = (float *) malloc((nno+1)*sizeof(float));
   TEMP = (float *) malloc((nno+1)*sizeof(float));
   PRES  = (float *) malloc((nno+1)*sizeof(float));
   VISC = (float *) malloc((nno+1)*sizeof(float));
   VEL = (float *) malloc((3*nno+1)*sizeof(float));
   VTYP  = (float *) malloc((nno+1)*sizeof(float));
   COMP  = (float *) malloc((nno+1)*sizeof(float));
   EDOT  = (float *) malloc((nno+1)*sizeof(float));
   STRESS = (float *) malloc((6*nno+1)*sizeof(float));
   TPGT = (float *) malloc((nsf+1)*sizeof(float));
   TPGB = (float *) malloc((nsf+1)*sizeof(float));
 
	
   /* Calculate local processor information */
 
   for (m=1; m<=nproc; m++) {
    	me = m - 1;
    	locz[m] = me % nprocz;
    	locy[m] = (me + 1) / nprocxz - (((me + 1) % nprocxz == 0) ? 1 : 0);
        locx[m] = (me + 1 - locy[m] * nprocxz) / nprocz - 
        	(((me + 1 - locy[m] * nprocxz) % nprocz == 0) ? 1 : 0);
	//  fprintf(stderr,"m %d locx %d locy %d locz %d\n",m,locx[m],locy[m],locz[m]);
    }

   /* Read in processor information and assign to global mesh */
   size1 = (nnol + 1)*sizeof(float);
   size2 = (3*nnol + 1)*sizeof(float);
   size3 = (6*nnol + 1)*sizeof(float);
   size4 = (nsfl + 1)*sizeof(float);
   //  fprintf(stderr,"size1 %d size2 %d size3 %d size4 %d\n",
   //	   size1,size2,size3,size4);

    for (m=1; m<=nproc; m++) {
	me = m - 1;
	fprintf(stderr,"Reading in processor, m =  %d\n",m);		
      	/* Read in files for this processor */      
       	sprintf(infile,"%s.x.%d",prefix,me);
       	fp = fopen(infile,"r");
       	fread(x,size1,1,fp);
       	fclose(fp);
	
       	sprintf(infile,"%s.y.%d",prefix,me);
       	fp = fopen(infile,"r");
	fread(y,size1,1,fp);
       	fclose(fp);
      
       	sprintf(infile,"%s.z.%d",prefix,me);
       	fp = fopen(infile,"r");
       	fread(z,size1,1,fp);
       	fclose(fp);

       	sprintf(infile,"%s.temp.%d.%s",prefix,me,timestep);
       	fp = fopen(infile,"r");
       	fread(temp,size1,1,fp);
       	fclose(fp);    
       	sprintf(infile,"%s.pres.%d.%s",prefix,me,timestep); 
       	fp = fopen(infile,"r"); 
       	fread(pres,size1,1,fp); 
       	fclose(fp); 
       	sprintf(infile,"%s.velo.%d.%s",prefix,me,timestep); 
       	fp = fopen(infile,"r"); 
       	fread(vel,size2,1,fp); 
       	fclose(fp); 
       	sprintf(infile,"%s.visc.%d.%s",prefix,me,timestep); 
       	fp = fopen(infile,"r"); 
       	fread(visc,size1,1,fp);
       	fclose(fp);     
	sprintf(infile,"%s.comp.%d.%s",prefix,me,timestep);
       	if ((fp = fopen(infile,"r")) == NULL) {
	  fcomp = 0;
	}
	else {
	  fcomp = 1;
	  fread(comp,size1,1,fp);
	  fclose(fp);
	}
	
       	sprintf(infile,"%s.vtyp.%d.%s",prefix,me,timestep); 
       	fp = fopen(infile,"r"); 
       	fread(vtyp,size1,1,fp);  
       	fclose(fp); 
	
       	sprintf(infile,"%s.edot.%d.%s",prefix,me,timestep);
       	fp = fopen(infile,"r"); 
       	fread(edot,size1,1,fp);
       	fclose(fp); 
	
       	sprintf(infile,"%s.dstress.%d.%s",prefix,me,timestep); 
       	fp = fopen(infile,"r"); 
       	fread(stress,size3,1,fp); 
       	fclose(fp); 
	       	
       	sprintf(infile,"%s.topo_t.%d.%s",prefix,me,timestep);  
       	fp = fopen(infile,"r");  
       	fread(tpgt,size4,1,fp);  
       	fclose(fp);  
       	sprintf(infile,"%s.topo_b.%d.%s",prefix,me,timestep);  
       	fp = fopen(infile,"r");  
       	fread(tpgb,size4,1,fp);  
       	fclose(fp);  

       	for (k=1;k<=noyl;k++) {	 
	  for (i=1;i<=noxl;i++) {
	    for (j=1;j<=nozl;j++) {
               	nodel = j + (i-1)*nozl + (k-1)*nozl*noxl;
               	jg = locz[m]*(nozl-1) + j;
               	ig = locx[m]*(noxl-1) + i;
               	kg = locy[m]*(noyl-1) + k;
               	nodeg = jg + (ig-1)*noz + (kg-1)*noz*nox;
		//  fprintf(stderr,"m %d i %d j %d k %d nodel %d ig %d jg %d kg %d nodeg %d\n",m,i,j,k,nodel,ig,jg,kg,nodeg);
	
               	/* Assign values from this node to global mesh */
               	X[nodeg-1] = x[nodel]; /* -1 to start at zeroth element of output vector */
               	Y[nodeg-1] = y[nodel];
               	Z[nodeg-1] = z[nodel];
               	TEMP[nodeg-1] = temp[nodel];
		PRES[nodeg-1] = pres[nodel];
		VISC[nodeg-1] = visc[nodel];
		if (fcomp == 1)
		  COMP[nodeg-1] = comp[nodel];            	
               	EDOT[nodeg-1] = edot[nodel];
	        VTYP[nodeg-1] = vtyp[nodel];    		
		for (d=1;d<=nsd;d++){
		    VEL[(nodeg-1)*nsd + d -1] = vel[(nodel-1)*nsd + d];
		}
		for (d=1;d<=2*nsd;d++){ 
		    STRESS[(nodeg-1)*2*nsd + d -1] = stress[(nodel-1)*2*nsd + d];
		}
               	/* For top or bottom mesh get topography */         
		if (locz[m] == 0) { /* bottom*/
		  nodesl = i + (k-1)*noxl;               
		  nodesg = ig + (kg-1)*nox;
		  TPGB[nodesg-1] = tpgb[nodesl];
		}
		if (locz[m] == nprocz-1 ) {  /* top */
		  nodesl = i + (k-1)*noxl;               
		  nodesg = ig + (kg-1)*nox;
		  TPGT[nodesg-1] = tpgt[nodesl];
		}                	
	    } /* i */
	  } /* j */
	} /* k */
    } /* me */

    if (argc > 3) { /* Dimensionalize */
      for (i=0; i<nno; i++) {
	X[i] = X[i]*r2d;
	Y[i] = Y[i]*r2d;
	Z[i] = Z[i]*erad/1000;      /* km */
	TEMP[i] = To + TEMP[i]*dT;
	PRES[i] = pcon*PRES[i];
	VISC[i] = eta*VISC[i];
	EDOT[i] = econ*EDOT[i];
      }
      for (i=0; i<nno*nsd; i++) {
	VEL[i] = vcon*VEL[i];
      }
      for (i=0; i<nno*2*nsd; i++) {	
	STRESS[i] = pcon*STRESS[i];
      }
      for (k=1;k<=noy;k++) {	 /* Topography in km */ 
	for (i=1;i<=nox;i++) {
	  j = 1; /* bottom */            	
	  nodeg = j + (i-1)*noz + (k-1)*noz*nox;
	  nodesg = i + (k-1)*nox;
	  TPGB[nodesg] = ((pcon/g/1000)/(drhotop-COMP[nodeg]*drhochem))*TPGB[nodesg];
	  j = noz; /* top */            	
	  nodeg = j + (i-1)*noz + (k-1)*noz*nox;
	  nodesg = i + (k-1)*nox;
	  TPGT[nodesg] = ((pcon/g/1000)/(drhobot=COMP[nodeg]*drhochem))*TPGT[nodesg];
	}
      }
    }


    /* Write Full Mesh Data to Combined files */
    fprintf(stderr,"Writing Full Mesh Data to Files\n");
    size1g = (nno)*sizeof(float);
    size2g = (3*nno)*sizeof(float);
    size3g = (6*nno)*sizeof(float);
    size4g = (nsf)*sizeof(float);
    //  fprintf(stderr,"size1g %d size2g %d size3g %d size4g %d\n",
    //	   size1g,size2g,size3g,size4g);

    sprintf(outfile,"comb_%s.x.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(X,size1g,1,fp);   
    fclose(fp);
    sprintf(outfile,"comb_%s.y.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(Y,size1g,1,fp);   
    fclose(fp);
    sprintf(outfile,"comb_%s.z.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(Z,size1g,1,fp);   
    fclose(fp);
    sprintf(outfile,"comb_%s.temp.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(TEMP,size1g,1,fp);   
    fclose(fp);
    sprintf(outfile,"comb_%s.pres.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(PRES,size1g,1,fp);   
    fclose(fp);
    sprintf(outfile,"comb_%s.visc.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(VISC,size1g,1,fp);   
    fclose(fp);  
    sprintf(outfile,"comb_%s.vtyp.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(VTYP,size1g,1,fp);   
    fclose(fp);    
    sprintf(outfile,"comb_%s.edot.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(EDOT,size1g,1,fp);   
    fclose(fp);   
    sprintf(outfile,"comb_%s.vel.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(VEL,size2g,1,fp);   
    fclose(fp);    
    sprintf(outfile,"comb_%s.stress.%s",prefix,timestep);
    fp = fopen(outfile,"w");
    fwrite(STRESS,size3g,1,fp);   
    fclose(fp);    
    sprintf(outfile,"comb_%s.topot.%s",prefix,timestep); 
    fp = fopen(outfile,"w"); 
    fwrite(TPGT,size4g,1,fp);    
    fclose(fp); 
    sprintf(outfile,"comb_%s.topob.%s",prefix,timestep); 
    fp = fopen(outfile,"w"); 
    fwrite(TPGB,size4g,1,fp);    
    fclose(fp); 
    
}
  
