/********************************************************************************************************************
 *
 *   Program   : epsurf																			
 *   Purpose   : calculate MEP extrema and electronic descriptors on a user defined isodensity surface from two grids
 *   Input     : density.cube potential.cube densityisosurface (format g98)
 *   Output    : .log and .xyz - files with appended extrema, written to cwd
 *   Usage     : ./epsurf densityfile potentialfile surface
 *   Compiling : gcc -o epsurf epsurf.c -lm -O5
 *   Modified  : 03/Aug/04
 * 
 ********************************************************************************************************************/
 
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define STEP     2.00          
#define MAXLINE  200
#define MAXAT    100
#define MAXCOUNT 10000
#define B2A      0.5291772
#define H2KC     627.50956
#define H2KJ     2625.5     /************************************************************
#define PSLIM    0.000       * +ve surface limit: +0.035, zero considers entire surface *
#define NSLIM    0.000       * -ve surface limit: -0.035                                *
#define DEBUG    0           ************************************************************/

int main(int argc, char *argv[])
{
 char line[MAXLINE], title[MAXLINE];
 char sym[18][2]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar"};
   
 int i, j, k, m, umin, num_umin, umax, num_umax, flag;
 int count;
 int numat, numx, numy, numz, numden, numpot, numsurf;
 int at[MAXAT];
 int no_full, no_left, minkept, maxkept;
 int *surfx, *surfy, *surfz;

 float d, x0, y0, z0, xstep, ystep, zstep, junk, dist, SURF;
 float xcar[MAXAT], ycar[MAXAT], zcar[MAXAT];
 float newep, oldep;
 float *density, *potential, *surfden, *surfpot;
 float *x, *y, *z;
 int *locmin, *locmax, minu[MAXCOUNT], maxu[MAXCOUNT];
 
 /*************************
  * statistical variables *
  *************************/
 
 float avepot_pos, avepot_neg, dev, varpos, varneg, balance;
 float negpot, pospot, avepot, numneg, numpos, minep, maxep;
 
 FILE *fden, *fpot, *flog, *fxyz;
 
 div_t ans;
 time_t t1,t2;
  
 /************
  * get time *
  ************/
 
 (void) time(&t1);
 
 /***************************
  * parse command line args *
  ***************************/
 
 if (argc != 4)
   {
    printf("\nThis program locates MEP extrema and calculates Politzer-like statistical\n");
    printf("descriptors of the MEP on a user defined electron density isosurface from\n");
    printf("two grids (files generated with the g98 'cube' command). Three arguments should be specified.\n\n");
    printf("Usage: %s densityfile potentialfile isodensitysurface\n\n", argv[0]);
    return 1;
   }
 
 if ( (fden = fopen( argv[1], "r" )) == NULL )
   {
    printf("Error opening density file\n");
    return 2;
   }

 if ( (fpot = fopen( argv[2], "r" )) == NULL )
   {
    printf("Error opening potential file\n");
    return 2;
   }
   
  /*****************
   * open log file *
   *****************/
 
  flog = fopen(strcat(argv[1],".log"),"w");
   
  /********************************************
   * get user defined density in numeric form *
   ********************************************/
   
 SURF = strtod(argv[3],NULL);
 fgets(title, MAXLINE, fden);
 fgets(line, MAXLINE, fden);

 fscanf(fden, "%d %f %f %f", &numat, &x0, &y0, &z0);

 /***************************************
  * read grid info - assume cuboid grid *
  ***************************************/

 fscanf(fden, "%d %f %f %f", &numx, &xstep, &junk, &junk);
 fscanf(fden, "%d %f %f %f", &numy, &junk, &ystep, &junk);
 fscanf(fden, "%d %f %f %f", &numz, &junk, &junk, &zstep);

 /******************
  * read atom info *
  ******************/

 for(i=0; i<numat; i++)
 {
  fscanf(fden, "%d %f %f %f %f", &at[i], &junk, &xcar[i], &ycar[i], &zcar[i]);
 }
 
 /*********************************************************
  * get no full and not full lines in cube (six per line) *
  *********************************************************/

 ans = div(numz,6);
 no_full = ans.quot;
 no_left = ans.rem;
 
 /************************************
  * open an array for density values *
  ************************************/

 density = (float *) malloc(numx*numy*numz*sizeof(float));
 if (!density)
 {
  fprintf(flog,"Can't allocate memory for density -exiting\n");
  return 4;
 }

 /******************************************************
  * read in density values in numx*numy blocks of numz *
  ******************************************************/

 for(k=0; k<numx*numy*numz; k+= numz)
 {
  for(i=0; i<no_full*6; i+=6)
  {
   	 fscanf(fden,"%f%f%f%f%f%f",&density[k+i],&density[k+i+1],&density[k+i+2],&density[k+i+3],&density[k+i+4],&density[k+i+5]);
  }
  	for(j=0; j<no_left; j++)
  	{
    	fscanf(fden,"%f", &density[k + i + j]);
  	}
 }  
 numden = k+i+j-numz;
 fprintf(flog,"\ntotal number of density points scanned   = %d\n", numden);

 /**********************************
  * Done with density file - close *
  **********************************/

 fclose(fden);
 
 /******************************************************************
  * do same for potential cube - assume parameters same as density *
  ******************************************************************/

 potential = (float *) malloc(numx*numy*numz*sizeof(float));
 if (!potential)
 {
  fprintf(flog,"Can't allocate memory for potential -exiting\n");
  return 5;
 }

 for(i=0; i<numat+6; i++)
  fgets(line, MAXLINE, fpot);

 for(k=0; k<numx*numy*numz; k+= numz)
 {
  for(i=0; i<no_full*6; i+=6)
  {
    fscanf(fpot,"%f%f%f%f%f%f",&potential[k+i],&potential[k+i+1],&potential[k+i+2],&potential[k+i+3],&potential[k+i+4],&potential[k+i+5]);
  }
  	for(j=0; j<no_left; j++)
  	{
    fscanf(fpot,"%f", &potential[k + i + j]);
  	}
 } 
 numpot = k+i+j-numz;
 fprintf(flog,"total number of potential points scanned = %d\n\n", numpot);
 fclose(fpot);

 /**********************************************************
  * Now run through to count points with density within 1% *
  * of defined surface and set up arrays accordingly       *
  **********************************************************/

 numsurf = 0;
 for(i=0; i<numden; i++)
 {
   if (((density[i] - SURF) < 0.01*SURF) && ((density[i] - SURF) > -0.01*SURF))
   {
     numsurf++;
   }
 }

 surfden = (float *) malloc(numsurf*sizeof(float));
 surfpot = (float *) malloc(numsurf*sizeof(float));
 locmin  = (int *) malloc(numsurf*sizeof(float));
 locmax  = (int *) malloc(numsurf*sizeof(float));

 surfx = (int *) malloc(numsurf*sizeof(float));
 surfy = (int *) malloc(numsurf*sizeof(float));
 surfz = (int *) malloc(numsurf*sizeof(float));
 x = (float *) malloc(numsurf*sizeof(float));
 y = (float *) malloc(numsurf*sizeof(float));
 z = (float *) malloc(numsurf*sizeof(float));

 /***********************************************
  * Now run through to find points with density * 
  * close to SURF au and save values            *
  ***********************************************/
  
 numsurf = 0;
 for(i=0; i<numden; i++)
 {
   if (((density[i] - SURF) < 0.01*SURF) && ((density[i] - SURF) > -0.01*SURF))
   {
     surfden[numsurf] = density[i];
     surfpot[numsurf] = potential[i];

     surfx[numsurf] = (int) (i/(numy*numz));
     surfy[numsurf] = (int) (i/numz - surfx[numsurf]*numy);
     surfz[numsurf] = (int) (i - surfy[numsurf]*numz - surfx[numsurf]*numy*numz);

     numsurf++;
   }
 }
 
 /*************************************************************
  * ?should report that we are scanning +/- 1% au of surface? * 
  *************************************************************/
 
 fprintf(flog,"molecular surface   = %s\n",argv[3]);  
 fprintf(flog,"abs ESP cutoff      = %4.3f\n",PSLIM);
 fprintf(flog,"num ESPs on surface = %d\n\n",numsurf);

 /*******************************************
  * no need to always print surface details *
  *******************************************/

 #if DEBUG
 
 fprintf(flog,"     x          y          z         rho           V\n");
 
 for(i=0; i<numsurf; i++)
 {
   x[i] = x0 + surfx[i]*xstep;
   y[i] = y0 + surfy[i]*ystep;
   z[i] = z0 + surfz[i]*zstep;

   fprintf(flog,"%10.6f %10.6f %10.6f  %10.6f %10.6f\n", 
                           x[i], y[i], z[i], surfden[i], surfpot[i]); 
 }

 for(i=0; i<numsurf; i++)
   {
     for(j=i+1;j<numsurf; j++)
     {
       dist = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
     }
   } 
   
 #endif

for(i=0; i<numsurf; i++)
 {
   x[i] = x0 + surfx[i]*xstep;
   y[i] = y0 + surfy[i]*ystep;
   z[i] = z0 + surfz[i]*zstep;
 }

 for(i=0; i<numsurf; i++)
   {
     for(j=i+1;j<numsurf; j++)
     {
       dist = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
     }
   }

 free(density);
 free(potential);

 /*************************
  * get maximum potential *
  *************************/

for (i=0; i<numsurf; i++)
{
	if(surfpot[i] > maxep)
	{
	 maxep = surfpot[i];
	}
}

 /*************************
  * get minimum potential *
  *************************/

for (i=0; i<numsurf; i++)
{
	if(surfpot[i] < minep)
	{
	 minep = surfpot[i];
	}
}

 /****************************
  * Now to find local minima *
  ****************************/

 for (i=0; i<numsurf; i++)
 {
   m = i;
   count = 0;
   do
   {
     for (j=0; j<=numsurf; j++)
     {
      d=sqrt((x[m]-x[j])*(x[m]-x[j])+(y[m]-y[j])*(y[m]-y[j])+(z[m]-z[j])*(z[m]-z[j]));
      if( (d < STEP) && (surfpot[j] < surfpot[m]))
      {
        oldep = surfpot[m];
        newep = surfpot[j];
        m = j;
      }
     }
     count++;
     if(count > MAXCOUNT)
       break;

   } while (fabs((oldep - newep)) > 0.000000000001);
   locmin[i] = m;
 }

 /************************************
  * Identify unique minima from list *
  ************************************/

  umin = 0;

  for(i=0; i<numsurf; i++)
  {
    flag=0;
    for(j=0; j<umin; j++)
    {
      if (locmin[i] == minu[j])
       flag=1;
    }
    if (flag == 0)
    {
      minu[umin] = locmin[i];
      umin++;
    }
  }

  minkept=0;
  fprintf(flog,"minima\n",NSLIM);
  for (i=0; i<umin; i++)
  {
  if (surfpot[minu[i]] < NSLIM)
   { 													   
     minkept++;
     fprintf(flog,"H%-4d %12.6f %12.6f %12.6f\t\tEP=%10.4f\n", numat+minkept, x[minu[i]]*B2A, y[minu[i]]*B2A, z[minu[i]]*B2A, surfpot[minu[i]]);
   } 
  }

 /****************************
  * Now to find local maxima *
  ****************************/

 for (i=0; i<numsurf; i++)
 {
   m = i;
   count = 0;
   do
   {
     for (j=0; j<=numsurf; j++)
     {
      d=sqrt((x[m]-x[j])*(x[m]-x[j])+(y[m]-y[j])*(y[m]-y[j])+(z[m]-z[j])*(z[m]-z[j]));
      if( (d < STEP) && (surfpot[j] > surfpot[m]))
      {
        oldep = surfpot[m];
        newep = surfpot[j];
        m = j;
      }
     }
     count++;
     if(count > MAXCOUNT)
       break;

   } while (fabs((oldep - newep)) > 0.000000000001);
   locmax[i] = m;
 }

 /************************************
  * Identify unique maxima from list *
  ************************************/

 umax = 0;

 for(i=0; i<numsurf; i++)
 {
   flag=0;
   for(j=0; j<umax; j++)
   {
     if (locmax[i] == maxu[j])
     flag=1;
   }
   if (flag == 0)
   {
     maxu[umax] = locmax[i];
     umax++;
   }
 }

 maxkept=0;
 fprintf(flog,"maxima\n",PSLIM);
 for (i=0; i<umax; i++)
 {
   if (surfpot[maxu[i]] > PSLIM)
   {												
     maxkept++;							                    
     fprintf(flog,"F%-3d %12.6f %12.6f %12.6f\t\tEP=%10.4f\n", numat+minkept+maxkept, x[maxu[i]]*B2A, y[maxu[i]]*B2A, z[maxu[i]]*B2A, surfpot[maxu[i]]);
   }
 }

  /************************************
  * get num umin we are interested in *
  ************************************/
   
   num_umin = 0;
   for (i=0; i<umin; i++)
   {
     if (surfpot[minu[i]] < NSLIM) 
     { 
       num_umin++;
     } 
   }
  
  /************************************
  * get num umax we are interested in *
  ************************************/
  
   num_umax = 0;
   for (i=0; i<umax; i++)
   {
     if (surfpot[maxu[i]] > PSLIM) 
     { 
       num_umax++;
     } 
   }
 
 /*****************
  * print headers *
  *****************/
 
 fprintf(flog,"\nsurface ESP properties\n",argv[3]);
 fprintf(flog,"     	                     au        kcal/mol      kJ/mol\n"); 

 /*********************************
  * report global maxep and minep *
  *********************************/

 fprintf(flog,"max ESP          =  %12.4f %12.4f %12.4f\n",maxep,maxep*H2KC,maxep*H2KJ);
 fprintf(flog,"min ESP          =  %12.4f %12.4f %12.4f\n",minep,minep*H2KC,minep*H2KJ);

 /***********************************************
  * get average -ve surface potential < cuttoff *
  ***********************************************/
  
 negpot = 0;
 for (i=0; i<numsurf; i++)
  { 
   if (surfpot[i] < NSLIM) 
   { 
     negpot += surfpot[i];
	 numneg++;
   } 
  }
 avepot_neg  = negpot / numneg;
 fprintf(flog,"-ve ESP average  =  %12.4f %12.4f %12.4f\n",avepot_neg,avepot_neg*H2KC,avepot_neg*H2KJ);
      
 /***********************************************
  * get average +ve surface potential < cuttoff *
  ***********************************************/

 pospot = 0;
 for (i=0; i<numsurf; i++)
  { 
   if (surfpot[i] > PSLIM)  
   { 
     pospot += surfpot[i];
	 numpos++;
   }
  }   
 avepot_pos = pospot / numpos;
 fprintf(flog,"+ve ESP average  =  %12.4f %12.4f %12.4f\n",avepot_pos,avepot_pos*H2KC,avepot_pos*H2KJ);

 /*************************************************************
  * get average MEP deviation over entire surface of interest *
  *************************************************************/


 avepot=0;
 for (i=0; i<numsurf; i++)
 {
    avepot += surfpot[i];
 }
 avepot = avepot/numsurf;

 dev = 0;
 for (i=0; i<numsurf; i++)
 {
    if(fabs(surfpot[i]) > PSLIM)
	{
		dev += fabs(surfpot[i] - avepot);
	}
 }															
 fprintf(flog,"    ESP ave dev  =  %12.4f %12.4f %12.4f\n\n",dev/numsurf,(dev/numsurf)*H2KC,(dev/numsurf)*H2KJ);
 
 /****************************************
  * print header for values with units^2 *
  ****************************************/
 
 fprintf(flog,"     	                     au^2      kcal/mol^2    kJ/mol^2\n"); 
 
 /*****************
  * get variances *
  *****************/
 
 varneg = 0;
 for (i=0; i<numsurf; i++)
  { 
   if (surfpot[i] < NSLIM) 
   { 
   	 varneg += ((surfpot[i] - avepot_neg) * (surfpot[i] - avepot_neg));
   } 
  }
 fprintf(flog,"-ve ESP variance =  %12.4f %12.4f %12.4f\n",varneg,varneg*H2KC,varneg*H2KJ);

 varpos = 0;
 for (i=0; i<numsurf; i++)
  { 
   if (surfpot[i] > PSLIM) 
   { 
   	 varpos += ((surfpot[i] - avepot_pos) * (surfpot[i] - avepot_pos));
   } 
  }
 fprintf(flog,"+ve ESP variance =  %12.4f %12.4f %12.4f\n",varpos,varpos*H2KC,varpos*H2KJ);
 fprintf(flog,"tot ESP variance =  %12.4f %12.4f %12.4f\n\n",varpos+varneg,(varpos+varneg)*H2KC,(varpos+varneg)*H2KJ);

 /******************
  * report balance *
  ******************/

 balance = (varpos*varneg)/((varpos+varneg)*(varpos+varneg));
 fprintf(flog,"    ESP balance  =  %12.4f (dimensionless)\n\n",balance);

 /*******************
  * debugging stuff *
  *******************/
  
 #if DEBUG
 fprintf(flog,"    numat        =  %d \n",numat);
 fprintf(flog,"    num_umin     =  %d \n",num_umin);
 fprintf(flog,"    num_umax     =  %d \n\n",num_umax);
 fprintf(flog,"    umin         =  %d \n",umin);
 fprintf(flog,"    umax         =  %d \n\n",umax);
 #endif
 
 /***********************************
  * get time, report total job time *
  ***********************************/
 
 (void) time(&t2);
 fprintf(flog,"    job time     =  %4.2d minutes, %2.0f seconds\n\n",((int)t2-t1)/60,fmod((int)t2-t1,60));
    
 /*************************************
  * close log file, open xyz outfile  *
  *************************************/
   
 fclose(flog);
 fxyz = fopen(strcat(argv[2],".xyz"),"w");
 
 /******************************
  * write headers xyz outfile, *
  * line1: numAtoms, required  *
  * line2: title, optional     *
  ******************************/
 
 fprintf(fxyz, "%d\n", numat+num_umin+num_umax);
 fprintf(fxyz, "%s", title);
	
 /******************************************************************************
  * write molecule specification in xyz format, convert to angstroms using B2A *
  ******************************************************************************/
  
 for (i=0; i<numat; i++)
 {
   fprintf(fxyz, "%4s %12.6f %12.6f %12.6f\n", sym[at[i]-1], xcar[i]*B2A,ycar[i]*B2A,zcar[i]*B2A);
 }

 /************************************
  * append atoms denoting minima (H) *
  ************************************/
	
  for (i=0; i<umin; i++)
  {
   if (surfpot[minu[i]] < NSLIM)
   { 
     fprintf(fxyz,"   H %12.6f %12.6f %12.6f\n", x[minu[i]]*B2A, y[minu[i]]*B2A, z[minu[i]]*B2A);
   } 
  }
	
 /************************************
  * append atoms denoting maxima (F) *
  ************************************/

  for (i=0; i<umax; i++)  
  { 
 	if (surfpot[maxu[i]] > PSLIM) 
   { 
     fprintf(fxyz,"   F %12.6f %12.6f %12.6f\n", x[maxu[i]]*B2A, y[maxu[i]]*B2A, z[maxu[i]]*B2A);
   } 
  }
	
 /******************************************
  * close xyz outfile, end main() function *
  ******************************************/
		
  fclose(fxyz);

} 

/* TODO ****************************************************************************************************
 * calculate distance and angles between MEP extrema, report geometric relation of local minima and maxima *                                 *
 ***********************************************************************************************************/

/* BUGS **************************************
 *                                           *
 *                                           *
 *********************************************/
