///////////////////////////////////
// DATA STRUCTURES //
///////////////////////////////////

/////////////////////////////////
// ADITIONAL CODE //
////////////////////////////////

float *U;

#include<../useful_routines/useful_routines.c>

////////////////////////////////////
// GLOBAL VARIABLES //
///////////////////////////////////
typedef struct
{
  char nameGalaxy[500];
  int nParticles;
  double pos[3];
  double vel[3];
  int flagRotation;
  double inclinationAngle;
  double positionAngle;
  int types[6];
}galaxies;

galaxies *gal;
particulas *mergerParticles;  

////////////////////////////////////
// ROUTINES //
///////////////////////////////////
int printAscii(int partTotal[]);

int main(int argc, char *argv[])
{
  
  int i, j, k;
  int nGalaxies, nTotal, type, centerIn;
  int totalParticlesPerType[6], stored[6];
  char *infile, outfile[500],headerFile[500];
  FILE *fParam, *fOutParam, *fLog;
  CM *massCenters;

  fLog = fopen("configure_merger.log","w");

  infile = argv[1];
  
  // Reading parameters file
  fParam = fopen(infile,"r");
  
  returnRead = fscanf(fParam,"%d",&nGalaxies);
  returnRead = fscanf(fParam,"%s",outfile);
  returnRead = fscanf(fParam,"%s",headerFile);
  returnRead = fscanf(fParam,"%d",&centerIn);
  
  int indexMinPerType[nGalaxies][6], indexMaxPerType[nGalaxies][6];

  gal = (galaxies *)malloc((size_t)nGalaxies*sizeof(galaxies));
  if(gal == NULL){
    printf("Allocation of gal failed\n");
    exit(0);
  }

  massCenters = (CM *)malloc((size_t)(nGalaxies+1)*sizeof(CM));
  if(massCenters == NULL){
    printf("Allocation of massCenters failed\n");
    exit(0);
  }

  nTotal = 0;
  for( type=0; type<6; type++)
    totalParticlesPerType[type] = 0;
  
  for( i=0; i<nGalaxies; i++ )
    {
      returnRead = fscanf(fParam,"%s",gal[i].nameGalaxy);
      returnRead = fscanf(fParam,"%lf %lf %lf",&gal[i].pos[X],&gal[i].pos[Y],&gal[i].pos[Z]);
      returnRead = fscanf(fParam,"%lf %lf %lf",&gal[i].vel[X],&gal[i].vel[Y],&gal[i].vel[Z]);
      returnRead = fscanf(fParam,"%d",&gal[i].flagRotation);
      returnRead = fscanf(fParam,"%lf %lf",&gal[i].inclinationAngle,&gal[i].positionAngle);
      returnRead = fscanf(fParam,"%d %d %d %d %d %d",
	     &gal[i].types[0],&gal[i].types[1],&gal[i].types[2],&gal[i].types[3],&gal[i].types[4],&gal[i].types[5]);

      // Reading Headers
      read_header_gadget(gal[i].nameGalaxy);

      for( type=0; type<6; type++)
	{

	  for( j=0; j<6; j++)
	    if( gal[i].types[j]==type )
	      totalParticlesPerType[type] = totalParticlesPerType[type] + Header.npartTotal[j];

	  indexMinPerType[i][type] = indexMaxPerType[i][type] = 0;
	   if( Header.npartTotal[type]>0 )
	     {
	       for( j=0; j<type; j++)
		 indexMinPerType[i][type] = indexMinPerType[i][type] + Header.npartTotal[j];
	       indexMaxPerType[i][type] = indexMinPerType[i][type] + Header.npartTotal[type];
	       nTotal = nTotal + Header.npartTotal[type];
	     }

	}

    }
  
  fprintf(fLog,"\nConfigure merger for : %s\n\n",outfile);
  for( type=0; type<6; type++)
    fprintf(fLog,"Total particles per type %d = %d\n",type,totalParticlesPerType[type]);
  fprintf(fLog,"\nTotal particles in initial conditions = %d\n",nTotal);

  // Printing used parameters
  fOutParam = fopen("used_parameters.output","w");
  
  fprintf(fOutParam,"Number of galaxies = %d\n",nGalaxies);
  fprintf(fOutParam,"Output file = %s\n",outfile);
  fprintf(fOutParam,"Final Header from = %s\n\n",headerFile);
  for( i=0; i<nGalaxies; i++ )
    {
      fprintf(fOutParam,"Galaxy %d = %s\n",i,gal[i].nameGalaxy);
      fprintf(fOutParam,"Position  = ( %lf, %lf, %lf )\n",gal[i].pos[X],gal[i].pos[Y],gal[i].pos[Z]);
      fprintf(fOutParam,"Velocity = ( %lf,  %lf, %lf )\n",gal[i].vel[X],gal[i].vel[Y],gal[i].vel[Z]);
      fprintf(fOutParam,"Rotation = %d\n",gal[i].flagRotation);
      fprintf(fOutParam,"Inclination angle = %lf\n",gal[i].inclinationAngle);
      fprintf(fOutParam,"Position angle = %lf\n",gal[i].positionAngle);
      fprintf(fOutParam,"Components in = %d %d %d %d %d %d\n\n",
	      gal[i].types[0],gal[i].types[1],gal[i].types[2],gal[i].types[3],gal[i].types[4],gal[i].types[5]);
    }
  fclose(fParam);
  fclose(fOutParam);

  for( i=0; i<nGalaxies; i++ )
    {    
      gal[i].inclinationAngle = gal[i].inclinationAngle*M_PI/180.0; 
      gal[i].positionAngle = gal[i].positionAngle*M_PI/180.0; 
    }

 mergerParticles = (particulas *)malloc((size_t)nTotal*sizeof(particulas));
 if(mergerParticles == NULL){
   printf("Allocation of mergerParticles failed\n");
   exit(0);
 }  
 
 U = (float *)malloc((size_t)totalParticlesPerType[0]*sizeof(float));
 if(U == NULL){
   printf("Allocation of U failed\n");
   exit(0);
 }

 stored[0] = 0;
 stored[1] = totalParticlesPerType[0];
 stored[2] = stored[1] + totalParticlesPerType[1];
 stored[3] = stored[2] + totalParticlesPerType[2];
 stored[4] = stored[3] + totalParticlesPerType[3];
 stored[5] = stored[4] + totalParticlesPerType[4];

  for( type=0; type<6; type++)
    fprintf(fLog,"Initial index per type %d = %d\n",type,stored[type]);

  massCenters[nGalaxies].cm[X] = 0.0;
  massCenters[nGalaxies].cm[Y] = 0.0;
  massCenters[nGalaxies].cm[Z] = 0.0;
  
  massCenters[nGalaxies].vcm[X] = 0.0;
  massCenters[nGalaxies].vcm[Y] = 0.0;
  massCenters[nGalaxies].vcm[Z] = 0.0;
  
  massCenters[nGalaxies].M = 0.0;

  for( i=0; i<nGalaxies; i++ )
    {

      printf("\nConfiguring galaxy: %s\n",gal[i].nameGalaxy);

      // Reading galaxy
      read_gadget1(gal[i].nameGalaxy);

      gal[i].nParticles = N_part_total;

      // Checking center of mass
      centerMass(NULL, gal[i].nParticles, 0, gal[i].nParticles, &cmdummy);
      printf("RCM : %e %e %e\n",cmdummy.cm[X],cmdummy.cm[Y],cmdummy.cm[Z]);
      printf("VCM : %e %e %e\n",cmdummy.vcm[X],cmdummy.vcm[Y],cmdummy.vcm[Z]);

      // Moving galaxy to zero
      totalTranslationMinus(&cmdummy, gal[i].nParticles);

      // Checking center of mass
      centerMass(NULL, gal[i].nParticles, 0, gal[i].nParticles, &cmdummy);
      printf("RCM : %e %e %e\n",cmdummy.cm[X],cmdummy.cm[Y],cmdummy.cm[Z]);
      printf("VCM : %e %e %e\n",cmdummy.vcm[X],cmdummy.vcm[Y],cmdummy.vcm[Z]);

      // Rotanting galaxies
      if( gal[i].flagRotation==1  )
	rotationInclinationPosition(gal[i].inclinationAngle, gal[i].positionAngle, gal[i].nParticles);

      // New position of center of mass
      cmdummy.cm[X] = gal[i].pos[X];
      cmdummy.cm[Y] = gal[i].pos[Y];
      cmdummy.cm[Z] = gal[i].pos[Z];

      // New velocity of center of mass
      cmdummy.vcm[X] = gal[i].vel[X];
      cmdummy.vcm[Y] = gal[i].vel[Y];
      cmdummy.vcm[Z] = gal[i].vel[Z];

      // Moving galaxy to new center of mass
      totalTranslationPlus(&cmdummy, gal[i].nParticles);

      // Checking the new center of mass
      centerMass(NULL, gal[i].nParticles, 0, gal[i].nParticles, &cmdummy);
      printf("NEW RCM : %e %e %e\n",cmdummy.cm[X],cmdummy.cm[Y],cmdummy.cm[Z]);
      printf("NEW VCM : %e %e %e\n",cmdummy.vcm[X],cmdummy.vcm[Y],cmdummy.vcm[Z]);

      massCenters[i].cm[X] = cmdummy.cm[X];
      massCenters[i].cm[Y] = cmdummy.cm[Y];
      massCenters[i].cm[Z] = cmdummy.cm[Z];

      massCenters[i].vcm[X] = cmdummy.vcm[X];
      massCenters[i].vcm[Y] = cmdummy.vcm[Y];
      massCenters[i].vcm[Z] = cmdummy.vcm[Z];

      massCenters[i].M = cmdummy.M;

      massCenters[nGalaxies].cm[X] = massCenters[nGalaxies].cm[X] + cmdummy.cm[X]*cmdummy.M; 
      massCenters[nGalaxies].cm[Y] = massCenters[nGalaxies].cm[Y] + cmdummy.cm[Y]*cmdummy.M; 
      massCenters[nGalaxies].cm[Z] = massCenters[nGalaxies].cm[Z] + cmdummy.cm[Z]*cmdummy.M; 

      massCenters[nGalaxies].vcm[X] = massCenters[nGalaxies].vcm[X] + cmdummy.vcm[X]*cmdummy.M; 
      massCenters[nGalaxies].vcm[Y] = massCenters[nGalaxies].vcm[Y] + cmdummy.vcm[Y]*cmdummy.M; 
      massCenters[nGalaxies].vcm[Z] = massCenters[nGalaxies].vcm[Z] + cmdummy.vcm[Z]*cmdummy.M; 

      massCenters[nGalaxies].M = massCenters[nGalaxies].M + cmdummy.M;

      // New types of particles
      for( type=0; type<6; type++ )
	for( j=0; j<gal[i].nParticles; j++ )
	  if( particles[j].type==type )
	    particles[j].type = gal[i].types[type];
      
      // Storing particles
      for( type=0; type<6; type++ )
	{
	  fprintf(fLog,"gal = %d - type = %d - indexmin = %d - indexmax = %d\n",i,type,indexMinPerType[i][type],indexMaxPerType[i][type]);
	  for( j=0; j<6; j++ )
	    {
	      if( gal[i].types[j]==type )
		{
		  for( k=indexMinPerType[i][type]; k<indexMaxPerType[i][type]; k++ )
		    {
		      mergerParticles[stored[type]] = particles[k];
		      mergerParticles[stored[type]].id = stored[type];  
		      
		      if( type==0 )
			U[stored[type]] = gaspro[k].U;
		      
		      stored[type]++;  
		    }
		}
	    }
	}

      free(particles);    
      free(gaspro);
    }

  fprintf(fLog,"\n");
  for( type=0; type<6; type++ )
    fprintf(fLog,"Final index per type %d = %d\n",type,stored[type]);
  
  particles = (particulas *)malloc((size_t)nTotal*sizeof(particulas));
  if(particles == NULL){
    printf("Allocation of particles failed\n");
    exit(0);
  }
  
  massCenters[nGalaxies].cm[X] = massCenters[nGalaxies].cm[X]/massCenters[nGalaxies].M;
  massCenters[nGalaxies].cm[Y] = massCenters[nGalaxies].cm[Y]/massCenters[nGalaxies].M; 
  massCenters[nGalaxies].cm[Z] = massCenters[nGalaxies].cm[Z]/massCenters[nGalaxies].M; 
  
  massCenters[nGalaxies].vcm[X] = massCenters[nGalaxies].vcm[X]/massCenters[nGalaxies].M;
  massCenters[nGalaxies].vcm[Y] = massCenters[nGalaxies].vcm[Y]/massCenters[nGalaxies].M;
  massCenters[nGalaxies].vcm[Z] = massCenters[nGalaxies].vcm[Z]/massCenters[nGalaxies].M;

  printf("Final centros:\n");
  for( i=0; i<(nGalaxies+1); i++ )
    {
      printf("\n %d %e %e %e\n",i,massCenters[i].cm[X],massCenters[i].cm[Y],massCenters[i].cm[Z]);
      printf(" %d %e %e %e\n",i,massCenters[i].vcm[X],massCenters[i].vcm[Y],massCenters[i].vcm[Z]);
    }    

  for( i=0; i<nTotal; i++ )
    {
      
      particles[i] = mergerParticles[i];
      
      particles[i].pos[X] = particles[i].pos[X] - massCenters[centerIn].cm[X];
      particles[i].pos[Y] = particles[i].pos[Y] - massCenters[centerIn].cm[Y];
      particles[i].pos[Z] = particles[i].pos[Z] - massCenters[centerIn].cm[Z];
      
      particles[i].vel[X] = particles[i].vel[X] - massCenters[centerIn].vcm[X];
      particles[i].vel[Y] = particles[i].vel[Y] - massCenters[centerIn].vcm[Y];
      particles[i].vel[Z] = particles[i].vel[Z] - massCenters[centerIn].vcm[Z];
      
    }

  printAscii( totalParticlesPerType );  

  write_gadget1(headerFile, outfile);

  free(gal);
  free(mergerParticles);
  free(U);
  free(particles);
  free(massCenters);
  fclose(fLog);    

  return 0;    
}


int printAscii(int partTotal[])
{
  int i, type;
  char outfiles[500];
  FILE *fAscii;

  N_min = 0;
  
  for( type=0; type<6; type++)
    {

      N_max = N_min + partTotal[type];	 

      if( partTotal[type] > 0 )
	{
	  sprintf(outfiles,"configured.%d",type);
	  //  FILE *fAscii;
	  fAscii = fopen(outfiles,"w");
	  if(fAscii==NULL) printf("No se pudo abrir %s\n",outfiles);  
	  
	  for(i=N_min;i<N_max;i++)
	    {
	      
#ifdef LONGIDS
	      fprintf(fAscii,"%lu",particles[i].id);
#else
	      fprintf(fAscii,"%u",particles[i].id);
#endif
	      
	      fprintf(fAscii," %f %f %f %f %f %f %f",
		      particles[i].pos[X],particles[i].pos[Y],particles[i].pos[Z],
		      particles[i].vel[X],particles[i].vel[Y],particles[i].vel[Z],
		      particles[i].mass);
	      /*
#ifdef LONGIDS
	      fprintf(fAscii,"%lu",mergerParticles[i].id);
#else
	      fprintf(fAscii,"%u",mergerParticles[i].id);
#endif
	      
	      fprintf(fAscii," %f %f %f %f %f %f %f",
		      mergerParticles[i].pos[X],mergerParticles[i].pos[Y],mergerParticles[i].pos[Z],
		      mergerParticles[i].vel[X],mergerParticles[i].vel[Y],mergerParticles[i].vel[Z],
		      mergerParticles[i].mass);
*/
	      
	      if( type == 0 &&   Header.npartTotal[0]>0 )
	      fprintf(fAscii," %f",gaspro[i].U);

	      
	      fprintf(fAscii,"\n");	
	      
	    }
	  fclose(fAscii);
	}

      N_min = N_max;

    }

  return 0;
}
