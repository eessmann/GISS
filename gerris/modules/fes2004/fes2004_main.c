// standard include //
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>   
#include <strings.h>
#include <malloc.h>

#include "fes2004_lib.h"



/*####################################################*/
/*                                                    */
/*    read an simple  file construct with             */
/*    lat lon time structs                            */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/
void load_input_file(FILE *in, int nb, double *lat ,double *lon, double *time)
{
  int i,val;
  char carac[256];

  for(i=0;i<nb;i++) 
  {
    val=fscanf(in, "%lf %lf %lf",&lat[i],&lon[i],&time[i]);
    if (val!=3) {sprintf(carac,"fatal error during the read of the input file at line : %d",i+1); print_error_4(carac);}
    if(lon[i]>180) lon[i]-=360;
  }

}

/*####################################################*/
/*                                                    */
/*    search the number of lines in filename          */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/
int  compter_lignes(char *filename)
{

  FILE *fichier;
  int lignes = 0;
  char line[1024];

  fichier=fopen(filename,"r");
  if(fichier==NULL) print_error_5(filename);
  while (!feof(fichier))
    {
      fgets(line,sizeof(line),fichier);
      if(!feof(fichier))lignes++;
    }
  printf("%s:  %d lignes\n", filename, lignes);
  fclose(fichier);
  return(lignes);

} 


/*####################################################*/
/*                                                    */
/*     vectors variables allocation                    */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

void alloc_vectors(int nb, double **lat,double **lon,double **time,double **prediction,double ***amplitude,double ***phase)
{
  int i;
  *lat=calloc(nb,sizeof(double));
  *lon=calloc(nb,sizeof(double));
  *time=calloc(nb,sizeof(double));
  *prediction=calloc(nb,sizeof(double));
  *amplitude=malloc(nb*sizeof(double *));
  *phase=malloc(nb*sizeof(double *));
  
  for(i=0;i<nb;i++)
    {
      (*amplitude)[i]=calloc(14,sizeof(double));
      (*phase)[i]=calloc(14,sizeof(double));
    }
  
}

/*####################################################*/
/*                                                    */
/*     generate a row output for the prediction       */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 07Oct2005              */
/*                                                    */
/*####################################################*/
int  prediction_output(int nb_position,double *prediction,double *lat,double *lon,double *time,int time_reference,char *out_filename)
{
  int cnt;
  FILE *out;
  char ref[56];
  out=fopen(out_filename,"w");
  if(time_reference==0) sprintf(ref,"01/01/1950");
  if(time_reference==1) sprintf(ref,"01/01/1958");
  if(time_reference==2) sprintf(ref,"01/01/1985");
  if(time_reference==3) sprintf(ref,"01/01/2000");
  fprintf(out,"# FES2004 prediction in meter\n# read : lat, lon, time, prediction\n#the time is the number of hour since the %s\n###############################################\n",ref);
  for(cnt=0;cnt<nb_position;cnt++) fprintf(out,"%12.6lf %12.6lf %12.6lf %12.6lf \n",lat[cnt],lon[cnt],time[cnt],prediction[cnt]);
  fclose(out);
}


/*####################################################*/
/*                                                    */
/*     generate a row output for the extraction       */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 07Oct2005              */
/*                                                    */
/*####################################################*/
int  extraction_output(int nb_position,double **amplitude,double **phase,double *lat,double *lon,char *out_filename)
{
  int cnt;
  FILE *out;
  out=fopen(out_filename,"w");
  fprintf(out,"###############################################\n#\n# FES2004 extraction -- amplitude in meter and phase lag in degree\n# read header : Pt_X -- lat lon\n#Wave_name  Amp  phi \n#\n###############################################\n");
  for(cnt=0;cnt<nb_position;cnt++)
  { 
    fprintf(out,"----------------------------------------------------------\n");
    fprintf(out,"Pt %6d     %12.6lf N %12.6lf E\n",cnt,lat[cnt],lon[cnt]);
    if(amplitude[cnt][2]==-9999.00){ fprintf(out,"Masked point\n");continue;}
    fprintf(out,"M2    %12.6lf %12.6lf\n",amplitude[cnt][3],phase[cnt][3]);
    fprintf(out,"S2    %12.6lf %12.6lf\n",amplitude[cnt][13],phase[cnt][13]);
    fprintf(out,"K1    %12.6lf %12.6lf\n",amplitude[cnt][1],phase[cnt][1]);
    fprintf(out,"K2    %12.6lf %12.6lf\n",amplitude[cnt][2],phase[cnt][2]);
    fprintf(out,"N2    %12.6lf %12.6lf\n",amplitude[cnt][9],phase[cnt][9]);
    fprintf(out,"2N2   %12.6lf %12.6lf\n",amplitude[cnt][0],phase[cnt][0]);
    fprintf(out,"O1    %12.6lf %12.6lf\n",amplitude[cnt][10],phase[cnt][10]);
    fprintf(out,"P1    %12.6lf %12.6lf\n",amplitude[cnt][11],phase[cnt][11]);
    fprintf(out,"Q1    %12.6lf %12.6lf\n",amplitude[cnt][12],phase[cnt][12]);
    fprintf(out,"M4    %12.6lf %12.6lf\n",amplitude[cnt][4],phase[cnt][4]);
    fprintf(out,"Mf    %12.6lf %12.6lf\n",amplitude[cnt][5],phase[cnt][5]);
    fprintf(out,"Mm    %12.6lf %12.6lf\n",amplitude[cnt][6],phase[cnt][6]);
    fprintf(out,"Mtm   %12.6lf %12.6lf\n",amplitude[cnt][8],phase[cnt][8]);
    fprintf(out,"Msqm  %12.6lf %12.6lf\n",amplitude[cnt][7],phase[cnt][7]);
 }
  fclose(out);
}


void usage()
{
  printf("\n\n");
  printf("################################################################################################################################\n#\n");
  printf("#  This is the FES2004 extraction and prediction software : version 1.6.0\n#\n");
  printf("#  Usage : FES2004.exe [parameters] [option] file ...\n");
  printf("#  \n");
  printf("#  ----PARAMETERS---- :\n");
  printf("#  ------------------ \n");
  printf("#  -t [argument]    Process type           : extraction - prediction  # Note : this option can be use two times ...\n");
  printf("#  -s [argument]    Data set for process   :    tide    -  loading    # Note : this option can be use two times ...\n");
  printf("#  -A [FILE]        Necessary (with -s tide) path and filename of the --TIDE-- FES netcdf data file \n");
  printf("#  -B [FILE]        Necessary (with -s load) path and filename of the --LOAD-- FES netcdf data file \n");
  printf("#  -----VARIABLE PARAMETER---- :\n");
  printf("#  --------------------------- \n");
  printf("#  -r [value]       Necessary for prediction :\n");
  printf("#                   This is the time refence of the input file.\n");
  printf("#                   Use : 0 --> 01Jan1950  0H00\n");
  printf("#                         1 --> 01Jan1958  0H00\n");
  printf("#                         2 --> 01Jan1985  0H00\n");
  printf("#                         3 --> 01Jan2000  0H00\n");
  printf("#                         DEFAULT is 0\n");
  printf("#  ----OPTIONAL PARAMETER---- :\n");
  printf("#  -------------------------- \n");
  printf("#  -o [file]        Optionnal give a specific root output file name\n");
  printf("#  \n");
  printf("#  \n");
  printf("#  \n");
  printf("#  ----EXEMPLES---- : \n");
  printf("#  ----------------  \n");
  printf("#  \n");
  printf("#    1- FES2004.exe -t extraction -s loading -B load.nc inputFILE             --> extraction of the loading data in out_file.loading.extract\n");
  printf("#  \n");
  printf("#    2- FES2004.exe -t extraction -t prediction -s tide -A tide.nc inputFILE  --> extraction of the tide data in out_file.tide.extract\n");
  printf("#                                                                                 and prediction of the tide data in out_file.tide.pred\n");
  printf("#  \n");
  printf("#    3- FES2004.exe -t prediction -s tide -s loading -A tide.nc -B load.nc inputFILE     --> prediction of both the tide and loading data\n");
  printf("#                                                                                            in out_file.tide.pred and out_file.loading.extract\n");
  printf("#  \n");
  printf("#    4- FES2004.exe -t prediction -A tide.nc -s tide -o my_output inputFILE   --> prediction of the tide data in my_output.tide.pred\n");
  printf("#  \n");
  printf("#    5- FES2004.exe -t prediction -s tide -A tide.nc -d /MYPATH/  inputFILE   --> prediction of the tide data (FES2004.nc) found in /MYPATH/\n");
  printf("#  \n");
  printf("################################################################################################################################\n\n");

  exit(9);
}

/*--------------------------------------------------------------------------*/
/*                                 MAIN 				    */ 
/*                                 ----  				    */
/*--------------------------------------------------------------------------*/


int main(int argc, char *argv[])
{

  int rstatus,n,cnt;
  double *lat,*lon,*time,*prediction;
  int nb_position,nb_CPU;
  FILE *in;
  char *keyword;
 
  int time_reference=99;
  char *netcdf_filename=NULL,*output_filename=NULL,*input_filename=NULL,*TIDE_cdf_filename=NULL,*LOAD_cdf_filename=NULL;
  int conf_extract=0, conf_predict=0 , conf_tide=0,conf_loading=0;

  double **amplitude,**phase;

  keyword=malloc(56*sizeof(char));
  netcdf_filename=malloc(256*sizeof(char));
  if(argc==1) usage();;

  n=1;
  while (n < argc)
    {
      keyword=strdup(argv[n]);
      switch (keyword[0])
	{
	case '-':
	  {
	    switch (keyword[1])
	      {
	      case 't' :
		if(strcmp(argv[n+1],"extraction")==0)conf_extract=1;
		else if(strcmp(argv[n+1],"prediction")==0)conf_predict=1;
		else {printf("=======>> unknown option %s\n",argv[n+1]);exit(99);}
		n++;
		n++;
		break;

	      case 's' :
		if(strcmp(argv[n+1],"tide")==0)conf_tide=1;
		else if(strcmp(argv[n+1],"loading")==0)conf_loading=1;
		else {printf("=======>> unknown option %s\n",argv[n+1]);exit(99);}
		n++;
		n++;
		break;

	      case 'o' :
		output_filename=strdup(argv[n+1]);
		n++;
		n++;
		break;

	      case 'r' :
		time_reference=atoi(argv[n+1]);
		n++;
		n++;
		break;

	      case 'A' :
		TIDE_cdf_filename= strdup(argv[n+1]);
		n++;
		n++;
		break;

	      case 'B' :
		LOAD_cdf_filename= strdup(argv[n+1]);
		n++;
		n++;
		break;

	      case 'h' :
		usage();
		break;

	      default:
		printf("=======>> unknown option %s\n",keyword);
		exit(98);
	      }
	    break;
	  }
	default :
	  {
	    input_filename= strdup(argv[n]);
	    n++;	
	    break;
	  }
	  free(keyword);
	}
    }
  if(input_filename==NULL)
    {
      printf("=======>> please give an input file \n");
      printf("=======>> use FES2004.EXE -help \n");
      exit(34);
    }

  if((conf_predict==1)&&(time_reference==-99) )
    {
      printf("=======>> You have not give a time reference use the -r option \n");
      printf("=======>> use FES2004.EXE -help \n");
    }

  if(output_filename==NULL)
    {
      printf("=======>> the outputs will be write in files named (type.process) (example tide.pred)\n");
      output_filename=calloc(56,sizeof(char));
      sprintf(output_filename,"out_file");
    }

 if(time_reference==99)
    {
      printf("=======>> No time ference was given --> take option 0 (01/01/1950) as reference\n");
      time_reference=0;
    }

  
  //  generality


  nb_position=compter_lignes(input_filename);
  printf("found %d lines in input file : %s\n",nb_position,input_filename);
  alloc_vectors( nb_position, &lat,&lon,&time,&prediction,&amplitude,&phase);
  printf("done --> allocation to load the input file\n");


  in=fopen(input_filename,"r");
  load_input_file(in,nb_position , lat ,lon, time);
  fclose(in);
  printf("done --> load the input file\n");


  //  conf dependent



  if(conf_predict==1)
    {
      if(conf_loading==1)
	{
	  if(LOAD_cdf_filename!=NULL)sprintf(netcdf_filename,"%s",LOAD_cdf_filename);
	  else {printf("=======>> The load Netcdf file was not given\n");usage();}
	  rstatus=fes2004_prediction(netcdf_filename,time_reference,nb_position,lat,lon,time,prediction,1);
	  sprintf(output_filename,"loading.pred");
	  rstatus=prediction_output(nb_position,prediction,lat,lon,time,time_reference,output_filename);
	}
      if(conf_tide==1)
	{
	  if(TIDE_cdf_filename!=NULL)sprintf(netcdf_filename,"%s",TIDE_cdf_filename);
	  else {printf("=======>> The tide Netcdf file was not given do not forget the -A option !!!!\n");usage();printf("=======>> The tide Netcdf file was not given do not forget the -A option !!!!\n");}
	  printf("entering in the prediction function ...\n");
	  rstatus=fes2004_prediction(netcdf_filename,time_reference,nb_position,lat,lon,time,prediction,1);
	  sprintf(output_filename,"tide.pred");
	  rstatus=prediction_output(nb_position,prediction,lat,lon,time,time_reference,output_filename);
	}
      if((conf_loading==0)&&(conf_tide==0))
	{
	  printf("=======>> you have not give the data set type for prediction, use the -s option\n");
	  printf("=======>> use FES2004.EXE -help \n");
	  exit(66);
	}
    }

  if(conf_extract==1)
    {
      if(conf_loading==1)
	{
	  if(LOAD_cdf_filename!=NULL)sprintf(netcdf_filename,"%s",LOAD_cdf_filename);
	  else {printf("=======>> The load Netcdf file was not given\n");usage();}
	  rstatus=fes2004_extraction(netcdf_filename,nb_position ,lat,lon,amplitude,phase,1);
	  sprintf(output_filename,"loading.extract");
	  rstatus=extraction_output(nb_position,amplitude,phase ,lat,lon,output_filename);	    
	}
      if(conf_tide==1)
	{
	  if(TIDE_cdf_filename!=NULL)sprintf(netcdf_filename,"%s",TIDE_cdf_filename);
	  else {printf("=======>> The tide Netcdf file was not given\n");usage();}
	  rstatus=fes2004_extraction(netcdf_filename,nb_position ,lat,lon,amplitude,phase,1);
	  sprintf(output_filename,"tide.extract");
	  rstatus=extraction_output(nb_position,amplitude,phase ,lat,lon,output_filename);	    
	}
      if((conf_loading==0)&&(conf_tide==0))
	{
	  printf("=======>> you have not give the data set type for extraction, use the -s option\n");
	  printf("=======>> use FES2004.EXE -help \n");
	  exit(66);
	}
    }

}
