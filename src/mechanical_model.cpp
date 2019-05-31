//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


/*!\file mechanical_model.cpp
  \brief This file is related to the mechanical model
*/

#include "neurite.h"
#include "mechanical_model.h"
  /*! \brief This function calculates the macroscopic strain to know where it takes the measurements
    @param epsilon_axial_membrane_t Microscopic axial strain
    @param epsilon_macro_t Macroscopic axial strain
  */

extern void calculate_epsilon_macro(TYPE &epsilon_axial_membrane_t, TYPE &epsilon_macro_t){

  TYPE epsilon_macro_t0=0;
  TYPE epsilon_axial_membrane_t0=0;
  TYPE epsilon_d_axial_membrane=0;
  TYPE constant=0;
  TYPE K=1; // DEFAULT value

  TYPE max=0, min=0;

  FILE *myfile;
  TYPE myvariable;
  int i;
  int j, HEIGHT =1, WIDTH=5, CASE=0;
  TYPE bandwidth[HEIGHT][WIDTH];
  char fname[250];
  char *execution_folder = getenv("working_directory");
  if (execution_folder != NULL){
    strcpy (fname, execution_folder);
    strcat (fname, "/input_data.txt");
    myfile=fopen(fname, "r");
  }else{
    _ERROR("ERROR. You did not define the working directory variable\n");
    exit(0);
  }
  if(myfile== NULL){
        _ERROR("ERROR: You forgot the input_data.txt file with the mechanical model parameters. See mechanical_model.cpp file. As all input files for Neurite, the input data file must be in your working directory\n");
	exit(0);	
  }

  for(i = 0; i < HEIGHT; i++)
  {
    for (j = 0 ; j < WIDTH; j++)
    {
      if(!fscanf(myfile,"%lf",&myvariable)){
	_ERROR("ERROR: The scanning process in the mechanical model file was wrong\n");
	exit(0);
      };
      printf("%.5f ",myvariable);
      bandwidth[i][j] = myvariable;
    }
    printf("\n");
  }
  // BMMB Paper
  /* CASE 1-> 25 SLOW
          2-> 50 SLOW
	  3-> 100 SLOW
	  4-> 25 FAST
	  5-> 50 FAST
	  6-> 100 FAST
   */
  CASE = int(bandwidth[0][0]+0.1);
  K =bandwidth[0][1];
  epsilon_macro_t0 = bandwidth[0][2];
  epsilon_axial_membrane_t0 = bandwidth[0][3];
  max =epsilon_axial_membrane_t0;
  epsilon_d_axial_membrane = bandwidth[0][4];
  min = epsilon_d_axial_membrane;

  if(CASE<1 || CASE>6){
    _ERROR("ERROR: CASE=%i is not a case. See mechanical_model.cpp\n",CASE);
    exit(0);	
  }

  if(K>1 || K<0){
    _ERROR("ERROR: K=%f must be between 0 and 1\n",K);
    exit(0);	
  }
  if(epsilon_axial_membrane_t> max || epsilon_axial_membrane_t<min){
    _ERROR("ERROR: case %i for the %f of macroscopic strain at t0. You are out of your bandwidth: epsilon_ma=%f and you max=%f and min=%f\n", CASE,epsilon_macro_t0*100,epsilon_axial_membrane_t,max,min);
    exit(0);
  }

  fclose(myfile);

  constant= epsilon_macro_t0- K*epsilon_axial_membrane_t0;
  epsilon_macro_t= K*epsilon_axial_membrane_t+constant;

}
