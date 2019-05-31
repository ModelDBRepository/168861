////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// $$    $  $$$$  $   $  $$$$$  $ $$$$$$$ $$$$  ///////////////////////////
/////////////////////////// $ $   $  $     $   $  $   $  $    $    $     ///////////////////////////
/////////////////////////// $  $  $  $$$$  $   $  $$$$$  $    $    $$$$  ///////////////////////////
/////////////////////////// $   $ $  $     $   $  $  $   $    $    $     ///////////////////////////
/////////////////////////// $    $$  $$$$  $$$$$  $   $  $    $    $$$$  ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * \mainpage NEURITE
 * \brief Neurite electrical signal propagation under mechanical loading in sequential and parallel simulations.
 * \author Julián Andrés García Grajales
 * \author Jose María Peña
 * \author Antoine Jérusalem
 * \copyright (C) 2014  This software is licenced under the terms stipulated in the license.txt file located in its root directory
 */

//
//
// File author(s):  <Julian Andres Garcia Grajales, Jose Maria Pena and Antoine Jerusalem>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!
 * \file Neurite.cpp
 * \brief This is the main file for Neurite.
 * \details The time looping is in this file.
 */
#include "configuration.h"
#include "neurite.h"  
#include "graphics.h"
#include "init.h"
#include "input.h"
#include "coupling.h"
#include "discretization.h"
#include "mechanical_model.h"
#include "space.h"
#include "GL/glut.h"
#include <time.h>
#include "solver.h"

#if defined GPU
#include "cudaException.h"
#include "solverEx-gpu.h"
#include "solverIm-gpu.h"

#elif defined CPU
#include "solverEx-cpu.h"
#include "solverIm-cpu.h"
#else
#error "You did not choose any processor"
#endif

#if (TYPE_ID == DOUBLE_ID)
#pragma message "Type ==> DOUBLE"

#elif (TYPE_ID == FLOAT_ID)
#pragma message "Type ==> FLOAT"

#else
#error "You are using an undefined TYPE_ID. Only available double or float"
#endif

// Global variables
int numStepsX=0,numStepsT=0;/*!< Total number of elements in the spatial and time discretizations respectively*/
TYPE xspan=0., tspan=0., dist_mes=0;
; /*!< Total length and total time respectively*/
const TYPE pi = 3.14159265;
int *terminals; /*!< This array contains the element numbers for the boundary conditions.*/
int num_terminals; /*!< This is the number of terminals of the tree. The first element (0) is also taken into account as a terminal element.*/
TYPE **potential, *potentialR; /*!< Two arrays to save the potential.*/
TYPE *open_GL_graphE,**open_GL_graphR,***open_GL_graph; /*!< Pointers related to the graphics part.*/
char *log_file, *output; /*!< Output files*/
FILE *standard, *study; /*!< Output files*/
TYPE unit_x=0.;
int open_GL_steps=1000; /*!< Frequency to save the data*/
char *simulation;/*!< This value defines the type of simulation*/
/*! \brief This is the main function of Neurite
 */
int main(int argc, char **argv)
{
  if (argc<5 || !strcmp(argv[1],"help"))
    { 
      if(strcmp(argv[1],"help")){
	_INFO("Your number of arguments is not enough. It has to be at least 5 (type man doc/./neurite.7 )");
      }else{
	_INFO("You have some examples for running a simulation. More information you have to type man ./Neurite");
      }
      _INFO("The solver:");
      _INFO("e: explicit mode");
      _INFO("i: implicit mode");
      _INFO("Examples:");
      _INFO("%s","./Neurite e dtf 0.9 IRE 10 NRE 5 out output.txt log simulation.log");
      _INFO("%s","./Neurite i out output.txt log simulation.log");
      exit(0);
    }
  
  //Building the simulation:
  TYPE  dT=0.,time_run=0., temp=0.;
  int i=0, k=0;
  char *args;
  input in;
  params_t params = {false};
  int GL=0, GL_cont=0, interval=0;
  time_t t1, t2;
  discretization **element; // Object oriented. Class mother = discretization with two daughters = Cable and Hodgkin-Huxley
  Hodgkin_Huxley **HH_elements; // One array for each kind of elements that Neurite has
  Cable **CT_elements;
  int num_HH_elements=0;
  int num_CT_elements=0;
  TYPE E_MAX=0;
  TYPE monitoring=0;
  int dist_mes_points=0;
  
  simulation = getenv("simulation");
  if (simulation == NULL){
    _ERROR("ERROR. You did not define the specific case to simulate\n");
  }
  t1=time(NULL);

  // Initializing some variables. Default values
  numStepsX = 0; // Initial Total number of elements;
  numStepsT = 0; // Initial Total number of time steps;
  in.IRE=0;
  in.NRE =0;
  xspan =0.;
  tspan=0.;
  in.epsilon = 0.;   // Axial strain default
  in.nu=0.00261813; // Default value for an equilibrated distribution of the strain
  in.dt_factor = 1.;  //dt_factor_default
  in.init_time = 0.;
  in.E_MAX = 1;   // Maximum microscopic strain

  // handle command-line options
  args = new char[argc]();
  inputArgs(argc, argv, args, in);
  E_MAX = in.E_MAX;
  study = fopen((const char *)&output[0],"w");
  standard = fopen((const char *)&log_file[0],"w");
  args[argc-1] = '\0';
  flags_init(args, &params); // Obtaining the flags
  test_logical_options(&params,in); // testing logical options
   
  /* Here the discretization is defined. You have to define your own discretization scheme. There are some examples in the space.cpp file*/
  inputInit(in, &params, element, HH_elements, num_HH_elements, CT_elements, num_CT_elements); 
 // Strain for the electrical model
  mechanical_model_initialization(in, element, dist_mes);
  // This way of extracting the results
  dist_mes_points=measurement_point(element,dist_mes);
  // Writing some information in the log file
  fprintf(standard,"in.nu=%f\n", in.nu);
  fprintf(standard,"in.ep=%f\n", in.epsilon);
  fprintf(standard," Total time %g\n", tspan);
  fprintf(standard,"REAL The total length is =%g, dis_mes=%g\n",xspan, dist_mes);
  for(k=1;k<numStepsX-1;k++){
    fprintf(standard,"element=%i type=%i REAL epsilon =%.10f, epsilon=%g\n",k,element[k]->get_kind_HH(),in.epsilon, element[k]->get_epsilon());    
    fprintf(standard," element size = %g \n", element[k]->get_dX());
  }

  // Allocating memory
  neuriteInit(&params);
  electrical_properties(HH_elements, num_HH_elements, E_MAX);
  setting_rest_pot(potential, element);
  hhIconds(potential[0],HH_elements, num_HH_elements, E_MAX);
  cable_init_constants(CT_elements, num_CT_elements, params.HH_model);
  critical_time_step_class(HH_elements, num_HH_elements,CT_elements, num_CT_elements);
  dT = critical_time_step_selection(element);
  fprintf(standard,"stability dT=%g\n",dT);
  dT = dT*in.dt_factor;
  fprintf(standard,"applied factor %f dT=%g\n",in.dt_factor,dT);
  fprintf(standard,"total_time t span %f\n",tspan);
  temp = tspan/dT;
  numStepsT = int(temp);
  _INFO("Total time Steps=%u",numStepsT);
 
  // For open_GL
  interval = int(numStepsT/open_GL_steps);
  fprintf(standard,"intervals =%i\n", interval);
  fprintf(standard,"TotalTime=%g,Total_length=%g\n",numStepsT*dT,xspan);
  fprintf(standard,"NumStepsT=%i,NumStepsX=%i\n",numStepsT,numStepsX);
  Solver *solver;
#if defined GPU
  try {  
    _INFO("The solver is being created");
    if (params.explicit_mode){
      solver = new SolverExGPU(dT,potential,num_terminals,terminals, in.input_current,
			       num_HH_elements, HH_elements, numStepsX, element);
    }else if (params.implicit_mode) {
      solver = new SolverImGPU(dT,potential,num_terminals,terminals, in.input_current,
				 num_HH_elements, HH_elements, numStepsX, element);
    }else{
      _ERROR("You did not define the processor");
      exit(-1); 
    }

#elif defined CPU
      if (params.explicit_mode){
	solver = new SolverExCPU(dT, potential, num_terminals, terminals, 
				  in.input_current, num_HH_elements, HH_elements,
				  numStepsX, element);
      }else if (params.implicit_mode){
	solver = new SolverImCPU(dT, potential, num_terminals, terminals, 
				  in.input_current, num_HH_elements, HH_elements,
				  numStepsX, element);
      }else{
	_ERROR("ERROR: Choose explicit or implicit)");
	exit(-1); 
      }
#else
	_ERROR("You did not define the processor");
	exit(-1); 
#endif

        
      _INFO("The simulation begins...");
      int printing=0;
      int times_printing = 20; // The times you monitor the simulation in the stdout 
      monitoring=numStepsT/times_printing;

      for(i = 1; i < numStepsT-1; i++) {

	printing++;
	time_run = i * dT;

	// Updating variables
	if ((time_run >= in.init_time) && (time_run <= in.end_input)) {
	  solver->update(true);
	}else{
	  solver->update(false);
	}

	solver->calculate();

	// OUTPUT of the simulations
	// You have to be careful about how extract the results
	// Copy for open_GL
	if(++GL_cont > interval){
	  GL_cont=0;
	  solver->snapshot(open_GL_graph[0][GL++]);
	  solver->snapshot(potential[1]);
	  if(strcmp(simulation,"Ion_Channel") == 0){
	    fprintf(study,"%g %g\n",i*dT, potential[1][1]);
	  }else {
	    fprintf(study,"%.20f %.20f\n",time_run, potential[1][dist_mes_points]);
       	  }
	}
	if(printing>monitoring){
	  _INFO_REP("%.2f%% [Time = %g]",double((100*i/numStepsT)), time_run);	   
	  printing=0;
	}
      }// DT end
      _INFO_REP("%.2f%% [Time = %g]",double((100*i/numStepsT)), time_run);	   
      t2 = time(NULL);
      _INFO("\n The simulation has ended. Real_time=%f seconds", difftime(t2,t1));
      
      // Initiate and Run Glut
      solver->snapshot(potential[1]);
      if(params.openGL){
	i=0;
	unit_x = xspan/10;
	glutInit(&argc, argv);
	runGL();
      }
            
      delete[] open_GL_graphE;
      delete[] open_GL_graphR;
      delete[] open_GL_graph;
      delete[] potential;
      delete[] potentialR;
      deleting_discretization(element, HH_elements, num_HH_elements, CT_elements, num_CT_elements);
      delete[] element;
      delete[] terminals;
      delete[] HH_elements;
      delete[] CT_elements;

      delete solver;      
      delete[] args;
      fclose(standard);

#if defined GPU
    } //try
    catch (cudaException &e) {
      e.print();
    }
#endif    
    return 0;
  }
