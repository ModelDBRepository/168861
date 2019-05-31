//
//
// File authors:  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file input.cpp
  \brief Several functions to modify the values given in the command line
*/
#include "input.h"
#include "configuration.h"

/*global variables*/
extern int numStepsX, numStepsT;
extern TYPE xspan, tspan;
extern char *log_file, *output, *out_epsilons;


static void inputArgFlags(int argc, char **argv, char *args, input &in, int &i, int j);
static void inputVals(int argc, char **argv, int &i, char *args, int &inpt);
static void inputVals(int argc, char **argv, int &i, char *args, TYPE &inpt);

extern void inputArgs(int argc, char **argv, char *args, input &in)
{
	for(int i=1; i<argc; i++)
	{
		if(argv[i][0] == '-')
			inputArgFlags(argc, argv, args, in, i, 1);
		else
			inputArgFlags(argc, argv, args, in, i, 0);
	}
}

static void inputArgFlags(int argc, char **argv, char *args, input &in, int &i, int f)
{
	args[i-1] = argv[i][f];
	if(argv[i][f] == 'x'){
	  if(argv[i][f+1] == 's'){
	    if(argv[i][f+2] == 'p'){
	      if(argv[i][f+3] == 'a'){
		if(argv[i][f+4] == 'n'){
		  inputVals(argc, argv, i, args, xspan);
		}	      
	      }
	    }
	  }
	}


	if(argv[i][f] == 'T')
	if(argv[i][f+1] == 't'){
		inputVals(argc, argv, i, args, tspan);
	}

	if(argv[i][f] == 'n')
	if(argv[i][f+1] == 'u'){
		inputVals(argc, argv, i, args, in.nu);
	}

	if(argv[i][f] == 's')
	  if(argv[i][f+1] == 't'){
	    if(argv[i][f+2] == 'r'){
	      inputVals(argc, argv, i, args, in.epsilon);
	    }
	  }

	if(argv[i][f] == 'I'){
	  if(argv[i][f+1] == 'R'){
	    if(argv[i][f+2] == 'E'){
	    inputVals(argc, argv, i, args, in.IRE);
	    }
	  }
	}
	if(argv[i][f] == 'N'){
	  if(argv[i][f+1] == 'R'){
	    if(argv[i][f+2] == 'E'){
	    inputVals(argc, argv, i, args, in.NRE);
	    }
	  }
	}
	if(argv[i][f] == 'l')
	  if(argv[i][f+1] == 'o'){

	    if(argv[i][f+2] == 'g'){

	      log_file = argv[i+1];  
	      args[i] = ' ';
	      i++;
	    }
	  }

	if(argv[i][f] == 'd')
	  if(argv[i][f+1] == 't'){
	    if(argv[i][f+2] == 'f'){
	      inputVals(argc, argv, i, args, in.dt_factor);
	    }
	  }

	if(argv[i][f] == 'o')
	  if(argv[i][f+1] == 'u'){
	    if(argv[i][f+2] == 't'){
	      output = argv[i+1];  
	      args[i] = ' ';
	      i++;
	    }
	  }

	if(argv[i][f] == 'm')
	  if(argv[i][f+1] == 'e'){
	    if(argv[i][f+2] == 's'){
	      inputVals(argc, argv, i, args, in.dist_mes);
	    }
	  }

	if(argv[i][f] == 'M'){
	  if(argv[i][f+1] == 'A'){
	    if(argv[i][f+2] == 'X'){
	      if(argv[i][f+3] == 'E'){
		inputVals(argc, argv, i, args, in.E_MAX);
	      }
	    }
	  }
	}


}

static void inputVals(int argc, char **argv, int &i, char *args, int &inpt)
{
	if((i != argc-1) && ((argv[i+1][0] < 58) && (argv[i+1][0] > 47)))
	{
		inpt = atoi(argv[i+1]); 
		args[i] = ' ';
		i++;
	}
}

static void inputVals(int argc, char **argv, int &i, char *args, TYPE &inpt)
{

  if((i != argc-1) && ((argv[i+1][0] < 58) && (argv[i+1][0] > 47)))
    {
      inpt = atof(argv[i+1]);                                                                                                                
      args[i] = ' ';
      i++;
    }
}

