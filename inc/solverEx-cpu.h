//
//
// File author(s):  <Julian Andres Garcia Grajales and Gabriel Rucabado>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _SOLVER_EX_CPU
#define _SOLVER_EX_CPU

#include "discretization.h" 
#include "solver.h"

class SolverExCPU : public Solver {
 private:
  TYPE dT;
  TYPE * __restrict__ currT, *  __restrict__ prevT;

  int numStepsX;
  discretization **element;

  int num_terminals;
  int *terminals;

  int num_HH_elements;
  Hodgkin_Huxley **HH_elements;

  TYPE * input_off;
  TYPE * input_on;
  TYPE * input;


 public:
  SolverExCPU(TYPE dT, TYPE **potential, 
	    int num_term, int *ind_terminals, TYPE input_val,
	    int num_HH_elements, Hodgkin_Huxley **HH_elements,
	    int num_elem, discretization **element);
  virtual void update (bool input_active) throw();
  virtual void calculate() throw();
  virtual void snapshot(TYPE *potential) throw();
  virtual ~SolverExCPU();
};

#endif
