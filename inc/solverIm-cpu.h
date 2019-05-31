//
//
// File author(s):  <Julian Andres Garcia Grajales and Gabriel Rucabado>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


#ifndef _SOLVER_IM_CPU
#define _SOLVER_IM_CPU

#include "discretization.h" 
#include "configuration.h" 
#include "solver.h"
#include "mkl_rci.h"
#include "mkl_blas.h"
#include <mkl_spblas.h>
#include "mkl_service.h"


class SolverImCPU : public Solver {
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

  // Implicit
  int nnz;
  TYPE * csr_Val;
  MKL_INT * csr_Col_Ind;
  MKL_INT * csr_Row_Ptr;

  // Solver for sparse matrix
  void fgmres(int N, TYPE *rhs, TYPE *yoacsr,
		     MKL_INT *yoja, MKL_INT *yoia);
  void createMatrix();
  
 public:
  SolverImCPU(TYPE dT, TYPE **potential, 
	       int num_term, int *ind_terminals, TYPE input_val,
	       int num_HH_elements, Hodgkin_Huxley **HH_elements,
	       int num_elem, discretization **element);
  virtual void update (bool input_active) throw();
  virtual void calculate() throw();
  virtual void snapshot(TYPE *potential) throw();
  virtual ~SolverImCPU();
};

#endif
