//
//
// File author(s):  <Gabriel Rucabado and Antonio Garcia Dopico>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


#ifndef _SOLVER_GPU
#define _SOLVER_GPU

#include "discretization.h" 
#include "solver.h" 

// Includes from cuda
#include <cuda_runtime.h>
#include <cublas_v2.h>
  
class SolverExGPU : public Solver {
  private:
  // CUDA
  dim3 dimGrid_all;
  dim3 dimBlock_all;
  dim3 dimGrid_term;
  dim3 dimBlock_term;
  dim3 dimGrid_branch;
  dim3 dimBlock_branch;
  dim3 dimGrid_HH;
  dim3 dimBlock_HH;
  cudaStream_t stream[2];
  cudaEvent_t evento[2];
  int idGPU;
  TYPE dT;
  
  // TERMINALS 
  int num_term;
  int * ind_term;
  int * ind_mo_term;

  // BRANCHING
  int num_branch;
  int * ind_DL;
  int * ind_branch;
  int * ind_mo_branch;

  // ALL
  int num_elem;
  TYPE * potential[2];
  TYPE * new_Pot, *  old_Pot;
  TYPE * dx;
  TYPE * c_m;
  TYPE * r_a;
  TYPE * W;
  TYPE * K;
  TYPE * input_off;
  TYPE * input_on;
  TYPE * input;

  // HH
  int num_HH;
  int * ind_HH;
  TYPE * leftShiftNa;
  TYPE * leftShiftK;
  TYPE * n;
  TYPE * m;
  TYPE * h;
  TYPE * E_Na;
  TYPE * G_Na;
  TYPE * g_Na;
  TYPE * E_L;
  TYPE * G_L;
  TYPE * E_K;
  TYPE * G_K;
  TYPE * g_K;
  TYPE * rest_pot;

 public:
  SolverExGPU (TYPE dT, TYPE **potential, 
	     int num_term, int *ind_terminals, TYPE intpu_val,
	     int num_HH_elements, Hodgkin_Huxley **HH_elements,
	     int num_elem, discretization **element) throw();
  virtual void update (bool input_active) throw();
  virtual void calculate() throw();
  virtual void snapshot(TYPE *potential) throw();
  ~SolverExGPU();
};

#endif
