//
//
// File author(s):  <Gabriel Rucabado and Antonio Garcia Dopico>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _SOLVER_IM_GPU
#define _SOLVER_IM_GPU

#include "configuration.h"
#include "discretization.h" 
#include "solver.h"

// Includes from cuda
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>

#define NUM_STREAM 1 // 1 Update 1 calculate 0 up-cal
#define NUM_EVENT 0  // 0 Update 0 calculate 0 up-cal

#define MAX_ITER this->num_elem
#define ITER_CHECK 30


class SolverImGPU : public Solver {
 private:
  //CUBLAS & CUSPARSE
  cublasHandle_t cublasHandle;
  cusparseHandle_t cusparseHandle;
  cusparseMatDescr_t descA;

  // vectors
  TYPE *R, *R0, *P, *AP, *AS;  	

  // escalars
  TYPE *p0, *p1, *apr,*beta, *res, *alpha,		//alpha[2] (0=>pos, 1=>neg)
    *oneNeg, *onePos, *zero, *w1, *w2, *w;	//w[2]     (0=>pos, 1=>neg)
  
  

  // CUDA
  dim3 dimGrid_HH;
  dim3 dimBlock_HH;
  dim3 dimGrid_ALL;
  dim3 dimBlock_ALL;
  dim3 dimGrid_DIV;
  dim3 dimBlock_DIV;
  dim3 dimGrid_DIVMUL;
  dim3 dimBlock_DIVMUL;
  cudaStream_t stream[NUM_STREAM];
  cudaEvent_t event[NUM_EVENT];
  int idGPU;

  TYPE dT;
  
  // TERMINAL
  int num_term;

  // ALL
  int num_elem;
  TYPE * pot;
  TYPE * dx;
  TYPE * c_m;
  TYPE * r_a;
  TYPE * K;
  TYPE * input_off;
  TYPE * input_on;
  TYPE * input;

  // HH
  int num_HH;
  int * indVec_DL_HH;
  int * indVec_HH;
  int * indMat_HH;
  short * type_HH;
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
  TYPE * W;

  // Implicit
  int nnz;
  int * csr_Col_Ind;
  int * csr_Row_Ptr;
  TYPE * csr_Val;

  void prevBoundaries (discretization **elements,
		       int *CPU_indMat_HH, int *CPU_indVec_DL_HH, short *CPU_type_HH,
		       TYPE *CPU_csr_Val, int *CPU_csr_Col, int * CPU_csr_Row);
  void coreMatrix (discretization **elements, 
		   int *CPU_indMat_HH, int *CPU_indVec_DL_HH, short *CPU_type_HH,
		   TYPE *CPU_csr_Val, int *CPU_csr_Col, int * CPU_csr_Row);
  void postBoundaries (discretization **elements,
		       int *CPU_indMat_HH, int *CPU_indVec_DL_HH, short *CPU_type_HH,
		       TYPE *CPU_csr_Val, int *CPU_csr_Col, int * CPU_csr_Row);
 public:
  SolverImGPU(TYPE dT, TYPE **potential, 
	      int num_term, int *ind_terminals, TYPE input_val,
	      int num_HH_elements, Hodgkin_Huxley **HH_elements,
	      int num_elem, discretization **element);
  virtual void update (bool input_active) throw();
  virtual void calculate() throw();
  virtual void snapshot(TYPE *potential) throw();
  virtual ~SolverImGPU();
};

#endif
