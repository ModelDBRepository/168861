//
//
// File author(s):  <Gabriel Rucabado and Antonio Garcia Dopico>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

//!\file solverIm-gpu.cpp
//  \brief The implicit gpu solver of Neurite.
#include <math.h>
#include "solverIm-gpu.h"
#include "discretization.h"
#include "cudaException.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
using namespace std;


// Includes from cuda
#include "cuda_runtime.h"
#include "cublas_v2.h"

#define TYPE double


#define NORMAL 1
#define BRANCH 2
#define TERM 4
#define NEXT_TERM 8
#define FIRST_NORMAL 16
#define FIRST_BRANCH 32
#define OTHER 64

//#######################################
//#### HEADERDS of CUDA KERNELS #########
//#######################################
__global__ void Update_HH (int numElem, TYPE *m, TYPE *h, TYPE *n, TYPE dT,
			   TYPE *oldPot, TYPE *restPot, int *indVec, 
			   TYPE *G_L, TYPE *E_Na, TYPE *E_K, TYPE *E_L, 
			   TYPE *g_Na, TYPE *g_K, TYPE *G_Na, TYPE *G_K, 
			   TYPE *leftShiftNa, TYPE *leftShiftK,
			   TYPE *W, TYPE *K);
__global__ void Update_Mat (int numElem, short *type, TYPE *csrVal, TYPE dT,
			    TYPE *r_a, TYPE *c_m, TYPE *dx, TYPE *W,
			    int *indVec_DL, int *indVec, int *indMat);
__global__ void Update_B (int numElem, TYPE *K, TYPE *input, TYPE dT,
			  TYPE *c_m, TYPE *dx, TYPE *oldPot, TYPE *b);
__global__ void div (TYPE *op1, TYPE *op2, TYPE *res);
__global__ void divMul (TYPE *num1, TYPE *denom1, TYPE *num2, TYPE *denom2, TYPE *res);



//#######################################
//#### BODY of CUDA KERNELS #############
//#######################################
__global__ void Update_HH (int numElem, TYPE *m, TYPE *h, TYPE *n, TYPE dT,
			   TYPE *pot, TYPE *restPot, int *indVec, 
			   TYPE *G_L, TYPE *E_Na, TYPE *E_K, TYPE *E_L, 
			   TYPE *g_Na, TYPE *g_K, TYPE *G_Na, TYPE *G_K, 
			   TYPE *leftShiftNa, TYPE *leftShiftK,
			   TYPE *W, TYPE *K){
  TYPE dt	= dT * 1000;
  int id	= blockIdx.x * blockDim.x + threadIdx.x;
  if (id < numElem) {
    TYPE temp_G_Na, temp_G_K;
    TYPE mloc, nloc, hloc;
    TYPE Am, Ah, An, Bm, Bh, Bn; 
    int indexVec  = indVec[id];
    TYPE V_prevNa = ((pot[indexVec] + leftShiftNa[id]) - restPot[id]) / SCALE; 
    TYPE V_prevK  = ((pot[indexVec] + leftShiftK[id]) - restPot[id]) / SCALE; 

    Am 	= (2.5 - 0.1 * V_prevNa) / (exp(2.5 - 0.1 * V_prevNa) - 1);
    Bm 	= 4 * exp(-V_prevNa / 18);
    Ah 	= 0.07 * exp(-V_prevNa / 20);
    Bh	= 1 / (exp(3 - 0.1 * V_prevNa) + 1);
    An	= (0.1 - 0.01 * V_prevK) / (exp(1 - 0.1 * V_prevK)-1);
    Bn	= 0.125 * exp(-V_prevK / 80);	

    mloc	= m[id];
    hloc	= h[id];
    nloc	= n[id];

    m[id] = mloc = mloc + dt * (Am * (1 - mloc) - Bm * mloc);	    
    h[id] = hloc = hloc + dt * (Ah * (1 - hloc) - Bh * hloc);	  
    n[id] = nloc = nloc + dt * (An * (1 - nloc) - Bn * nloc);

    G_Na[id]  = temp_G_Na = g_Na[id] * mloc * mloc * mloc * hloc;
    G_K[id]   = temp_G_K  = g_K[id] * nloc * nloc * nloc * nloc;
    
    W[id] 	= -(temp_G_Na + temp_G_K + G_L[id]);
    K[indexVec] = temp_G_Na * E_Na[id] + temp_G_K * E_K[id] + G_L[id] * E_L[id];
  }
}



//! \brief Update values related to HH nodes of implicit Matrix 
__global__ void Update_Mat (int numElem, short *type, TYPE *csrVal, TYPE dT,
			    TYPE *r_a, TYPE *c_m, TYPE *dx, TYPE *W,
			    int *indVec_DL, int *indVec, int *indMat){
  double alpha, betaR, betaL, betaMe, beta, delta, gamma;
  // Figure out id and index of thread 
  int id	= blockIdx.x * blockDim.x + threadIdx.x;

  // Update Matrix (HH elements)
  if (id < numElem) {
    int indexVec = indVec[id];
    int indexMat = indMat[id]; 

    // Update Matrix
    alpha = -1/(r_a[indexVec] * dx[indexVec]);
    betaR = 1/(r_a[indexVec+1] * dx[indexVec+1]);
    betaMe = -(W[id] * dx[indexVec]) + 1/(r_a[indexVec] * dx[indexVec]);
    beta = (c_m[indexVec] * dx[indexVec])/dT + betaMe + betaR;
    delta = -betaR;

    switch (type[id]){
    case (NEXT_TERM):
      csrVal[indexMat]   = alpha;
      csrVal[indexMat+1] = beta + delta;
      break;
    case (BRANCH) :
      betaL = 1/(r_a[indVec_DL[id]] * dx[indVec_DL[id]]);
      gamma = -betaL;
    
      csrVal[indexMat]   = alpha;
      csrVal[indexMat+1] = beta + betaL;
      csrVal[indexMat+2] = delta;
      csrVal[indexMat+3] = gamma;
      break;
    case (NORMAL) :
      csrVal[indexMat]   = alpha;
      csrVal[indexMat+1] = beta;
      csrVal[indexMat+2] = delta;
      break;
    case (FIRST_NORMAL) :
      csrVal[indexMat]   = alpha + beta;
      csrVal[indexMat+1] = delta;
      break;
    case (FIRST_BRANCH) :
      betaL = 1/(r_a[indVec_DL[id]] * dx[indVec_DL[id]]);
      gamma = -betaL;
      csrVal[indexMat]   = alpha + beta + betaL;
      csrVal[indexMat+1] = delta;
      csrVal[indexMat+2] = gamma;
      break;
    default:
      break;
    }
  }
}




//! \brief Update values from B vector
__global__ void Update_B (int numElem, TYPE *K, TYPE *input, TYPE dT,
			  TYPE *c_m, TYPE *dx, TYPE *pot, TYPE *b){
  //figure out id of thread
  int id	= blockIdx.x * blockDim.x + threadIdx.x;

  //figure out B
  if (id < numElem) {
    b[id] = c_m[id] * dx[id] * pot[id]/dT + K[id] * dx[id] + input[id];
  }
}


//! \brief Divide kernel
__global__ void div (TYPE *num, TYPE *denom, TYPE *res){
  register double a = *num, b = *denom;
  if (threadIdx.x == 0){
    res[0] = +(a/b);
  } else if (threadIdx.x == 1){
    res[1] = -(a/b);
  }else{
    _INFO_CUDA ("It should not be here!!!\n");
  }
}


//! \brief Divide and multiple kernel
__global__ void divMul (TYPE *num1, TYPE *denom1, TYPE *num2, TYPE *denom2, TYPE *res){
  if (threadIdx.x == 0){
    res[0] = (num1[0] / denom1[0]) * (num2[0] / denom2[0]);
  }else{
    _INFO_CUDA ("It should not be here!!!\n");
  }  
}

//#################################################
//#### BODY of PRIVATE METHODS ####################
//#################################################
void SolverImGPU::prevBoundaries (discretization **elements,
				    int *CPU_indMat_HH, int *CPU_indVec_DL_HH, short *CPU_type_HH,
				    TYPE *CPU_csr_Val, int *CPU_csr_Col, int * CPU_csr_Row){
  discretization *elem;
  double alpha, betaR, betaL, betaMe, beta, delta, gamma;

  // 0 row
  if (elements[0]->get_kind_HH()){
    CPU_indMat_HH[this->num_HH] = 0;
    CPU_indVec_DL_HH[this->num_HH] = -1;
    CPU_type_HH[this->num_HH++] = OTHER;
  }
  CPU_csr_Row[0] = 0;
  CPU_csr_Val[this->nnz] = 1;
  CPU_csr_Col[this->nnz++] = 0;

  //1 row
  elem = elements[1];

  alpha = -1/(elem->get_r_a() * elem->get_dX());
  betaR = 1/(elements[elem->get_DR()]->get_r_a()*elements[elem->get_DR()]->get_dX());
  betaL = elem->get_branching()/
    (elements[elem->get_DL()]->get_r_a()*elements[elem->get_DL()]->get_dX());
  betaMe = -elem->get_W() * elem->get_dX();
  beta = (elem->get_c_m() * elem->get_dX())/this->dT + betaMe - alpha + betaL + betaR;
  delta = -betaR;

  CPU_csr_Row[1] = 1;
  CPU_csr_Val[this->nnz] = alpha + beta;
  CPU_csr_Col[this->nnz++] = 1;
  
  CPU_csr_Val[this->nnz] = delta;
  CPU_csr_Col[this->nnz++] = 2;

  if (elements[1]->get_DL() != 0) {
    gamma = -betaL;
    CPU_csr_Val[this->nnz] = gamma;
    CPU_csr_Col[this->nnz++] = elements[1]->get_DL();
      
    if (elements[1]->get_kind_HH()){
      CPU_indMat_HH[this->num_HH] = 1;
      CPU_indVec_DL_HH[this->num_HH] = 1;
      CPU_type_HH[this->num_HH++] = FIRST_BRANCH;
    }
  } else{
    if (elements[1]->get_kind_HH()){
      CPU_indMat_HH[this->num_HH] = 1;
      CPU_indVec_DL_HH[this->num_HH] = -1;
      CPU_type_HH[this->num_HH++] = FIRST_NORMAL;
    }
  }
}



void SolverImGPU::coreMatrix (discretization **elements,
				int *CPU_indMat_HH, int *CPU_indVec_DL_HH, short *CPU_type_HH,
				TYPE *CPU_csr_Val, int *CPU_csr_Col, int * CPU_csr_Row){
  // 2<-->(n-2) Rows
  int i;
  discretization *elem, *nextElem;
  double alpha, betaR, betaL, betaMe, beta, delta, gamma;

  for (i=2; i<this->num_elem-2; i++){
    elem = elements[i];
    nextElem = elements[i+1];
    
    if (elem->get_kind_HH()){
      CPU_indMat_HH[this->num_HH] = this->nnz;
    }

    CPU_csr_Row[i] = this->nnz;

    if(elem->get_DR() == i){ 
      CPU_csr_Val[this->nnz] 	= 1;
      CPU_csr_Col[this->nnz++] = i;
      if (elem->get_kind_HH()){
	CPU_indVec_DL_HH[this->num_HH] = -1;
	CPU_type_HH[this->num_HH++] = TERM;
      }

    } else{
      alpha = -1/(elem->get_r_a() * elem->get_dX());
      betaR = 1/(elements[elem->get_DR()]->get_r_a()*elements[elem->get_DR()]->get_dX());
      betaL = elem->get_branching()/
	(elements[elem->get_DL()]->get_r_a()*elements[elem->get_DL()]->get_dX());
      betaMe = -elem->get_W() * elem->get_dX();
      beta = (elem->get_c_m() * elem->get_dX())/this->dT + betaMe - alpha + betaL + betaR;
      delta = -betaR;

      CPU_csr_Val[this->nnz] = alpha;
      CPU_csr_Col[this->nnz++] = elem->get_mother();

      if(nextElem->get_DR() == i+1){
	CPU_csr_Val[this->nnz] = beta + delta;
	CPU_csr_Col[this->nnz++] = i;

	if (elem->get_kind_HH()){
	  CPU_indVec_DL_HH[this->num_HH] = -1;
	  CPU_type_HH[this->num_HH++] = NEXT_TERM;
	}
      }else{
	CPU_csr_Val[this->nnz] = beta;
	CPU_csr_Col[this->nnz++] = i;

	CPU_csr_Val[this->nnz] = delta;
	CPU_csr_Col[this->nnz++] = elem->get_DR();

	if (elem->get_DL() != 0) {
	  gamma = -betaL;
	  CPU_csr_Val[this->nnz] = gamma;
	  CPU_csr_Col[this->nnz++] = elem->get_DL();

	  if (elem->get_kind_HH()){
	    CPU_indVec_DL_HH[this->num_HH] = elem->get_DL();
	    CPU_type_HH[this->num_HH++] = BRANCH;
	  }
	} else {
	  if (elem->get_kind_HH()){
	    CPU_indVec_DL_HH[this->num_HH] = -1;
	    CPU_type_HH[this->num_HH++] = NORMAL;
	  }
	}
      }
    }
  }
}




void SolverImGPU::postBoundaries (discretization **elements,
				  int *CPU_indMat_HH, int *CPU_indVec_DL_HH, short *CPU_type_HH,
				  TYPE *CPU_csr_Val, int *CPU_csr_Col, int * CPU_csr_Row){
  discretization *elem;
  double alpha, betaR, betaMe, beta, delta;

  // n-2 row
  elem = elements[this->num_elem-2];

  alpha = -1/(elem->get_r_a() * elem->get_dX());
  betaR = 1/(elements[elem->get_DR()]->get_r_a()*elements[elem->get_DR()]->get_dX());
  betaMe = -elem->get_W() * elem->get_dX();
  beta = (elem->get_c_m() * elem->get_dX())/this->dT + betaMe -alpha + betaR;
  delta = -betaR;


  if (elem->get_kind_HH()){
    CPU_indMat_HH[this->num_HH] = this->nnz;
    CPU_indVec_DL_HH[this->num_HH] = -1;
    CPU_type_HH[this->num_HH++] = NEXT_TERM;
  }
  CPU_csr_Row[this->num_elem-2] = this->nnz;
  CPU_csr_Val[this->nnz] = alpha;
  CPU_csr_Col[this->nnz++] = elem->get_mother();

  CPU_csr_Val[this->nnz] = beta + delta;
  CPU_csr_Col[this->nnz++] = this->num_elem-2;

  //n-1 row
  if (elements[this->num_elem-1]->get_kind_HH()){
    CPU_indMat_HH[this->num_HH] = this->nnz;
    CPU_indVec_DL_HH[this->num_HH] = -1;
    CPU_type_HH[this->num_HH++] = OTHER;
  }

  CPU_csr_Row[this->num_elem-1] = this->nnz;
  CPU_csr_Val[this->nnz] = 1;
  CPU_csr_Col[this->nnz++] = this->num_elem - 1;
  CPU_csr_Row[this->num_elem] = this->nnz;
}








//#################################################
//#### BODY of PUBLIC METHODS #####################
//#################################################
SolverImGPU::SolverImGPU(double dT, double **potential,
			   int num_term, int *ind_terminals, double input_val,
			   int num_HH_elements, Hodgkin_Huxley **HH_elements,
			   int num_elem, discretization **elements) {
  int i=0;
  cudaError_t error;
  TYPE onePos=+1, oneNeg=-1, cero=0;
  char *env = getenv("idGPU");
  cudaDeviceProp deviceProp;

  try {
    this->idGPU =(env != NULL)? atoi(env):CUDA_ID_DEF;
    if ((error = cudaSetDevice (this->idGPU)) != OK_CUDA){
      std::stringstream params;
      params << "idGPU=" << this->idGPU;
      throw cudaException(__FILE__, __FUNCTION__, __LINE__,	       
			  "Identifying the GPU", error,  params.str());     
    }

    cudaGetDeviceProperties(&deviceProp, this->idGPU);
    _INFO("GPU used ==> Device %d: \"%s\"", this->idGPU, deviceProp.name);

    // COMMON
    //////////////////////////////////////////////////
    _DEBUG("STARTING TO CREATE STRUCTURES OF SOLVER");
    DIM3(dimBlock_ALL, TH_ALL,1,1);
    DIM3(dimBlock_HH, TH_HH,1,1);
    DIM3(dimBlock_DIV, 2,1,1);
    DIM3(dimBlock_DIVMUL, 1,1,1);

    DIM3(dimGrid_DIV, 1,1,1);
    DIM3(dimGrid_DIVMUL, 1,1,1);
    this->dT 	= dT;

    // Create Events
    for (i=0; i<NUM_EVENT; i++){
      if ((error = cudaEventCreateWithFlags(&this->event[i],
					    cudaEventDisableTiming)) != OK_CUDA){
	std::stringstream params, msg;
	params << "event=" << &event[i];
	msg    << "creating event: " + i;
	throw cudaException(__FILE__, __FUNCTION__, __LINE__,
			    msg.str(), error,  params.str());
      }
    }
    // Create Streams
    for (i=0; i<NUM_STREAM; i++){
      cudaError_t error;
      if ((error = cudaStreamCreate(&stream[i]))!= OK_CUDA){
	std::stringstream params,msg;
	params << "stream=" << &stream[i];
	msg << "Creating the stream: " << i;
	throw cudaException(__FILE__, __FUNCTION__, __LINE__,	       
			    msg.str(), error,  params.str());				
      }
    }

    // Create CUBLAS context 
    cublasStatus_t cublasStatus;
    if ((cublasStatus = cublasCreate(&cublasHandle)) != OK_CUBLAS){
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,
			  "creating context cublas", cublasStatus, "");
    }
    
    if ((cublasStatus = cublasSetPointerMode(this->cublasHandle, 
					     CUBLAS_POINTER_MODE_DEVICE)) != OK_CUBLAS){
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,
			  "configurating pointers for cublas", cublasStatus, "");
    }
    
    // Create CUSPARSE context
    cusparseStatus_t cusparseStatus;
    if ((cusparseStatus = cusparseCreate(&cusparseHandle)) != OK_CUSPARSE){
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,
			  "Creating cusparse context", cusparseStatus, "");
    }
    if ((cusparseStatus = cusparseSetPointerMode(this->cusparseHandle, 
						 CUSPARSE_POINTER_MODE_DEVICE)) != OK_CUSPARSE){
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,
			  "Creating pointers for cusparse", cusparseStatus, "");
    }
    if ((cusparseStatus = cusparseCreateMatDescr(&this->descA)) != OK_CUSPARSE){
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,
			  "Creating A matrix for cusparse", cusparseStatus, "");
    }
    if ((cusparseStatus = cusparseSetMatIndexBase(this->descA, CUSPARSE_INDEX_BASE_ZERO)) != OK_CUSPARSE){
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,
			  "Creating A Matrix for cusparse", cusparseStatus, "");
    }
    if ((cusparseStatus = cusparseSetMatType(this->descA, CUSPARSE_MATRIX_TYPE_GENERAL)) != OK_CUSPARSE){
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,
			  "Creating A Matrix for cusparse", cusparseStatus, "");
    }

   
    CUDA_ALLOC(&(this->AP), num_elem * sizeof(TYPE));
    CUDA_ALLOC(&(this->AS), num_elem * sizeof(TYPE));
    CUDA_ALLOC(&(this->P), num_elem * sizeof(TYPE));
    CUDA_ALLOC(&(this->R), num_elem * sizeof(TYPE));
    CUDA_ALLOC(&(this->R0), num_elem * sizeof(TYPE));
    CUDA_SET (this->AP, 0, num_elem * sizeof(TYPE), stream[0]);
    CUDA_SET (this->AS, 0, num_elem * sizeof(TYPE), stream[0]);
    CUDA_SET (this->P, 0, num_elem * sizeof(TYPE), stream[0]);
    CUDA_SET (this->R, 0, num_elem * sizeof(TYPE), stream[0]);
    CUDA_SET (this->R0, 0, num_elem * sizeof(TYPE), stream[0]);


    CUDA_ALLOC(&(this->p0), 1 * sizeof(TYPE));
    CUDA_ALLOC(&(this->p1), 1 * sizeof(TYPE));
    CUDA_ALLOC(&(this->apr), 1 * sizeof(TYPE));
    CUDA_ALLOC(&(this->alpha), 2 * sizeof(TYPE));
    CUDA_ALLOC(&(this->w1), 1 * sizeof(TYPE));
    CUDA_ALLOC(&(this->w2), 1 * sizeof(TYPE));
    CUDA_ALLOC(&(this->w), 2 * sizeof(TYPE));
    CUDA_ALLOC(&(this->beta), 1 * sizeof(TYPE));
    CUDA_ALLOC(&(this->onePos), 1 * sizeof(TYPE));
    CUDA_ALLOC(&(this->oneNeg), 1 * sizeof(TYPE));
    CUDA_ALLOC(&(this->zero), 1 * sizeof(TYPE));
    CUDA_ALLOC(&(this->res), 1 * sizeof(TYPE));

    CUDA_CPY (this->onePos, &onePos, 1*sizeof(TYPE), stream[0]);
    CUDA_CPY (this->oneNeg, &oneNeg, 1*sizeof(TYPE), stream[0]);
    CUDA_CPY (this->zero,   &cero,   1*sizeof(TYPE), stream[0]);
    CUDA_SET (this->p0,   0, 1*sizeof(TYPE), stream[0]);
    CUDA_SET (this->p1,   0, 1*sizeof(TYPE), stream[0]);
    CUDA_SET (this->apr,  0, 1*sizeof(TYPE), stream[0]);
    CUDA_SET (this->w1,   0, 1*sizeof(TYPE), stream[0]);
    CUDA_SET (this->w2,   0, 1*sizeof(TYPE), stream[0]);
    CUDA_SET (this->beta, 0, 1*sizeof(TYPE), stream[0]);
    CUDA_SET (this->alpha,0, 2*sizeof(TYPE), stream[0]);
    CUDA_SET (this->w,    0, 2*sizeof(TYPE), stream[0]);


    // TERMINALS
    //////////////////////////////////////////////////
    this->num_term 	= num_term;


    // ALL
    //////////////////////////////////////////////////
    this->num_elem	= num_elem;
    //   POTENTIAL
    CUDA_ALLOC(&this->pot, num_elem * sizeof(TYPE));
    CUDA_CPY(this->pot, potential[0], sizeof(TYPE)*num_elem, stream[0]);

    //   ATTRIBUTES
    TYPE *CPU_dx	= new TYPE[num_elem];
    TYPE *CPU_c_m	= new TYPE[num_elem];
    TYPE *CPU_r_a	= new TYPE[num_elem];
    TYPE *CPU_K		= new TYPE[num_elem];
    TYPE *CPU_input_on	= new TYPE[num_elem];
    TYPE *CPU_input_off	= new TYPE[num_elem];
    for (i=0; i<num_elem; i++){
      CPU_dx[i]		= elements[i]->get_dX();
      CPU_c_m[i]	= elements[i]->get_c_m();
      CPU_r_a[i]	= elements[i]->get_r_a();
      CPU_K[i]		= elements[i]->get_K();
      CPU_input_off[i]	= 0;
      CPU_input_on[i]	= elements[i]->get_input_current();
    }

    DIM3(dimGrid_ALL, (this->num_elem + TH_ALL - 1)/TH_ALL,1,1);
    int num_all_cuda = dimGrid_ALL.x * dimBlock_ALL.x;

    CUDA_ALLOC(&(this->dx), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->c_m), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->r_a), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->K), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->input_on), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->input_off), num_all_cuda * sizeof(TYPE));
    CUDA_CPY(this->dx, CPU_dx, this->num_elem * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->c_m, CPU_c_m, this->num_elem * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->r_a, CPU_r_a, this->num_elem * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->K, CPU_K, this->num_elem * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->input_on, CPU_input_on, this->num_elem * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->input_off, CPU_input_off, this->num_elem * sizeof(TYPE), stream[0]);

    // HH
    //////////////////////////////////////////////////
    this->num_HH	= num_HH_elements;
    if (this->num_HH > 0){
      int *CPU_indVec_HH	= new int [num_HH_elements];
      TYPE *CPU_W		= new TYPE [num_HH_elements];
      TYPE *CPU_Shift_Na	= new TYPE [num_HH_elements];
      TYPE *CPU_Shift_K	= new TYPE [num_HH_elements];
      TYPE *CPU_n		= new TYPE [num_HH_elements];
      TYPE *CPU_m		= new TYPE [num_HH_elements];
      TYPE *CPU_h		= new TYPE [num_HH_elements];
      TYPE *CPU_E_Na	= new TYPE [num_HH_elements];
      TYPE *CPU_G_Na	= new TYPE [num_HH_elements];
      TYPE *CPU_g_Na	= new TYPE [num_HH_elements];
      TYPE *CPU_E_L	= new TYPE [num_HH_elements];
      TYPE *CPU_G_L	= new TYPE [num_HH_elements];
      TYPE *CPU_E_K	= new TYPE [num_HH_elements];
      TYPE *CPU_G_K	= new TYPE [num_HH_elements];
      TYPE *CPU_g_K	= new TYPE [num_HH_elements];
      TYPE *CPU_rest_pot	= new TYPE [num_HH_elements];
      for (i=0; i<num_HH_elements; i++){
	CPU_indVec_HH[i]	= HH_elements[i]->get_element();
	CPU_W[i]		= HH_elements[i]->get_W();
	CPU_Shift_Na[i]	= HH_elements[i]->get_Left_Shift_Na();
	CPU_Shift_K[i]	= HH_elements[i]->get_Left_Shift_K();
	CPU_n[i]		= HH_elements[i]->get_n();
	CPU_m[i]		= HH_elements[i]->get_m();
	CPU_h[i]		= HH_elements[i]->get_h();
	CPU_E_Na[i]	= HH_elements[i]->get_E_Na();
	CPU_G_Na[i]	= HH_elements[i]->get_G_Na();
	CPU_g_Na[i]	= HH_elements[i]->get_g_Na();
	CPU_E_L[i]	= HH_elements[i]->get_E_L();
	CPU_G_L[i]	= HH_elements[i]->get_G_L();
	CPU_E_K[i]	= HH_elements[i]->get_E_K();
	CPU_G_K[i]	= HH_elements[i]->get_G_K();
	CPU_g_K[i]	= HH_elements[i]->get_g_K();
	CPU_rest_pot[i]	= HH_elements[i]->get_rest_pot();
      }

      DIM3(dimGrid_HH, (this->num_HH + TH_HH - 1)/TH_HH,1,1);
      int num_HH_cuda = dimGrid_HH.x * dimBlock_HH.x;

      CUDA_ALLOC(&(this->indVec_HH), num_HH_cuda * sizeof(int));
      CUDA_ALLOC(&(this->W), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->leftShiftNa), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->leftShiftK), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->n), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->m), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->h), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->E_Na), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->G_Na), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->g_Na), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->E_L), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->G_L), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->E_K), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->G_K), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->g_K), num_HH_cuda * sizeof(TYPE));
      CUDA_ALLOC(&(this->rest_pot), num_HH_cuda * sizeof(TYPE));

      CUDA_CPY(this->indVec_HH, CPU_indVec_HH, this->num_HH * sizeof(int), stream[0]);
      CUDA_CPY(this->W, CPU_W, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->leftShiftNa, CPU_Shift_Na, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->leftShiftK, CPU_Shift_K, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->n, CPU_n, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->m, CPU_m, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->h, CPU_h, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->E_Na, CPU_E_Na, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->G_Na, CPU_G_Na, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->g_Na, CPU_g_Na, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->E_L, CPU_E_L, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->G_L, CPU_G_L, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->E_K, CPU_E_K, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->G_K, CPU_G_K, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->g_K, CPU_g_K, this->num_HH * sizeof(TYPE), stream[0]);
      CUDA_CPY(this->rest_pot, CPU_rest_pot, this->num_HH * sizeof(TYPE), stream[0]);
    }
  }catch (cudaException &e){
    e.print();
    throw e;
  }

  // Creation temporal structures that related to the matrix
  //////////////////////////////////////////////////
  TYPE *CPU_csr_Val 	= new TYPE [this->num_elem*4];
  int *CPU_csr_Col 	= new int [this->num_elem*4];
  int *CPU_csr_Row 	= new int [this->num_elem+1];
  int *CPU_indMat	= new int [this->num_HH];
  int *CPU_indVec_DL	= new int [this->num_HH];
  short *CPU_type	= new short [this->num_HH];

  this->nnz=0;
  this->num_HH=0;
  this->prevBoundaries(elements, CPU_indMat, CPU_indVec_DL,
		       CPU_type, CPU_csr_Val, CPU_csr_Col, CPU_csr_Row);
  _DEBUG("prev->coreMatrix"); 
  this->coreMatrix(elements, CPU_indMat, CPU_indVec_DL,
		   CPU_type, CPU_csr_Val, CPU_csr_Col, CPU_csr_Row);
  _DEBUG("coreMatrix->post"); 
  this->postBoundaries(elements, CPU_indMat, CPU_indVec_DL,
		       CPU_type, CPU_csr_Val, CPU_csr_Col, CPU_csr_Row);
  
  // Transferring the matrix from CPU to GPU
  CUDA_ALLOC (&(this->csr_Val), sizeof(TYPE)*this->nnz);
  CUDA_CPY(this->csr_Val, CPU_csr_Val, sizeof(TYPE)*this->nnz, stream[0]);
  CUDA_ALLOC (&(this->csr_Col_Ind), sizeof(int)*this->nnz);
  CUDA_CPY(this->csr_Col_Ind, CPU_csr_Col, sizeof(int)*this->nnz, stream[0]);
  CUDA_ALLOC (&(this->csr_Row_Ptr), sizeof(int)*(this->num_elem+1));
  CUDA_CPY(this->csr_Row_Ptr, CPU_csr_Row, sizeof(int)*(this->num_elem+1), stream[0]);

  CUDA_ALLOC (&(this->type_HH), sizeof(short)*this->num_HH);
  CUDA_CPY(this->type_HH, CPU_type, sizeof(short)*this->num_HH, stream[0]);
  CUDA_ALLOC (&(this->indMat_HH), sizeof(int)*this->num_HH);
  CUDA_CPY(this->indMat_HH, CPU_indMat, sizeof(int)*this->num_HH, stream[0]);
  CUDA_ALLOC (&(this->indVec_DL_HH), sizeof(int)*this->num_HH);
  CUDA_CPY(this->indVec_DL_HH, CPU_indVec_DL, sizeof(int)*this->num_HH, stream[0]);

  // Free allocated temporal matrix (CPU)
  free(CPU_csr_Val);
  free(CPU_csr_Col);
  free(CPU_csr_Row);
  free(CPU_indMat);
  free(CPU_indVec_DL);

  cudaDeviceSynchronize();
}


  //!\brief This function updates the electrical variables based on the previous time step information
  //  @param input_active Flag used to indicate if it has to use  the input current to update the values
void SolverImGPU::update (bool input_active) throw(){
  
  try{
    // Set GPU card to use
    cudaSetDevice (this->idGPU);
    
    if (input_active){
      input	= input_on;
    }else{
      input	= input_off;
    }

    if (num_HH > 0){  
      KERNEL(Update_HH, dimGrid_HH, dimBlock_HH, stream[0], 
     	     "%d, %p, %p, %p, %f, %p, %p, %p, %p,"
	     "%p, %p, %p, %p, %p, %p, %p, %p, %p, %p, %p",
     	     num_HH, m, h, n, dT, pot, rest_pot, indVec_HH, 
     	     G_L, E_Na, E_K, E_L, g_Na, g_K, G_Na, G_K,
	     leftShiftNa, leftShiftK, W, K);

      KERNEL(Update_Mat, dimGrid_HH, dimBlock_HH, stream[0], 
     	     "%d, %p, %p, %f, %p, %p, %p, %p, %p, %p, %p",
     	     num_HH, type_HH, csr_Val, dT, r_a, c_m, dx, W,
	     indVec_DL_HH, indVec_HH, indMat_HH);
    }

    KERNEL(Update_B, dimGrid_ALL, dimBlock_ALL, stream[0], 
	   "%d, %p, %p, %f, %p, %p, %p, %p",
	   num_elem, K, input, dT, c_m, dx, this->pot, this->pot);
  } catch (cudaException &e){
    e.print();
    throw e;
  }
}



//! \brief Calculate one iteration of explicit finite difference discretization solver
void SolverImGPU::calculate() throw(){
  /*
R0==B
csrmv (A,x)->AP;
daxpy (R0 -AP)->R0;
dot (R0, R0)->p0;
copy (R0,P)
copy (R0,R)

while (true){
=================================
|  	STREAM 0 		|
=================================
| csrmv(A,P)->AP 		|
| dot(AP,R0)->apr 		|
| alpha=p0/apr			|
| daxpy(R-alpha*AP)->R		|
| csrmv(A, R)->As  		|
| dot(As,R)->w1			|
| dot (As, As)->w2	  	|
| w=w1/w2			|
| daxpy (Xi+ alpha*P)->Xi+1	|
| daxpy (Xi+1 + w*R)->Xi+1	|
| daxpy (R-w*As)->R		| 
| dot (R, R0)->p1		|
| beta=alpha/w * p1/p0		|
| copy(P->AS)			|
| daxpy(AS-w*AP)->AS		|
| copy(R->P)			|
| daxpy(P + beta*AS)->P		|
=================================
p0=p1 (flip)
}
  */
  int i,j;
  TYPE result, *flip, limit = RANGE_ERROR*RANGE_ERROR;

  try{
    // Set GPU card to use
    cudaSetDevice (this->idGPU);

    cublasSetStream(this->cublasHandle, stream[0]);
    cusparseSetStream(this->cusparseHandle, stream[0]);

    CUDA_CPY_SYNC(this->R0, this->pot, sizeof(TYPE)*this->num_elem);
    cusparseXcsrmv(this->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
    		   this->num_elem, this->num_elem, this->nnz, this->onePos,
    		   this->descA, this->csr_Val, this->csr_Row_Ptr, this->csr_Col_Ind,
    		   this->pot, this->zero, this->AP);
    cublasXaxpy(this->cublasHandle, this->num_elem, this->oneNeg, this->AP, 1, this->R0, 1);
    CUDA_CPY_SYNC(this->R, this->R0, sizeof(TYPE)*this->num_elem);
    CUDA_CPY_SYNC(this->P, this->R0, sizeof(TYPE)*this->num_elem);
    cublasXot(this->cublasHandle, this->num_elem, this->R0, 1, this->R0, 1, this->p0);

    for (i=0; i<MAX_ITER; i++){
      for (j=i; j<i+ITER_CHECK; j++){ 
	cusparseXcsrmv(this->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		       this->num_elem, this->num_elem, this->nnz, this->onePos,
		       this->descA, this->csr_Val, this->csr_Row_Ptr, this->csr_Col_Ind,
		       this->P, this->zero, this->AP);
	cublasXot(this->cublasHandle, this->num_elem, this->AP, 1, this->R0, 1, this->apr);
	KERNEL(div, dimGrid_DIV, dimBlock_DIV, stream[0], 
	       "%p, %p, %p", this->p0, this->apr, this->alpha);
	cublasXaxpy(this->cublasHandle, this->num_elem, &this->alpha[1], this->AP, 1, this->R, 1);
	cusparseXcsrmv(this->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		       this->num_elem, this->num_elem, this->nnz, this->onePos,
		       this->descA, this->csr_Val, this->csr_Row_Ptr, this->csr_Col_Ind,
		       this->R, this->zero, this->AS);
	cublasXot(this->cublasHandle, this->num_elem, this->AS, 1, this->R, 1, this->w1);
	cublasXot(this->cublasHandle, this->num_elem, this->AS, 1, this->AS, 1, this->w2);

	KERNEL(div, dimGrid_DIV, dimBlock_DIV, stream[0], 
	       "%p, %p, %p", this->w1, this->w2, this->w);

	cublasXaxpy(this->cublasHandle, this->num_elem, &(this->alpha[0]), this->P, 1, this->pot, 1);
	cublasXaxpy(this->cublasHandle, this->num_elem, &(this->w[0]), this->R, 1, this->pot, 1);

	cublasXaxpy(this->cublasHandle, this->num_elem, &this->w[1], this->AS, 1, this->R, 1);

	cublasXot(this->cublasHandle, this->num_elem, this->R, 1, this->R0, 1, this->p1);

	KERNEL(divMul, dimGrid_DIVMUL, dimBlock_DIVMUL, stream[0], 
	       "%p, %p, %p, %p, %p",
	       &(this->alpha[0]), &(this->w[0]), this->p1, this->p0, this->beta);
	CUDA_CPY_SYNC (this->AS, this->P, sizeof(TYPE)*this->num_elem);
	cublasXaxpy(this->cublasHandle, this->num_elem, &this->w[1], this->AP, 1, this->AS, 1);

	CUDA_CPY_SYNC (this->P, this->R, sizeof(TYPE)*this->num_elem);
	cublasXaxpy(this->cublasHandle, this->num_elem, this->beta, this->AS, 1, this->P, 1);

	flip 	= this->p0;
	this->p0= this->p1;
	this->p1= flip;
      }

      i=j-1;
      // SYNC HOST-DEVICE to check error
      cublasXot(this->cublasHandle, this->num_elem, this->R, 1, this->R,1,this->res);
      //COPY SYNC
      CUDA_CPY_SYNC(&result, this->res, sizeof(TYPE));
      if (result <= limit){
	break;
      }
    }


    if (result <= limit){
      _DEBUG("Finalize by precision at iteration %d [%.16e <= %.16e]",
	     i-1, result, RANGE_ERROR);
    } else {
      _DEBUG("Finalize by iterations [%.16e > %.16e]", result, RANGE_ERROR);
    }
    checkError(__PRETTY_FUNCTION__, __BASE_FILE__, __LINE__);

  }catch (cudaException &e){
    e.print();
    throw e;
  }
}

   //! \brief Copy the current potential to the pointer of the parameter
  //  @param potential Pointer on which the potential will be copied 
void SolverImGPU::snapshot(double *potential) throw(){
  try{
    CUDA_CPY(potential, this->pot, sizeof(TYPE) * this->num_elem, stream[0]);
  }catch (cudaException &e){
    e.print();
    throw e;
  }
}


//! \brief Free all memory used by solver
SolverImGPU::~SolverImGPU(){
  CUDA_FREE(this->csr_Val);
  CUDA_FREE(this->csr_Col_Ind);
  CUDA_FREE(this->csr_Row_Ptr);

  CUDA_FREE(this->pot);
  CUDA_FREE(this->dx);
  CUDA_FREE(this->c_m);
  CUDA_FREE(this->r_a);
  CUDA_FREE(this->K);
  CUDA_FREE(this->input_on);
  CUDA_FREE(this->input_off);

  CUDA_FREE(this->indVec_DL_HH);
  CUDA_FREE(this->indVec_HH);
  CUDA_FREE(this->indMat_HH);
  CUDA_FREE(this->type_HH);
  CUDA_FREE(this->n);
  CUDA_FREE(this->m);
  CUDA_FREE(this->h);
  CUDA_FREE(this->E_Na);
  CUDA_FREE(this->G_Na);
  CUDA_FREE(this->g_Na);
  CUDA_FREE(this->E_L);
  CUDA_FREE(this->G_L);
  CUDA_FREE(this->E_K);
  CUDA_FREE(this->G_K);
  CUDA_FREE(this->g_K);
  CUDA_FREE(this->rest_pot);
  CUDA_FREE(this->W);
}
