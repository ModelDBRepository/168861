//
//
// File author(s):  <Gabriel Rucabado and Antonio Garcia Dopico>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

//!\file solverIm-gpu.cpp
//  \brief The explicit gpu solver of Neurite.


#include "solverEx-gpu.h"
#include "discretization.h"
#include "cudaException.h" 
#include "configuration.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;

// Includes from cuda
#include "cuda_runtime.h"



//#######################################
//#### HEADERS of PRIVATE FUNCTIONS #####
//#######################################
void printVector (char *name, int iter, TYPE *vector, int numElem);

//#######################################
//#### BODY of PRIVATE FUNCTIONS ########
//#######################################
void printVector (char *name, int iter, TYPE *vector, int numElem){
  TYPE *buff	= new TYPE [numElem];
  char value[4096];
  ofstream myfile;
  int i;

  CUDA_CPY_SYNC(buff, vector, sizeof(TYPE)*numElem);

  myfile.open (name, ios::out | ios::app);
  if (myfile.is_open()){
    sprintf (value, "Iter:: %05d ==>", iter);
    myfile << value;
    for (i=0; i<numElem; i++){
      sprintf (value, " [%4d]%.30e ", i, buff[i]);
      myfile << value;
    }
    myfile << endl;
  } else{
    _ERROR("ERROR OPENING FILE (%s)\n", name);
  }
  myfile.close();
}



//#######################################
//#### HEADERS of CUDA FUNCTIONS ########
//#######################################
__global__ void explicit_solver_normal (int num_elem, TYPE *new_Pot, TYPE *old_Pot,
					TYPE dT, TYPE *c_m, TYPE *dx, TYPE *r_a, 
					TYPE *W, TYPE *K, TYPE *input);
__global__ void explicit_solver_branch (int num_elem, TYPE *new_Pot, TYPE *old_Pot,
					int * ind_mo_bra, int *ind_branch, int *ind_dl,
					TYPE dT, TYPE *c_m, TYPE *dx, TYPE *r_a, 
					TYPE *W, TYPE *K, TYPE *input);
__global__ void explicit_solver_terminal (int num_elem, int *ind_term,
					  int * ind_mo_term, TYPE *new_Pot);
__global__ void Update_HH (int num_elem, TYPE *old_Pot, TYPE *rest_pot, int *ind_HH, 
			   TYPE *m, TYPE *h, TYPE *n, TYPE dT, 
			   TYPE *G_L, TYPE *E_Na, TYPE *E_K, TYPE *E_L, 
			   TYPE *g_Na, TYPE *g_K, TYPE *G_Na, TYPE *G_K, 
			   TYPE *leftShiftNa, TYPE *leftShiftK,
			   TYPE *W, TYPE *K);



//###########################################
//#### BODY of CUDA FUNTIONS ################
//###########################################
//! \brief Update electrical properties of HH nodes
__global__ void Update_HH (int num_elem, TYPE *old_Pot, TYPE *rest_pot, int *ind_HH, 
			   TYPE *m, TYPE *h, TYPE *n, TYPE dT, 
			   TYPE *G_L, TYPE *E_Na, TYPE *E_K, TYPE *E_L, 
			   TYPE *g_Na, TYPE *g_K, TYPE *G_Na, TYPE *G_K, 
			   TYPE *leftShiftNa, TYPE *leftShiftK,
			   TYPE *W, TYPE *K){
  TYPE dt	= dT * 1000;
  int id	= blockIdx.x * blockDim.x + threadIdx.x;
  if(id < num_elem){
    TYPE temp_G_Na, temp_G_K;
    TYPE mloc, nloc, hloc;
    TYPE Am, Ah, An, Bm, Bh, Bn; 
    int index	  = ind_HH[id];

    TYPE V_prevNa = ((old_Pot[index] + leftShiftNa[id]) - rest_pot[id]) / SCALE; 
    TYPE V_prevK  = ((old_Pot[index] + leftShiftK[id]) - rest_pot[id]) / SCALE; 

    mloc	= m[id];
    hloc	= h[id];
    nloc	= n[id];

    Am 	= (2.5 - 0.1 * V_prevNa) / (exp(2.5 - 0.1 * V_prevNa) - 1);
    Bm 	= 4 * exp(-V_prevNa / 18);
    Ah 	= 0.07 * exp(-V_prevNa / 20);
    Bh	= 1 / (exp(3 - 0.1 * V_prevNa) + 1);
    An	= (0.1 - 0.01 * V_prevK) / (exp(1 - 0.1 * V_prevK)-1);
    Bn	= 0.125 * exp(-V_prevK / 80);	

    m[id] = mloc = mloc + dt * (Am * (1 - mloc) - Bm * mloc);	    
    h[id] = hloc = hloc + dt * (Ah * (1 - hloc) - Bh * hloc);	  
    n[id] = nloc = nloc + dt * (An * (1 - nloc) - Bn * nloc);

    G_Na[id]  = temp_G_Na = g_Na[id] * mloc * mloc * mloc * hloc;
    G_K[id]   = temp_G_K  = g_K[id] * nloc * nloc * nloc * nloc;

    W[index] = -(temp_G_Na + temp_G_K + G_L[id]);
    K[index] = temp_G_Na * E_Na[id] + temp_G_K * E_K[id] + G_L[id] * E_L[id];
  }
}



//! \brief Figure out all nodes like if they were normal nodes for one iteration of explicit solver
__global__ void explicit_solver_normal (int num_elem, TYPE *new_Pot, TYPE *old_Pot,
					TYPE dT, TYPE *c_m, TYPE *dx, TYPE *r_a, 
					TYPE *W, TYPE *K, TYPE *input){
  register TYPE general, dr, mother, me;
  register int id = blockIdx.x * blockDim.x + threadIdx.x;

  if (id > 0 && id < num_elem-1){
    general	= dT / (c_m[id] * dx[id]);
    dr		= (old_Pot[id+1] - old_Pot[id]) / (r_a[id+1] * dx[id+1]);
    mother	= (old_Pot[id-1] - old_Pot[id]) / (r_a[id] * dx[id]);
    me		= W[id] * dx[id] * old_Pot[id] + K[id] * dx[id] + input[id];

    new_Pot[id] = old_Pot[id] + general * (dr + mother + me); 
  }
}



//! \brief Figure out only the branches nodes for one iteration of explicit solver
__global__ void explicit_solver_branch (int num_elem, TYPE *new_Pot, TYPE *old_Pot,
					int * ind_mo_bra, int *ind_branch, int *ind_dl,
					TYPE dT, TYPE *c_m, TYPE *dx, TYPE *r_a, 
					TYPE *W, TYPE *K, TYPE *input){
  register int index, index_dl, id;
  register TYPE general_B, dl_B, dr_B, mother_B, me_B;
  register TYPE general_D, dr_D, mother_D, me_D;
  
  id	 = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < num_elem){
    index	= ind_branch[id];
    index_dl	= abs(ind_dl[id]);
  
    general_B	= dT / (c_m[index] * dx[index]);
    dl_B	= (old_Pot[index_dl] - old_Pot[index])/ (r_a[index_dl] * dx[index_dl]);
    dr_B	= (old_Pot[index+1] - old_Pot[index]) / (r_a[index+1] * dx[index+1]);
    mother_B	= (old_Pot[ind_mo_bra[id]] - old_Pot[index])/ (r_a[index] * dx[index]);
    me_B	= W[index] * dx[index] * old_Pot[index] + K[index] * dx[index] + input[index];
    
    new_Pot[index]	= old_Pot[index] + general_B * (dl_B + dr_B + mother_B + me_B);

    // my LD is branch type
    if ( ind_dl[id] >0 ) {
      general_D	= dT / (c_m[index_dl] * dx[index_dl]);
      dr_D	= (old_Pot[index_dl+1] - old_Pot[index_dl]) / (r_a[index_dl+1] * dx[index_dl+1]);
      mother_D	= (old_Pot[index] - old_Pot[index_dl]) / (r_a[index_dl] * dx[index_dl]);
      me_D	= W[index_dl] * dx[index_dl] * old_Pot[index_dl] + K[index_dl] * dx[index_dl] + input[index_dl]; 
      new_Pot[index_dl] 	= old_Pot[index_dl] + general_D * (dr_D + mother_D + me_D); 
    }
  }
}



//! \brief Figure out only the terminal nodes for one iteration of explicit solver
__global__ void explicit_solver_terminal (int num_elem, int *ind_term, 
					  int *ind_mo_term, TYPE *new_Pot){
  int id	= blockIdx.x * blockDim.x + threadIdx.x;
  int index 	= ind_term[id];

  if (id == 0)
    new_Pot[index] = new_Pot[index + 1];
  else if (id < num_elem){
    new_Pot[index] = new_Pot[ind_mo_term[id]];
  }
}






//###################################################
//#### BODY of PUBLIC METHODS #######################
//###################################################
SolverExGPU::SolverExGPU(TYPE dT, TYPE **potential,
			 int num_term, int *ind_terminals, TYPE input_val,
			 int num_HH_elements, Hodgkin_Huxley **HH_elements,
			 int num_elem, discretization **elements) throw(){
  int i;
  cudaError_t error;
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
    DIM3(dimBlock_all, TH_ALL,1,1);
    DIM3(dimBlock_term, TH_TERM,1,1);
    DIM3(dimBlock_branch, TH_BRANCH,1,1);
    DIM3(dimBlock_HH,TH_HH,1,1);
    this->dT 	= dT;

    if ((error = cudaStreamCreate(&stream[0]))!= OK_CUDA){
      std::stringstream params;
      params << "stream=" << &stream[0];
      throw cudaException(__FILE__, __FUNCTION__, __LINE__,	       
			  "creating the stream 0", error,  params.str());
    }
    if ((error = cudaStreamCreate(&stream[1]))!= OK_CUDA){
      std::stringstream params;
      params << "stream=" << &stream[1];
      throw cudaException(__FILE__, __FUNCTION__, __LINE__,
			  "creating the 1", error,  params.str());
    }


    if ((error = cudaEventCreateWithFlags(&evento[0], cudaEventDisableTiming)) != OK_CUDA){
      std::stringstream params;
      params << "event=" << &evento[0];
      throw cudaException(__FILE__, __FUNCTION__, __LINE__,
			  "Creating event", error,  params.str());
    }
    if ((error = cudaEventCreateWithFlags(&evento[1], cudaEventDisableTiming)) != OK_CUDA){
      std::stringstream params;
      params << "event=" << &evento[1];
      throw cudaException(__FILE__, __FUNCTION__, __LINE__,
			  "Creating event", error,  params.str());
    }

    cudaFuncSetCacheConfig(explicit_solver_normal, cudaFuncCachePreferL1 );
    cudaFuncSetCacheConfig(explicit_solver_branch, cudaFuncCachePreferL1 );
    cudaFuncSetCacheConfig(explicit_solver_terminal, cudaFuncCachePreferL1 );
    cudaFuncSetCacheConfig(Update_HH, cudaFuncCachePreferL1 );

  
    // TERMINALS
    //////////////////////////////////////////////////
    this->num_term 	= num_term;
    int *CPUind_mo_term	= new int[this->num_term];
    for (i=0; i<num_term; i++){
      CPUind_mo_term[i]= elements[ind_terminals[i]]->get_mother();
    }

    DIM3(dimGrid_term, (this->num_term + TH_TERM - 1)/TH_TERM,1,1);
    int num_term_cuda = dimGrid_term.x * dimBlock_term.x;
    CUDA_ALLOC (&(this->ind_term), num_term_cuda * sizeof (int));
    CUDA_ALLOC (&(this->ind_mo_term), num_term_cuda * sizeof (int));

    CUDA_CPY   (this->ind_term, ind_terminals, sizeof(int) * num_term, stream[0]);
    CUDA_CPY   (this->ind_mo_term, CPUind_mo_term, sizeof(int) * num_term, stream[0]);



    // BRANCHING
    //////////////////////////////////////////////////
    _DEBUG("Before creating structure temp of BRANCH");
    this->num_branch 	= num_term -2;    
    int *CPUind_branch	= new int[this->num_branch];
    int *CPUind_mo_bran	= new int[this->num_branch];
    int *CPUind_DL	= new int[this->num_branch];
    for (i=0, this->num_branch=0; i<num_elem; i++){
      if (elements[i]->get_branching()){
	CPUind_branch[this->num_branch]	= i;
	CPUind_mo_bran[this->num_branch]= elements[i]->get_mother();
	if (elements[elements[i]->get_DL()]->get_branching()){
	  CPUind_DL[this->num_branch]	= -elements[i]->get_DL();
	}else{
	  CPUind_DL[this->num_branch]	= +elements[i]->get_DL();
	}
	this->num_branch++;
      }
    }

    DIM3(dimGrid_branch, (this->num_branch + TH_BRANCH - 1)/TH_BRANCH,1,1);
    int num_branch_cuda = dimGrid_branch.x * dimBlock_branch.x;

    CUDA_ALLOC(&(this->ind_branch), num_branch_cuda * sizeof(int));
    CUDA_ALLOC(&(this->ind_mo_branch), num_branch_cuda * sizeof(int));
    CUDA_ALLOC(&(this->ind_DL), num_branch_cuda * sizeof(int));

    CUDA_CPY(this->ind_branch, CPUind_branch, this->num_branch * sizeof(int), stream[0]);
    CUDA_CPY(this->ind_mo_branch, CPUind_mo_bran, this->num_branch * sizeof(int), stream[0]);
    CUDA_CPY(this->ind_DL, CPUind_DL, this->num_branch * sizeof(int), stream[1]);

    // ALL
    //////////////////////////////////////////////////
    this->num_elem	= num_elem;
    CUDA_ALLOC(&(this->potential[0]), num_elem * sizeof(TYPE));
    CUDA_ALLOC(&(this->potential[1]), num_elem * sizeof(TYPE));
    CUDA_CPY(this->potential[0], potential[0], sizeof(TYPE)*num_elem, stream[0]);
    CUDA_CPY(this->potential[1], potential[1], sizeof(TYPE)*num_elem, stream[1]);
    
    this->new_Pot		= this->potential[1];
    this->old_Pot		= this->potential[0];

    //   ATTRIBUTES
    TYPE *CPU_dx	= new TYPE[num_elem];
    TYPE *CPU_c_m	= new TYPE[num_elem];
    TYPE *CPU_r_a	= new TYPE[num_elem];
    TYPE *CPU_W		= new TYPE[num_elem];
    TYPE *CPU_K		= new TYPE[num_elem];
    TYPE *CPU_input_on	= new TYPE[num_elem];
    TYPE *CPU_input_off	= new TYPE[num_elem];
    for (i=0; i<num_elem; i++){
      CPU_dx[i]		= elements[i]->get_dX();
      CPU_c_m[i]	= elements[i]->get_c_m();
      CPU_r_a[i]	= elements[i]->get_r_a();
      CPU_W[i]		= elements[i]->get_W();
      CPU_K[i]		= elements[i]->get_K();
      CPU_input_off[i]	= 0;
      CPU_input_on[i]	= elements[i]->get_input_current();
    }

    DIM3(dimGrid_all, (this->num_elem + TH_ALL - 1)/TH_ALL,1,1);
    int num_all_cuda = dimGrid_all.x * dimBlock_all.x;

    _DEBUG("Before creating structures of ALL(%d-%d)", 
	   this->num_elem, num_all_cuda);
    CUDA_ALLOC(&(this->dx), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->c_m), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->r_a), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->W), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->K), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->input_on), num_all_cuda * sizeof(TYPE));
    CUDA_ALLOC(&(this->input_off), num_all_cuda * sizeof(TYPE));
    CUDA_CPY(this->dx, CPU_dx, this->num_elem * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->c_m, CPU_c_m, this->num_elem * sizeof(TYPE), stream[1]);
    CUDA_CPY(this->r_a, CPU_r_a, this->num_elem * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->W, CPU_W, this->num_elem * sizeof(TYPE), stream[1]);
    CUDA_CPY(this->K, CPU_K, this->num_elem * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->input_on, CPU_input_on, this->num_elem * sizeof(TYPE), stream[1]);
    CUDA_CPY(this->input_off, CPU_input_off, this->num_elem * sizeof(TYPE), stream[0]);

    _DEBUG("Before creating structures of HH");
    // HH
    //////////////////////////////////////////////////
    this->num_HH	= num_HH_elements;
    int *CPU_ind_HH	= new int [num_HH_elements];
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
      CPU_ind_HH[i]= HH_elements[i]->get_element();
      CPU_Shift_Na[i]	= HH_elements[i]->get_Left_Shift_Na();
      CPU_Shift_K[i]	= HH_elements[i]->get_Left_Shift_K();
      CPU_n[i]	= HH_elements[i]->get_n();
      CPU_m[i]	= HH_elements[i]->get_m();
      CPU_h[i]	= HH_elements[i]->get_h();
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

    CUDA_ALLOC(&(this->ind_HH), num_HH_cuda * sizeof(int));
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

    CUDA_CPY(this->ind_HH, CPU_ind_HH, this->num_HH * sizeof(int), stream[0]);
    CUDA_CPY(this->leftShiftNa, CPU_Shift_Na, this->num_HH * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->leftShiftK, CPU_Shift_K, this->num_HH * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->n, CPU_n, this->num_HH * sizeof(TYPE), stream[1]);
    CUDA_CPY(this->m, CPU_m, this->num_HH * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->h, CPU_h, this->num_HH * sizeof(TYPE), stream[1]);
    CUDA_CPY(this->E_Na, CPU_E_Na, this->num_HH * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->G_Na, CPU_G_Na, this->num_HH * sizeof(TYPE), stream[1]);
    CUDA_CPY(this->g_Na, CPU_g_Na, this->num_HH * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->E_L, CPU_E_L, this->num_HH * sizeof(TYPE), stream[1]);
    CUDA_CPY(this->G_L, CPU_G_L, this->num_HH * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->E_K, CPU_E_K, this->num_HH * sizeof(TYPE), stream[1]);
    CUDA_CPY(this->G_K, CPU_G_K, this->num_HH * sizeof(TYPE), stream[0]);
    CUDA_CPY(this->g_K, CPU_g_K, this->num_HH * sizeof(TYPE), stream[1]);
    CUDA_CPY(this->rest_pot, CPU_rest_pot, this->num_HH * sizeof(TYPE), stream[0]);
  }catch (cudaException &e){
    e.print();
    throw e;
  }
}

  //!\brief This function updates the electrical variables based on the previous time step information
  //  @param input_active Flag used to indicate if it has to use  the input current to update the values
void SolverExGPU::update (bool input_active) throw(){
  cudaError_t error;

  try{
    // Set GPU card to use
    error = cudaSetDevice (this->idGPU);

    // Set INPUTS to extrems (only first one)
    if (input_active){
      input	= input_on;
    }else{
      input	= input_off;
    }

    if (num_HH > 0){  
      if ((error = cudaStreamWaitEvent (stream[0], evento[1], 0)) != OK_CUDA){
	std::stringstream params;
	params << "event=" << &evento[1] << " stream="<< stream[0];
	throw cudaException(__FILE__, __FUNCTION__, __LINE__,
			    "Synchronizing Update_HH", error,  params.str());
      }
      KERNEL(Update_HH, dimGrid_HH, dimBlock_HH, stream[0], 
     	     "%d, %p, %p, %p, %p, %p, %p, %f, %p, %p, %p,"
	     "%p, %p, %p, %p, %p, %p, %p, %p, %p",
     	     num_HH, old_Pot, rest_pot, ind_HH, m, h, n, dT, 
     	     G_L, E_Na, E_K, E_L, g_Na, g_K, G_Na, G_K,
	     leftShiftNa, leftShiftK, W, K);
      
    }
  }catch (cudaException &e){
    e.print();
    throw e;
  }
}



//! \brief Calculate one iteration of explicit finite difference discretization solver
void SolverExGPU::calculate() throw(){
  cudaError_t error;

  try{
    // Set GPU card to use
    error = cudaSetDevice (this->idGPU);

    KERNEL(explicit_solver_normal, dimGrid_all, dimBlock_all, stream[0], 
	   "%d, %p, %p, %f, %p, %p, %p, %p, %p, %p",
	   num_elem, new_Pot, old_Pot, dT, c_m, dx, r_a, W, K, input);

    if ((error = cudaEventRecord(evento[0], stream[0])) != OK_CUDA){
      std::stringstream params;
      params << "event=" << &evento << " stream="<< &stream[0];
      throw cudaException(__FILE__, __FUNCTION__, __LINE__,
			  "Tracking event", error,  params.str());
    }
    if ((error = cudaStreamWaitEvent (NULL, evento[0], 0)) != OK_CUDA){
      std::stringstream params;
      params << "event=" << &evento << " stream=NULL";
      throw cudaException(__FILE__, __FUNCTION__, __LINE__,
			  "Synchronizing the solver", error,  params.str());
    }

    if (num_branch > 0){
           KERNEL(explicit_solver_branch, dimGrid_branch, dimBlock_branch, stream[0], 
	     "%p, %d, %p, %p, %p, %p, %f, %p, %p, %p, %p, %p, %p",
	     num_branch, new_Pot, old_Pot, ind_mo_branch, ind_branch, ind_DL, dT, c_m, dx, r_a, W, K, input);
    }
    KERNEL(explicit_solver_terminal, dimGrid_term, dimBlock_term, stream[0], 
	   "%d, %p, %p, %p", 
	   num_term, ind_term, ind_mo_term, new_Pot);

    if ((error = cudaEventRecord(evento[1], NULL)) != OK_CUDA){
      std::stringstream params;
      params << "event=" << &evento[1] << " stream=NULL";
      throw cudaException(__FILE__, __FUNCTION__, __LINE__,
			  "Tracking the event", error,  params.str());
    }
  }catch (cudaException &e){
    e.print();
    throw e;
  }


  // flip act potential with new potential
  TYPE *flip 	= this->new_Pot;
  this->new_Pot	= this->old_Pot;
  this->old_Pot	= flip;
}

  //! \brief Copy the current potential to the pointer of the parameter
  //  @param potential Pointer on which the potential will be copied 
void SolverExGPU::snapshot(TYPE *potential) throw(){
  try{
    CUDA_CPY(potential, this->old_Pot, sizeof(TYPE) * this->num_elem, stream[0]);
  }catch (cudaException &e){
    e.print();
    throw e;
  }
}



//! \brief Free all memory used by solver
SolverExGPU::~SolverExGPU(){
  try{
    CUDA_FREE(this->ind_term);
    CUDA_FREE(this->ind_mo_term);
  
    CUDA_FREE(this->ind_DL);
    CUDA_FREE(this->ind_branch);
    CUDA_FREE(this->ind_mo_branch);

    CUDA_FREE(this->potential[0]);
    CUDA_FREE(this->potential[1]);
    CUDA_FREE(this->dx);
    CUDA_FREE(this->c_m);
    CUDA_FREE(this->r_a);
    CUDA_FREE(this->W);
    CUDA_FREE(this->K);
    CUDA_FREE(this->input_on);
    CUDA_FREE(this->input_off);

    CUDA_FREE(this->ind_HH);
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
  }catch (cudaException &e){
    e.print();
    throw e;
  }
}
  
