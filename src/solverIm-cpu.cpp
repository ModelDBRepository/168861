//
//
// File author(s):  <Julian Andres Garcia Grajales and Gabriel Rucabado>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

//!\file solverIm-cpu.cpp
//  \brief The implicit cpu solver of Neurite
#include "coupling.h"
#include "solverIm-cpu.h"
#include "discretization.h"
#include "mkl_rci.h"
#include "mkl_blas.h"
#include <mkl_spblas.h>
#include "mkl_service.h"

#include <iostream>
#include <fstream>
using namespace std;



//#################################################
//#### BODY of PRIVATE METHODS ####################
//#################################################
//! \brief This function create the main matrix
void SolverImCPU::createMatrix(){
  int j;
  TYPE alpha=0, beta=0, gamma=0, delta=0;
  register discretization *Relement, *NextRelement;

  this->nnz=1;
  bool flag_salto= false;
  int max_size =0;
 
  // I know how many zeros I have in my matrix = this->numStepsX-2(rows with at most 4 elements)* 
  // 4 (elements = mother, me, DR and DL) + 1 (fictitius element) + 1 (fictitius element)
  max_size =this->numStepsX*4 +1 +1;

  this->csr_Col_Ind = (MKL_INT *)calloc(max_size,sizeof(MKL_INT));
  this->csr_Val = (TYPE *)calloc(max_size,sizeof(TYPE));
  this->csr_Row_Ptr = (MKL_INT *)calloc(this->numStepsX+1,sizeof(MKL_INT));

  this->csr_Val[0] = 1;
  this->csr_Col_Ind[0] = 0;
  this->csr_Row_Ptr[0] = 0;

  Relement = this->element[1];
  gamma = -Relement->get_branching()/(this->element[Relement->get_DL()]->get_r_a()*this->element[Relement->get_DL()]->get_dX());
  beta = (Relement->get_c_m()*Relement->get_dX())/this->dT - Relement->get_W()*Relement->get_dX() + Relement->get_branching()/(this->element[Relement->get_DL()]->get_r_a()*this->element[Relement->get_DL()]->get_dX()) + 1/(this->element[Relement->get_DR()]->get_r_a()*this->element[Relement->get_DR()]->get_dX()) + 1/(Relement->get_r_a()*Relement->get_dX());
  alpha = -1/(Relement->get_r_a()*Relement->get_dX());
  delta = -1/(this->element[Relement->get_DR()]->get_r_a()*this->element[Relement->get_DR()]->get_dX());

  // BOUNDARY CONDITIONS
  beta = beta + alpha;
  this->csr_Val[this->nnz] = beta;
  this->csr_Col_Ind[this->nnz] = 1;
  this->csr_Row_Ptr[1] = this->nnz;
  this->nnz++;
  if(Relement->get_DL() == 0){
    this->csr_Val[this->nnz] = delta;
    this->csr_Col_Ind[this->nnz] = Relement->get_DR();
    this->nnz++;
  }else{   
    if(Relement->get_DR()<Relement->get_DL()){
      this->csr_Val[this->nnz] = delta;
      this->csr_Col_Ind[this->nnz] = Relement->get_DR();
      this->nnz++;
      this->csr_Val[this->nnz] = gamma;
      this->csr_Col_Ind[this->nnz] = Relement->get_DL();
      this->nnz++;
    } else{
      this->csr_Val[this->nnz] = gamma;
      this->csr_Col_Ind[this->nnz] = Relement->get_DL();
      this->nnz++;
      this->csr_Val[this->nnz] = delta;
      this->csr_Col_Ind[this->nnz] = Relement->get_DR();
      this->nnz++;
    }
  }
  this->currT[1] = Relement->get_c_m()*Relement->get_dX()*this->prevT[1]/this->dT + Relement->get_K()*Relement->get_dX() + this->input[1];
 
  // Whole neurite
  int jedit=2;
  for(j=2; j<this->numStepsX-2; j++)
    {
      Relement = this->element[j];
      NextRelement = this->element[j+1];

      if(flag_salto){
	flag_salto = false;
	this->csr_Val[this->nnz] = 1;
	this->csr_Col_Ind[this->nnz] = j;
	this->csr_Row_Ptr[jedit] = this->nnz;
	this->nnz++;
	jedit++;
      }else{
	if(NextRelement->get_DR() == j+1){
	  flag_salto = true;
	  gamma = -Relement->get_branching()/(this->element[Relement->get_DL()]->get_r_a()*this->element[Relement->get_DL()]->get_dX());
	  beta = (Relement->get_c_m()*Relement->get_dX())/this->dT - Relement->get_W()*Relement->get_dX() + Relement->get_branching()/(this->element[Relement->get_DL()]->get_r_a()*this->element[Relement->get_DL()]->get_dX()) + 1/(this->element[Relement->get_DR()]->get_r_a()*this->element[Relement->get_DR()]->get_dX()) + 1/(Relement->get_r_a()*Relement->get_dX());
	  alpha = -1/(Relement->get_r_a()*Relement->get_dX());
	  delta = -1/(this->element[Relement->get_DR()]->get_r_a()*this->element[Relement->get_DR()]->get_dX());

	  this->csr_Val[this->nnz] = alpha;
	  this->csr_Col_Ind[this->nnz] = Relement->get_mother();
	  this->csr_Row_Ptr[jedit] = this->nnz;
	  this->nnz++;
	  this->csr_Val[this->nnz] = beta+delta;
	  this->csr_Col_Ind[this->nnz] = j;
	  this->nnz++;

	  this->currT[j] = Relement->get_c_m()*Relement->get_dX()*this->prevT[j]/this->dT + Relement->get_K()*Relement->get_dX() + this->input[j];
	}else{
	  gamma = -Relement->get_branching()/(this->element[Relement->get_DL()]->get_r_a()*this->element[Relement->get_DL()]->get_dX());
	  beta = (Relement->get_c_m()*Relement->get_dX())/this->dT - Relement->get_W()*Relement->get_dX() + Relement->get_branching()/(this->element[Relement->get_DL()]->get_r_a()*this->element[Relement->get_DL()]->get_dX()) + 1/(this->element[Relement->get_DR()]->get_r_a()*this->element[Relement->get_DR()]->get_dX()) + 1/(Relement->get_r_a()*Relement->get_dX());
	  alpha = -1/(Relement->get_r_a()*Relement->get_dX());
	  delta = -1/(this->element[Relement->get_DR()]->get_r_a()*this->element[Relement->get_DR()]->get_dX());
	  this->csr_Val[this->nnz] = alpha;
	  this->csr_Col_Ind[this->nnz] = Relement->get_mother();
	  this->csr_Row_Ptr[jedit] = this->nnz;
	  this->nnz++;
	  this->csr_Val[this->nnz] = beta;
	  this->csr_Col_Ind[this->nnz] = j;
	  this->nnz++;
	  // HACKKKKKKKKK if you are a terminal your dougther right is yourself. The neurite never ends in a branching
	  if(Relement->get_DL() == 0){
	    this->csr_Val[this->nnz] = delta; // because gamma is 0
	    this->csr_Col_Ind[this->nnz] = Relement->get_DR();
	    this->nnz++;
	  }else{
	    this->csr_Val[this->nnz] = delta;
	    this->csr_Col_Ind[this->nnz] = Relement->get_DR();
	    this->nnz++;
	    this->csr_Val[this->nnz] = gamma;
	    this->csr_Col_Ind[this->nnz] = Relement->get_DL();
	    this->nnz++;
	  }
	  this->currT[j] = Relement->get_c_m()*Relement->get_dX()*this->prevT[j]/this->dT + Relement->get_K()*Relement->get_dX() + this->input[j];
	}
      	jedit++;
      }
    }

  //  Boundary conditions:
  Relement = this->element[this->numStepsX-2];
  gamma = -Relement->get_branching()/(this->element[Relement->get_DL()]->get_r_a()*this->element[Relement->get_DL()]->get_dX());
  beta = (Relement->get_c_m()*Relement->get_dX())/this->dT - Relement->get_W()*Relement->get_dX() + Relement->get_branching()/(this->element[Relement->get_DL()]->get_r_a()*this->element[Relement->get_DL()]->get_dX()) + 1/(this->element[Relement->get_DR()]->get_r_a()*this->element[Relement->get_DR()]->get_dX()) + 1/(Relement->get_r_a()*Relement->get_dX());
  alpha = -1/(Relement->get_r_a()*Relement->get_dX());
  delta = -1/(this->element[Relement->get_DR()]->get_r_a()*this->element[Relement->get_DR()]->get_dX());


  // BOUNDARY CONDITIONS
  beta = beta + delta;
  this->csr_Val[this->nnz] = alpha;
  this->csr_Col_Ind[this->nnz] = Relement->get_mother();
  this->csr_Row_Ptr[jedit] = this->nnz;
  jedit++;
  this->nnz++;
  this->csr_Val[this->nnz] = beta;
  this->csr_Col_Ind[this->nnz] = j;
  this->nnz++;
  if(Relement->get_DL() != 0){
    this->csr_Val[this->nnz] = gamma;
    this->csr_Col_Ind[this->nnz] = Relement->get_DL();
    this->nnz++;	
  }
  this->currT[this->numStepsX-2] = Relement->get_c_m()*Relement->get_dX()*this->prevT[this->numStepsX-2]/this->dT + Relement->get_K()*Relement->get_dX() + this->input[this->numStepsX-2];

  this->csr_Val[this->nnz] = 1;
  this->csr_Col_Ind[this->nnz] = this->numStepsX-1;
  this->csr_Row_Ptr[jedit] = this->nnz;
  jedit++;
  this->nnz++;
  this->csr_Row_Ptr[jedit]= this->nnz;


  // In order to re-size these vectors once we know how many non-zeros we have in out matrix
  this->csr_Col_Ind = (MKL_INT  *)realloc(this->csr_Col_Ind,sizeof(MKL_INT)*this->csr_Row_Ptr[jedit]);
  this->csr_Val     = (TYPE *)  realloc(this->csr_Val,    sizeof(TYPE) *this->csr_Row_Ptr[jedit]);
}




/*!  Flexible generalized minimal residual method for solving the linear system of equations (Ax=b) resultant of the implicit scheme for the electrical simulation

  For information about this function, please go to the mkl intel libraries documentation
*/
void SolverImCPU::fgmres(int N, TYPE *rhs, TYPE *yoacsr, MKL_INT *yoja, MKL_INT *yoia) {
  // Please check the mkl documentation
  MKL_INT size=128;
  MKL_INT ipar[size];
  TYPE dpar[size], *tmp;
  // Initial guess
  TYPE *computed_solution;
 
  computed_solution = (TYPE *)calloc(N,sizeof(TYPE));

  MKL_INT itercount; 
  MKL_INT RCI_request, ivar; 
  char cvar;

  TYPE tol_exit=RANGE_ERROR; // This tolerance because my criteria for exiting the solver is based on the euclidean norm of the residual vector.

  ivar=N; 
  cvar='N'; 

  // This part has changed respect to the original mkl example. There are bugs inside that example. Be carefull.
  ipar[14] = 150;
  if(N<ipar[14]){
    ipar[14] = N;
  }
      
  tmp = (TYPE *)calloc(N*(2*ipar[14]+1)+(ipar[14]*(ipar[14]+9))/2+1,sizeof(TYPE));
  dfgmres_init(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);

  // Initializing the solver values
  if (RCI_request!=0) goto FAILED;
  // Modifing in our particular simulations
  // ipar[0] = N*5;
  ipar[4]=N; // Editing the maximum number of iterations. default min(150,Nx)
  ipar[7]=1; // 1 = Check the maximum number of iterarions
  ipar[8]=0; // 1 =Residual stopping test. This option does not work in my particular case, because the exit tolerance = tolerance*initial_residual_eucli_norm + absolute tolerance. In my case initial_residual_eucli_norm is almost 0 (8E-10). I defined my own criteria. I test if the current euclidean norm of the residual is close to zero. This criteria has been validated against the gauss results.
  ipar[9]=1; // 1 = User defined stopping test
  ipar[11]=1; // 1 = test for zero norm of the solution
  ipar[12]=1; // Put the final result in rhs instead of in computing solution. Be carefull with this value.
 
  // Checking what it changed
  dfgmres_check(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);

  if (RCI_request!=0) goto FAILED; 

  // If not, we solve the system
 ONE: 
  dfgmres(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);

  if (RCI_request==0) goto COMPLETE; // With my criteria this way is impossible, because we test the tolerances manually
  
  if (RCI_request==1) {
    // A new iteration is necessary
    //mkl_cspblas_dcsrgemv(&cvar, &ivar, acsr, ia, ja, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
    mkl_cspblas_dcsrgemv(&cvar, &ivar, yoacsr, yoia, yoja, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
    
    if(ipar[3]==ipar[4]){ // Checking the maximum number of iterations
      _INFO("Error you exceed the maximum number of iterations %i in FGMRES\n Try decreasing the time step\n The current residual norm is %g but the tolerance criterium is %g",ipar[4],dpar[4], tol_exit);
      exit(0);
    }
    goto ONE; // returns the control to dfgmres
  } 
  if(RCI_request==2){ // Cheking the stopping test manually
    if(dpar[4]>tol_exit){
      goto ONE; // returns the control to dfgmres
    }else{
      goto COMPLETE; // if not, the system has been solved.
    }
  }else { 
    goto FAILED; 
  } 

 COMPLETE:
  dfgmres_get(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
  free(computed_solution);
  free(tmp);
  MKL_Free_Buffers();
  return;

 FAILED: 
  printf("This example FAILED as the solver has returned the ERROR code %d", RCI_request); exit(0); 
}


//#################################################
//#### BODY of PUBLIC METHODS #####################
//#################################################
SolverImCPU::SolverImCPU(TYPE dT, TYPE **potential, 
			   int num_term, int *ind_terminals, TYPE input_val,
			   int num_HH_elements, Hodgkin_Huxley **HH_elements,
			   int num_elem, discretization **element){
  this->dT	= dT;
  this->currT	= new TYPE [num_elem];
  memcpy (this->currT, potential[1], sizeof(TYPE)*num_elem);
  this->prevT	= new TYPE [num_elem];
  memcpy (this->prevT, potential[0], sizeof(TYPE)*num_elem);
  this->num_terminals 	= num_term;
  this->terminals	= ind_terminals;
  this->num_HH_elements	= num_HH_elements;
  this->HH_elements	= HH_elements;
  this->numStepsX	= num_elem;
  this->element		= element;
  this->csr_Col_Ind = NULL;
  this->csr_Row_Ptr = NULL;
  this->csr_Val = NULL;

  this->input_on= new TYPE[num_elem];
  this->input_off= new TYPE[num_elem];
  for (int i=0; i<num_elem; i++){
    this->input_off[i]= 0;
    this->input_on[i]= element[i]->get_input_current();
  }
}

  //!\brief This function updates the electrical variables based on the previous time step information
  //  @param input_active Flag used to indicate if it has to use  the input current to update the values
void SolverImCPU ::update (bool input_active) throw(){
  int i;
  rate_const alpha, beta;
  TYPE temp_G_Na, temp_G_K;
  TYPE m=0,n=0,h=0;
  Hodgkin_Huxley *HH;
  const  TYPE dt=this->dT*1000;
  // alpha.i, beta.i from hhprob are in 1/msec, we put the time, that was in sec, to msec. HH original paper 


  // Set INPUTS 
  if (input_active){
    this->input = this->input_on;
  } else {
    this->input = this->input_off;
  }


  // Updating HH elements
  for(i=0; i<this->num_HH_elements; i++){
    HH = this->HH_elements[i];
    // This is to take into account the mechanical model
    hhRate_Const(alpha, beta, HH->get_rest_pot(), this->prevT[HH->get_element()],
	   HH_elements[i]->get_Left_Shift_Na(), HH_elements[i]->get_Left_Shift_K());

    m = HH->get_m() + dt * (alpha.m*(1- HH->get_m()) - beta.m*HH->get_m());
    h = HH->get_h() + dt * (alpha.h*(1- HH->get_h()) - beta.h*HH->get_h());
    n = HH->get_n() + dt * (alpha.n*(1- HH->get_n()) - beta.n*HH->get_n());

    temp_G_Na = HH->get_g_Na()*m*m*m*h;
    temp_G_K = HH->get_g_K()*n*n*n*n;
    HH->set_G_Na(temp_G_Na);
    HH->set_G_K(temp_G_K);

    // Saving for the next time step
    HH->set_m(m);
    HH->set_h(h);
    HH->set_n(n);

    HH->set_W();
    HH->set_K();

  }


  // Implicit matrix creation
  if (this->csr_Col_Ind != NULL){
    free(this->csr_Col_Ind);
    free(this->csr_Row_Ptr);
    free(this->csr_Val);
  }
  this->createMatrix();
}


//! \brief Calculate one iteration of explicit finite difference discretization solver
void SolverImCPU ::calculate() throw(){
  this->fgmres(this->numStepsX, this->currT, this->csr_Val, this->csr_Col_Ind,this->csr_Row_Ptr);
  
  //Flip prev and actual potential vectors
  TYPE *flip 	= this->currT;
  this->currT	= this->prevT;
  this->prevT	= flip;
}

  //! \brief Copy the current potential to the pointer of the parameter
  //  @param potential Pointer on which the potential will be copied 
void SolverImCPU::snapshot(TYPE *potential) throw(){
  memcpy (potential, this->prevT, sizeof(TYPE)*this->numStepsX);
}


//! \brief Free all memory used by solver
SolverImCPU::~SolverImCPU(){
  delete [] this->input_on;
  delete [] this->input_off;
}

