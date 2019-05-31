//
//
// File author(s):  <Julian Andres Garcia Grajales and Gabriel Rucabado>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file solverEx-cpu.cpp
  \brief The Explicit cpu solver of Neurite
*/
#include "coupling.h"
#include "solverEx-cpu.h"
#include "discretization.h"

//#################################################
//#### BODY of PUBLIC METHODS #####################
//#################################################
SolverExCPU::SolverExCPU(TYPE dT, TYPE **potential, 
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

  this->input_on= new TYPE[num_elem];
  this->input_off= new TYPE[num_elem];
  for (int i=0; i<num_elem; i++){
    this->input_off[i]= 0;
    this->input_on[i]= element[i]->get_input_current();
  }
}

  //!\brief This function updates the electrical variables based on the previous time step information
  //  @param input_active Flag used to indicate if it has to use  the input current to update the values
void SolverExCPU::update (bool input_active) throw(){
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
    hhRate_Const(alpha, beta, HH->get_rest_pot(), this->prevT[HH->get_element()],
	   HH_elements[i]->get_Left_Shift_Na(), HH_elements[i]->get_Left_Shift_K());

    m = HH->get_m() + dt * (alpha.m*(1- HH->get_m()) - beta.m*HH->get_m());
    h = HH->get_h() + dt * (alpha.h*(1- HH->get_h()) - beta.h*HH->get_h());
    n = HH->get_n() + dt * (alpha.n*(1- HH->get_n()) - beta.n*HH->get_n());

    temp_G_Na = HH->get_g_Na()*m*m*m*h;
    temp_G_K = HH->get_g_K()*n*n*n*n;
    HH->set_G_Na(temp_G_Na);
    HH->set_G_K(temp_G_K);
    HH->set_W();
    HH->set_K();

    // Saving for the next time step
    HH->set_m(m);
    HH->set_h(h);
    HH->set_n(n);
  }
}


 //!\brief Calculate one iteration of the explicit finite difference discretization solver
void SolverExCPU::calculate() throw(){
  register TYPE RprevT;
  register discretization *Relement;
  // ALL elements

  for(int j=1; j<numStepsX-1; j++){
    RprevT 	= this->prevT[j];
    Relement 	= this->element[j];
    
    this->currT[j] = RprevT + (this->dT/(Relement->get_c_m()*Relement->get_dX()))
      *(Relement->get_branching()
	*((this->prevT[Relement->get_DL()]-RprevT)
	  /(this->element[Relement->get_DL()]->get_r_a()
	    *this->element[Relement->get_DL()]->get_dX())) 
	+ (this->prevT[Relement->get_DR()]-RprevT)
	/(this->element[Relement->get_DR()]->get_r_a()
	  *this->element[Relement->get_DR()]->get_dX())
	+ (this->prevT[Relement->get_mother()]-RprevT)
	/(Relement->get_r_a()*Relement->get_dX())
	+ Relement->get_W()*Relement->get_dX()*RprevT + Relement->get_K()*Relement->get_dX()
	+ this->input[j]);
  }

  // TERMINALS
  // First element is different  // First fictitius element
  this->currT[this->terminals[0]] = this->currT[this->terminals[0]+1];
  // For the other terminals.
  for(int k=1; k<this->num_terminals; k++){
    this->currT[this->terminals[k]]=this->currT[this->element[this->terminals[k]]->get_mother()];
  }

  //Flip prev and current potential vectors
  TYPE *flip 	= this->currT;
  this->currT	= this->prevT;
  this->prevT	= flip;

}

  //! \brief Copy current potential to the pointer of the parameter
  //  @param potential Pointer on which the potential will be copied 
void SolverExCPU::snapshot(TYPE *potential) throw(){
  memcpy (potential, this->prevT, sizeof(TYPE)*this->numStepsX);
}

//! \brief Free all memory used by solver
SolverExCPU::~SolverExCPU(){
  delete [] this->input_on;
  delete [] this->input_off;
}
