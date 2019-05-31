//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file discretization.cpp
  \brief In this file the functions related to the space discretization are defined
*/

#include "discretization.h"
#include "configuration.h"
#include "neurite.h"

/*!Default constructor for discretization class
*/
discretization::discretization(){
  // Default value
  element=0;
  epsilon=0;
  epsilon_surface =0;
  DL = 0; 
  DR = 0;
  mother = 0;
  HH = false;
  branching = false;
  diameter = 3*pow(10,-6); 
  thickness = 4*pow(10,-9);
  W = 0;
  K = 0;
  input_current = 0;
  ro_a = 1.87;
  ro_m = 2.5E9; 
  per_m = 4E-11;
  rest_pot = -0.0655;
  critical_dT=0;
}

/*! Constructor for Cable class.

  @param ele Number of element
  @param dr Daughter right of the element
  @param mom Mother of the element
  @param deltaX Element size

*/
Cable::Cable(int ele,int dr, int mom, TYPE deltaX){

  diameter_my = 5.56*pow(10,-6);  
  thickness_my = 18E-9;
  per_my = 1.08E-10;
  ro_my = 4.44E6; 
  layers_my = 45; 
  element = ele;
  DR = dr;
  mother = mom;
  dX = deltaX;

}

void Cable::set_W(){
  W = - 1/r_m;
}
void Cable::set_K(){
  K = rest_pot/r_m;
}

/*! Constructor for Hodgkin and Huxley class.

  @param ele Number of element
  @param dr Daughter right of the element
  @param mom Mother of the element
  @param deltaX Element size

*/
Hodgkin_Huxley::Hodgkin_Huxley(int ele,int dr, int mom, TYPE deltaX){
  HH = true;
  E_Na = 0.115 + rest_pot;
  E_Na0 = E_Na;
  EMAX = 0;
  E_K = -0.012 + rest_pot;
  E_K0=E_K;
  E_L=0;
  E_L0= E_L;
  sigma_Na = 4.8E-6;
  sigma_K = 1.44E-6;
  sigma_L= 1.2E-8; 
  n_damage = 2;
  element = ele;
  DR = dr;
  mother = mom;
  dX = deltaX;
  Left_Shift_Na=0;
  Left_Shift_K=0;
}

void Hodgkin_Huxley::set_W(){
  W = -(G_Na + G_K + G_L);
}
void Hodgkin_Huxley::set_K(){
  K = G_Na*E_Na + G_K*E_K + G_L*E_L;
}
