//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file coupling.cpp
  \brief This file contains the functions related to the electrical properties that change during the time.
*/

#include "coupling.h"
extern int numStepsX;
extern TYPE tspan;
extern TYPE time_run;

/*! \brief This function calculates the probabilities of the Hodgkin and Huxley model
  @param V_prev Previous potential
  @param alpha Structure with alpha probabilities
  @param beta Structure with beta probabilities
  @param rest_pot Resting potential

*/
extern void hhRate_Const(rate_const &alpha, rate_const &beta, TYPE rest_pot,TYPE potential,TYPE Left_Shift_Na, TYPE Left_Shift_K)  
{

  TYPE V_prevNa=0,V_prevK=0;
  // Original from HH paper
  TYPE scale = 0.001; // In order to have the potential in mV for the HH equations

  // Applying the LS
  V_prevNa = potential + Left_Shift_Na;
  V_prevK = potential + Left_Shift_K;

  V_prevNa = (V_prevNa-rest_pot)/scale; // particular for HH dimensions
  V_prevK = (V_prevK-rest_pot)/scale; // particular for HH dimensions
  // Na gate values
  alpha.m = (2.5-0.1*V_prevNa)/(exp(2.5 - 0.1*V_prevNa)-1);
  beta.m = 4*exp(-V_prevNa/18);

  alpha.h = 0.07*exp(-V_prevNa/20);
  beta.h = 1/(exp(3 - 0.1*V_prevNa) + 1);

  // K gate value
  alpha.n=(0.1 - 0.01*V_prevK)/(exp(1 - 0.1*V_prevK)-1);
  beta.n=0.125*exp(-V_prevK/80);
}



