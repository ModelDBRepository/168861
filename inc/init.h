//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


#ifndef _init_H_
#define _init_H_

#include "discretization.h"
#include "configuration.h"
#include "neurite.h"

extern void setting_rest_pot(TYPE **&potential,discretization **&element);
extern void critical_time_step_class(Hodgkin_Huxley **&HH_elements, int &num_HH_elements, 
				     Cable **&CT_elements, int &num_CT_elements);
extern TYPE critical_time_step_selection(discretization **&element);
extern void flags_init(char *&args,params_t *params);
extern void test_logical_options(params_t *params,input &in);
extern void deleting_discretization(discretization **&element, Hodgkin_Huxley **&HH_elements,
				    int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements);
extern void inputInit(input &in ,  params_t *params,discretization **&element, Hodgkin_Huxley **&HH_elements,
		      int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements);
extern void neuriteInit( params_t *params);
extern void electrical_properties(Hodgkin_Huxley **&HH_elements,int  &num_HH_elements, TYPE &E_MAX);
extern void hhIconds(TYPE *&neurite,Hodgkin_Huxley **&HH_elements,int  &num_HH_elements,TYPE &E_MAX);
extern void cable_init_constants(Cable **&CT_elements, int &num_CT_elements, bool flag);
extern void mechanical_model_initialization(input &in,discretization **&element,TYPE &dist_mes);
extern int measurement_point(discretization **&element, TYPE dist_mes);

#endif
