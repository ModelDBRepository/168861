//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


//////////////// Space discretization //////////////////////////////////
/*
 The user needs to define the discretization. Here there are some examples for models without branching (CT and CT+HH),
 branching model with only Y shape (CT and CT-HH), a symmetric tree with n terminals and with cable theory elements, or an isolated HH element.
*/
//////////////// Space discretization //////////////////////////////////
#ifndef _space_H_
#define _space_H_
#include "neurite.h"
#include "discretization.h"

/*! Structure to load a segmented neuron*/
typedef struct  neuron_input{
  int me,mother,dR, dL,branching;
  double diameter, length;    
} neuron_input;


extern void loading_neuron_passive_dendritic_tree(input &in, discretization **&element, Cable **&CT_elements, int &num_CT_elements);

// This function gives us the discretization for the isolated HH element
extern void Only_HH_Model_class(input &in, discretization **&element, Hodgkin_Huxley **&HH_elements,
				  int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements);

// This function gives us the discretization for the model without branching (CT)
extern void without_branching_CT(discretization **&element,input &in, Cable **&CT_elements, int &num_CT_elements);

// This function gives us the discretization for the model without branching (CT+HH)
extern void without_branching_CT_HH(input &in, discretization **&element, Hodgkin_Huxley **&HH_elements,
				    int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements);

// This function gives us the discretization for the Y shape with CT and HH elements and its corresponding internodal distances and node size
extern void branching_Y_CT_HH(input &in, discretization **&element, Hodgkin_Huxley **&HH_elements,
			      int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements);

// This function gives us the discretization for the Y shape with only cable theory elements
extern void branching_Y_CT(input &in, discretization **&element, Cable **&CT_elements, int &num_CT_elements);

// This function gives us the discretization of a general symmetric passive tree
extern void general_symmetric_tree_CT(input &in, discretization **&element, Cable **&CT_elements, int &num_CT_elements);


#endif
