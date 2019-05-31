//
//
// File author(s):  <Julian Andres Garcia Grajales, Jose Maria Pena and Antoine Jerusalem>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


#ifndef _neurite_H_
#define _neurite_H_

#include "configuration.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

/*! Structure with the flags of the simulation */
typedef struct params {
  bool HH_model;
  bool explicit_mode;
  bool implicit_mode;
  bool openGL;
} params_t;

/*! Structure with the variables likely given in the command line */
typedef struct input {
  TYPE input_current, input_time, init_time, end_input, epsilon, nu, eta_ct, eta_hh, dt_factor, dist_mes, node_size, inter_node, E_MAX;
  int IRE, NRE;
} input;

/*! Structure with the probabilities of the ion channels for the Hodgkin and Huxley model */
typedef struct rate_const {
	TYPE m, h, n;
} rate_const;

#endif
