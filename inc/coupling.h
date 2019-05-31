//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _coupl_H_
#define _coupl_H_
#include "neurite.h"
#include "configuration.h"
#include "discretization.h"

extern void hhRate_Const(rate_const &alpha, rate_const &beta, TYPE rest_pot, TYPE potential,TYPE Left_Shif_Na,TYPE Left_Shif_K);

#endif
