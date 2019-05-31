//
//
// File author(s):  <Julian Andres Garcia Grajales and Gabriel Rucabado>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _SOLVER
#define _SOLVER

class Solver {
 public:
  //!\brief This function updates the electrical variables based on the previous time step information
  //  @param input_active Flag used to indicate if it has to use  the input current to update the values
  virtual void update (bool input_active) = 0;
  //!\brief Calculate one iteration of the explicit finite difference discretization solver
  virtual void calculate() = 0;
   //! \brief Copy the current potential to the pointer of the parameter
  //  @param potential Pointer on which the potential will be copied 
  virtual void snapshot(TYPE *potential) = 0;
  //! \brief Free all memory used by the solver
  virtual ~Solver(){};
};

#endif
