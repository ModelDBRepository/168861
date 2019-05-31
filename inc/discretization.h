//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _discretization_H_
#define _discretization_H_
#include "configuration.h"
#include "math.h"

extern const TYPE pi;
/*! Discretization class
*/
class discretization{

 protected:
  int element; /**< \brief Who are you?. Extra information */
  bool HH; /*!< \brief What kind of element are you? */
  bool branching;/*!< \brief Are you at a branching point? */
  TYPE epsilon;/*!< \brief  Strain for the element */
  TYPE epsilon_surface;/*!< \brief  Surface strain for the element */
  TYPE dX; /*!< \brief Element size*/
  int DL; /*!< \brief Left child. If the element is not a branch, only has right child and DL = 0. !!!IMPORTANT!!!! */
  int DR; /*!< \brief Right child. This element always exists*/
  int mother; /*!< \brief All elements have a mother */
  TYPE diameter; /*!< \brief Neurite diameter */
  TYPE thickness; /*!< \brief Membrane thickness */
  TYPE per_m; /*!< \brief Membrane capacitance (F/m) */
  TYPE ro_m; /*!< \brief Membrane resistivity (ohm.m) */
  TYPE  ro_a;/*!< \brief Intracellular (axial) resistivity (ohm-m) */
  TYPE r_a;/*!< \brief Total axial resistance without dX (dX is implicitly in the discretization) and without myelin */
  TYPE r_m;/*!< \brief Total membrane resistance without dX (dX is implicitly in the discretization) and without myelin */
  TYPE c_m;/*!< \brief Total membrane capacitance without dX (dX is implicitly in the discretization) and without myelin */
  TYPE rest_pot; /*!< \brief Resting potential */
  TYPE W; /*!< \brief W value for each kind of element (Please see final new set of equations) */
  TYPE K; /*!< \brief K value for each kind of element (Please see final new set of equations) */
  TYPE input_current;/*!< \brief  Input current. Only the current is applied to the first element (for the other elements input current =  0) so far. But in the future we can put the input wherever we want*/
  TYPE critical_dT;/*!< \brief Critical time step for the element */
 public:

  // Constructor
  discretization();


  // Setters 
  inline void set_element(int i){element = i;};
  inline void set_epsilon(TYPE eps){
    epsilon = eps;
    epsilon_surface = sqrt(1+epsilon)-1;
  };
  inline void set_DL(int daughter){DL = daughter;};
  inline void set_DR(int daughter){DR = daughter;};
  inline void set_mother(int mom){mother = mom;};
  inline void set_diameter(TYPE temp){diameter = temp;};
  inline void set_thickness(TYPE temp){thickness = temp;};
  /*!\brief Total membrane properties Without dX and without myelin
   */
  inline void set_membrane_axial(){
    r_a =  4*ro_a/(pi*diameter*diameter);
    r_m =  (ro_m*thickness)/(pi*diameter);
    c_m = (per_m*pi*diameter)/(thickness);
  };
  inline void set_branching(bool temp){branching = temp;};
  inline void set_input_current(TYPE temp){input_current = temp;};
  inline void set_dX(TYPE temp){dX = temp;};
  inline void set_critical_dT(TYPE temp){critical_dT = temp;};




  inline int get_element(){return element;};
  inline TYPE get_epsilon(){return epsilon;};
  inline double get_epsilon_surface(){return epsilon_surface;};
  inline int get_DL(){return DL;};
  inline  int get_DR(){return DR;};
  inline int get_mother(){return mother;};
  inline bool get_kind_HH(){return HH;};
  inline TYPE get_diameter(){return diameter;};
  inline TYPE get_thickness(){return thickness;};
  inline TYPE get_r_a(){return r_a;};
  inline TYPE get_r_m(){return r_m;};
  inline TYPE get_c_m(){return c_m;};
  inline bool get_branching(){return branching;};
  inline TYPE get_input_current(){return input_current;};
  inline TYPE get_rest_pot(){return rest_pot;};
  inline TYPE get_dX(){return dX;};// Common
  inline TYPE get_critical_dT(){return critical_dT;};
  /*! \brief The getter for W constant is common for all elements. But the setter is particular of each child, because the values for each kind of element are different
   */
  inline TYPE get_W(){return W;}; // Common
  /*! \brief The getter for K constant is common for all elements. But the setter is particular of each child, because the values for each kind of element are different
   */
  inline TYPE get_K(){return K;};// Common

 
  /*! Destructor
   */
  ~discretization(){};
};


/*! Cable class. This is a daughter of the discretization class. Cable is one type of element possible in Neurite
*/
class Cable : public discretization{

 private:
  TYPE diameter_my;  /*!<\brief  Diameter with myelin NOT USED */
  TYPE thickness_my; /*!<\brief  Myelin thickness*/
  TYPE per_my; /*!<\brief  Electrical properties of each myelin layer. Capacitance*/
  TYPE ro_my; /*!<\brief  Electrical properties of each myelin layer. Resistivity*/
  int layers_my; /*!<\brief  Number of myelin layers*/
  TYPE r_my; /*!<\brief  Electrical properties of each myelin layer without dX. Resistance*/
  TYPE c_my; /*!<\brief  Electrical properties of each myelin layer without dX. Capacitance*/
  
 public:
  /*! Default constructor
   */
  Cable();
  Cable(int ele,int dr, int mom, TYPE deltaX); 

  /*!\brief  Total myelin properties Without dX
   */
  inline void set_myelin(){
    TYPE c_my_T=0;
    TYPE temp=0;
    if(layers_my != 0){
      for(int i=1;i<=layers_my;i++){
	temp = temp + 1/(diameter+thickness*2+2*(i-1)*thickness_my);
      }
      c_my_T = thickness_my*temp/(per_my*pi);
      c_my = 1/c_my_T;
      
      r_my = (ro_my*thickness_my*temp)/pi;
    }else{
      _INFO("Your number of myelin layers is zero. Are you sure about that?");
    }
  };

  /*!\brief Overwrite the myelin properties
   */
  inline void set_final_internodal(){
    if(layers_my != 0){
      r_m = r_m + r_my;
      c_m = (c_m*c_my)/(c_my + c_m);
    }else{
      _INFO("Your number of myelin layers is zero. Are you sure about that?");
    }
  };

  /*! \brief Setting W constant for cable theory elements
   */
  void set_W();
  /*! \brief Setting K constant for cable theory elements
   */
  void set_K();
  /*! Destructor
   */
  ~Cable(){};
};


/*! Hodgkin and Huxley class. This is a daughter of the discretization class.  Hodgkin and Huxley is type of element possible in Neurite
 */
class Hodgkin_Huxley : public discretization{

 private:
  TYPE E_Na;/*!<\brief Sodium reversal potential*/
  TYPE E_Na0;/*!<\brief Sodium reversal potential reference*/
  TYPE E_K0;
  TYPE EMAX;/*!<\brief Parameter for the damage at the Na reversal potential*/
  TYPE E_K;/*!<\brief  Potassium reversal potential*/
  TYPE E_L;/*!<\brief  Leakage reversal potential*/
  TYPE E_L0;/*!<\brief  Leakage reversal potential without strain*/
  TYPE sigma_Na;/*!<\brief  Sodium conductance*/
  TYPE sigma_K; /*!<\brief  Potassium conductance*/
  TYPE n_damage;/*!<\brief  Exponent for the damage in the Na reversal potential*/



  TYPE sigma_L;
  TYPE Left_Shift_Na, Left_Shift_K ;
  TYPE g_L; /*!<\brief  Leakage conductance without dX*/
  TYPE g_Na; /*!<\brief Sodium conductance without dX nor probabilities. Constant*/
  TYPE g_K; /*!<\brief Potassium conductance without dX nor probabilities. Constant*/
  TYPE m;/*!<\brief Sodium probability*/
  TYPE h;/*!<\brief Sodium probability*/
  TYPE n;/*!<\brief Potassium probability*/ 
  TYPE G_Na;/*!<\brief Total Sodium conductance. It depends on the time and potential but without dX, that is implicitly in the discretization scheme*/
  TYPE G_K;/*!<\brief Total Potassium conductance. It depends on the time and potential but without dX, that is implicitly in the discretization scheme*/
  TYPE G_L;/*!<\brief Total leakage conductance (constant) without dX, that is implicitly in the discretization scheme*/
  
 public:
  /*! Default constructor
   */
  Hodgkin_Huxley();
  Hodgkin_Huxley(int ele,int dr, int mom, TYPE deltaX);

  //Setters
  inline void set_E_L(TYPE temp){
    E_L = temp;
    E_L0 = temp;
  };
  /*! \brief Setting the conductances without dX nor probabilities. Constant values.
   */
  inline void set_conductances(){
    //g_L = 1/r_m;
    g_L = sigma_L*pi*diameter/thickness;
    g_Na = sigma_Na*pi*diameter/thickness;
    g_K = sigma_K*pi*diameter/thickness;
  };
  /*! \brief Damaging the reversal potentials
   */
  inline void modifying_reversal_potentials(TYPE temp){
    EMAX = temp;
    if(epsilon_surface<=temp){
      E_Na = E_Na0-pow((epsilon_surface/temp),n_damage)*(E_Na0);
      E_K = E_K0-pow((epsilon_surface/temp),n_damage)*(E_K0);
    }else{
      _INFO("WARNING!!! CAREFUL, your EMAX (%f) is smaller than your surface strain (%f)",temp,epsilon_surface);
      // This way is because it is impossible to have an epsilon surface bigger than the EMAX. At the limit, the potentials are equal to the resting potential of the membrane.
      E_Na = 0;
      E_K = 0;
    }
  };
  /*! \brief This function calculates the final left-shifts
   */
  inline void calculating_Left_Shift(){
    if(epsilon_surface<=EMAX){
      Left_Shift_Na = pow(epsilon_surface/EMAX,n_damage)*(E_Na0);
      Left_Shift_K = pow(epsilon_surface/EMAX,n_damage)*(E_K0);
    }else{
      Left_Shift_Na = E_Na0;
      Left_Shift_K = E_K0;
    }    
  };
  inline void set_g_L(TYPE temp){g_L = temp;};
  inline void set_g_Na(TYPE temp){g_Na = temp;};
  inline void set_g_K(TYPE temp){g_K = temp;};
  inline void set_m(TYPE temp){m = temp;};
  inline void set_h(TYPE temp){h = temp;};
  inline void set_n(TYPE temp){n = temp;};
  inline void set_G_Na(TYPE temp){G_Na = temp;};
  inline void set_G_K(TYPE temp){G_K = temp;};
  inline void set_G_L(TYPE temp){G_L = temp;};

  // Getters
  inline TYPE get_m(){return m;};
  inline TYPE get_h(){return h;};
  inline TYPE get_n(){return n;};
  inline TYPE get_E_Na(){return E_Na;};
  inline TYPE get_E_Na0(){return E_Na0;};
  inline TYPE get_E_K0(){return E_K0;};
  inline TYPE get_EMAX(){return EMAX;};
  inline TYPE get_n_damage(){return n_damage;};
  inline TYPE get_E_K(){return E_K;};
  inline TYPE get_E_L(){return E_L;};
  inline TYPE get_g_L(){return g_L;};
  inline TYPE get_g_Na(){return g_Na;};
  inline TYPE get_g_K(){return g_K;};
  inline TYPE get_G_Na(){return G_Na;};
  inline TYPE get_G_K(){return G_K;};
  inline TYPE get_G_L(){return G_L;};
  inline TYPE get_Left_Shift_Na(){return Left_Shift_Na;};
  inline TYPE get_Left_Shift_K(){return Left_Shift_K;};

  /*! \brief Setting W constant for Hodgkin and Huxley elements
   */
  void set_W();
  /*! \brief Setting K constant for Hodgkin and Huxley elements
   */
  void set_K();
 /*! Destructor
   */
  ~Hodgkin_Huxley(){  };

};

#endif

