//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file init.cpp
  \brief Several functions to initialize all variables
*/

#include "init.h"
#include "coupling.h"
#include "space.h"
#include "mechanical_model.h"
#include "configuration.h"


#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

/* global variables */
extern int numStepsX, open_GL_steps;
extern unsigned numStepsT;
extern TYPE xspan, tspan, dist_mes;
extern TYPE *neuriteE, **neuriteR, ***neurite, *open_GL_graphE,**open_GL_graphR,***open_GL_graph;
extern TYPE **potential, *potentialR;
extern TYPE ***act_trans, **props, *ini;
extern int tot_santa,alloc_santa, num_props;
extern int **n_micro, *n_microE;
extern int *terminals, num_terminals;
extern char *simulation;

static void print_electrical_discretization(discretization **&element);

/*!  This function defines the type of simulation that you are running
  @param args Arguments from the command line
  @param params Structure with the flags of the simulation
*/
extern void flags_init(char *&args,params_t *params){

  if (memchr(args, 'e', strlen(args)) != NULL)
    {
      //EXPL_MODE
      params->explicit_mode = true;
    }
  if (memchr(args, 'i', strlen(args)) != NULL)
    {
      //IMPL_MODE
      params->implicit_mode = true;
    }
  if (memchr(args, 'g', strlen(args)) != NULL)
    {
      //openGL
      params->openGL = true;
    }

}

/*!  This function checks that the input flags make sense
  @param params Structure with the flags of the simulation
*/
extern void test_logical_options(params_t *params, input &in)
{
  if(params->implicit_mode && params->explicit_mode){
    _ERROR("ERROR. You are using both solvers at the same time\n");
    exit(0);
  }
  if(!params->implicit_mode && !params->explicit_mode){
    _ERROR("ERROR. You did not take any solver\n");
    exit(0);
  }
}

/*!  This function defines the spatial discretization. The user must put the corresponding function that defines his spatial discretization following the Neurite rules. You must put the function at the corresponding place inside this function.
  @param in Structure with parameters defined in the command line
  @param params Structure with the flags of the simulation
  @param element Array of discretization classes
  @param HH_elements This is an array of classes that contains the Hodgkin and Huxley elements in the discretization
  @param num_HH_elements Size of the HH_elements array
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array
*/
extern void inputInit(input &in, params_t *params, discretization **&element, Hodgkin_Huxley **&HH_elements, int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements)
{

  if(strcmp(simulation,"Only_HH_Model") == 0){
    _INFO("WARNING. You are using only the Hodgkin Huxley option. This option is to simulate the behavior of ONLY the HH model.");
    params->HH_model = true;
    if(params->implicit_mode){
      _ERROR("The implicit solver is not available for the isolated ion channel\n");
      exit(0);
    }
    // Special function to run only the HH model surrounded by two fictitious elements
    Only_HH_Model_class(in, element, HH_elements, num_HH_elements,CT_elements, num_CT_elements);
    params->HH_model = true;
    print_electrical_discretization(element);
  }else{ // Electrical model

      if (strcmp(simulation,"HH_Y_branching") == 0){
	// CT+HH with branching
	branching_Y_CT_HH(in, element, HH_elements, num_HH_elements,CT_elements, num_CT_elements); 
	params->HH_model = true;
      }else if(strcmp(simulation,"HH_axon") == 0){
	// CT+HH without branching
	without_branching_CT_HH(in,element, HH_elements, num_HH_elements, CT_elements, num_CT_elements);
	params->HH_model = true;
      }else if (strcmp(simulation,"CT_Y_branching") == 0){
	branching_Y_CT(in, element,CT_elements, num_CT_elements);  // CT with branching
      	HH_elements = NULL;
	num_HH_elements = 0;
      }else if(strcmp(simulation,"CT_cable") == 0){
	without_branching_CT(element, in,CT_elements, num_CT_elements); // CT without branching
      	HH_elements = NULL;
	num_HH_elements = 0;
      }else if(strcmp(simulation,"CT_symmetric_tree") == 0){
	general_symmetric_tree_CT(in, element, CT_elements, num_CT_elements);
	HH_elements = NULL;
	num_HH_elements = 0;
      }else if(strcmp(simulation,"CT_segmented_neuron") == 0){
	loading_neuron_passive_dendritic_tree(in, element, CT_elements, num_CT_elements);
	HH_elements = NULL;
	num_HH_elements = 0;
      }else{
	_ERROR("ERROR. You use a passive cable but did not define the specific case to simulate: branching Y, cable, symmetric tree or real segmented neuron\n");
	exit(0);
      }
  }
  print_electrical_discretization(element);
}

/*! This function calculates the critical time step for each element according to what kind of element you are
  @param HH_elements This is an array of classes that contains the Hodgkin and Huxley elements in the discretization
  @param num_HH_elements Size of the HH_elements array
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array
*/
extern void critical_time_step_class(Hodgkin_Huxley **&HH_elements, int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements){

  for(int i=0;i<num_HH_elements;i++){
      HH_elements[i]->set_critical_dT((2*HH_elements[i]->get_r_a()*HH_elements[i]->get_c_m()*HH_elements[i]->get_dX()*HH_elements[i]->get_dX())/4);
  }
  for(int i=0;i<num_CT_elements;i++){
      CT_elements[i]->set_critical_dT((2*CT_elements[i]->get_r_a()*CT_elements[i]->get_c_m()*CT_elements[i]->get_dX()*CT_elements[i]->get_dX())/(4 + (CT_elements[i]->get_r_a()*CT_elements[i]->get_dX()*CT_elements[i]->get_dX())/CT_elements[i]->get_r_m()));
  }
}

/*!  Deleting the whole discretization
  @param element Array of discretization classes
  @param HH_elements This is an array of classes that contains the Hodgkin and Huxley elements in the discretization
  @param num_HH_elements Size of the HH_elements array
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array
*/
extern void deleting_discretization(discretization **&element, Hodgkin_Huxley **&HH_elements, int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements){
    for(int i=0;i<numStepsX;i++){
      delete element[i];
    }
}

/*!  This function Chooses the most critical time step between the elements
  @param element Array of discretization classes
  @return temp The critical time step calculated
*/
TYPE critical_time_step_selection(discretization **&element){
  TYPE temp = 0;
  temp = element[1]->get_critical_dT();
  for(int i=2;i<numStepsX-1;i++){
    if(temp > element[i]->get_critical_dT()){
      temp = element[i]->get_critical_dT();
    }
  }
  return temp;
}

/*!  This function Initializes the cable theory electrical properties
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array
*/
extern void cable_init_constants(Cable **&CT_elements, int &num_CT_elements, bool flag){

  for(int i=0;i<num_CT_elements;i++){
      // Is there strain?
      CT_elements[i]->set_diameter(CT_elements[i]->get_diameter()/(sqrt(1+CT_elements[i]->get_epsilon())));
      // Setting the electrical properties for the membrane and cytoplasm   
      CT_elements[i]->set_membrane_axial();
      // Setting the electrical properties for the myelin
      // Is there myelin?
      if(flag == true){
	// Total properties for the cable theory part-> Membrane + myelin-> Internodal part
	// Within internodal distance-> R_axon=R_mem+d->layers*R_myelin,  1/cc=1/c_m+layers/c_my
	CT_elements[i]->set_myelin();
	CT_elements[i]->set_final_internodal(); 
      }
      CT_elements[i]->set_W();
      CT_elements[i]->set_K();
  }
}

/*!  This function allocates the memory needed for the simulation
  @param params Structure with the flags of the simulation
*/
extern void neuriteInit(params_t *params)
{
  int h, i;
  int old_new=2;
  int numCurves=1;

      potentialR = new TYPE[old_new*numStepsX]();
      potential = new TYPE*[old_new]();

     for(i = 0; i < old_new; i++){
       potential[i] = potentialR + (i*numStepsX);
     }

     if(numStepsX > 300000){
       open_GL_steps =2;
     }
     open_GL_graphE = new TYPE[numCurves*open_GL_steps*numStepsX]();
     open_GL_graphR = new TYPE*[numCurves*open_GL_steps]();  
     open_GL_graph = new TYPE**[numCurves]();  
     for(h = 0; h < numCurves; h++){
	  for(i = 0; i < open_GL_steps; i++){
	    open_GL_graphR[i+h*open_GL_steps] =   open_GL_graphE + (i*numStepsX) + (h*2*numStepsX);
	  }
	  open_GL_graph[h] = open_GL_graphR + (h*open_GL_steps);
     }  
}


/*!  This function defines the electrical properties for Hodgkin and Huxley elements
  @param HH_elements This is an array of classes that contains the Hodgkin and Huxley elements in the discretization
  @param num_HH_elements Size of the HH_elements array
*/
extern void electrical_properties(Hodgkin_Huxley **&HH_elements,int  &num_HH_elements, TYPE &E_MAX){
  for(int i=0; i<num_HH_elements; i++){

      // Is there strain?
      HH_elements[i]->set_diameter(HH_elements[i]->get_diameter()/(sqrt(1+HH_elements[i]->get_epsilon())));
      // Setting the membrane properties and cytoplasm with out myelin
      HH_elements[i]->set_membrane_axial();
      //modifying the reversal potentials
      HH_elements[i]->modifying_reversal_potentials(E_MAX);
      // Conductances without dX
      HH_elements[i]->set_conductances();
     // Calculating the left shift. EMAX is already inside
      HH_elements[i]->calculating_Left_Shift();
      
      // Is there strain? In order to normalize the conductances. Independent on the strain
      HH_elements[i]->set_g_Na(HH_elements[i]->get_g_Na()/(sqrt(1+HH_elements[i]->get_epsilon())));
      HH_elements[i]->set_g_K(HH_elements[i]->get_g_K()/(sqrt(1+HH_elements[i]->get_epsilon())));

  }
}

/*!  This function sets the resting potential for everybody
   @param potential Array with the current/old potential for all elements
   @param element Array of discretization classes
*/
extern void setting_rest_pot(TYPE **&potential,discretization **&element)
{
  int  j=0;
  for(j = 0; j < numStepsX; j++){
    potential[0][j] = element[j]->get_rest_pot();
    potential[1][j] = element[j]->get_rest_pot();
  }
}

/*!  This function initializes the probabilities for Hodgkin and Huxley elements
  @param HH_elements This is an array of classes that contains the Hodgkin and Huxley elements in the discretization
  @param num_HH_elements Size of the HH_elements array
*/
extern void hhIconds(TYPE *&neurite,Hodgkin_Huxley **&HH_elements,int  &num_HH_elements, TYPE &E_MAX)
{
  int i;
  TYPE m, h, n;
  rate_const alpha, beta;
  TYPE V_lost;

  for(i = 0; i < num_HH_elements; i++){
 
    hhRate_Const(alpha, beta, HH_elements[i]->get_rest_pot(),HH_elements[i]->get_rest_pot(),HH_elements[i]->get_Left_Shift_Na(),HH_elements[i]->get_Left_Shift_K());

    // The first time the rate is equal to 0
    m = alpha.m/(alpha.m+beta.m);
    h = alpha.h/(alpha.h+beta.h);
    n = alpha.n/(alpha.n+beta.n);

    HH_elements[i]->set_m(m);
    HH_elements[i]->set_h(h);
    HH_elements[i]->set_n(n);

    V_lost = (HH_elements[i]->get_g_Na()*m*m*m*h * (HH_elements[i]->get_rest_pot()-HH_elements[i]->get_E_Na()) + HH_elements[i]->get_g_K() *n*n*n*n* (HH_elements[i]->get_rest_pot()-HH_elements[i]->get_E_K())+ HH_elements[i]->get_g_L()*HH_elements[i]->get_rest_pot())/HH_elements[i]->get_g_L();

    HH_elements[i]->set_E_L(V_lost);
    HH_elements[i]->set_G_L(HH_elements[i]->get_g_L()); // This does not change during the time
  }
}

/*!  Printing the discretization of the electrical simulation
  @param element Array of discretization classes
*/
static void print_electrical_discretization(discretization **&element){
    for(int i=0; i<numStepsX;i++){
      _DEBUG ("I am %i. My mother is %i and my daughters DR %i DL %i HH=%i",
	      element[i]->get_element(), element[i]->get_mother(),
	      element[i]->get_DR(), element[i]->get_DL(),element[i]->get_kind_HH());
    }
    _DEBUG ("The daughter right of each terminal element i must be i (himself). MANDATORY");
    _DEBUG ("Initial fictitious element %i DR %i", terminals[0],element[terminals[0]]->get_DR());
    for(int i=1; i<num_terminals;i++){
      _DEBUG ("Terminal %i DR %i", terminals[i],element[terminals[i]]->get_DR());
    }
    _INFO ("Total steps in space=%i", numStepsX);
}

/*!  Modify the discretization accordingly to the microscopic strain
  @param element Array of discretization classes
*/
extern void mechanical_model_initialization(input &in,discretization **&element, TYPE &dist_mes)
{
  TYPE epsilon_macro=0;
  if(in.epsilon > 0){
    calculate_epsilon_macro(in.epsilon, epsilon_macro);
    dist_mes= dist_mes*(1+epsilon_macro);
    _INFO("Epsilon_axial_membrane %g epsilon macroscopic %g and the measurement distance is %g",in.epsilon, epsilon_macro,dist_mes);
  }
  // strain for the model. Each element will have its own strain
  for(int k=0;k<numStepsX;k++){
    element[k]->set_dX(element[k]->get_dX()*(1+element[k]->get_epsilon()));
  }
  xspan = xspan*(1+in.epsilon); // Epsilon membrane axial, because all calculations were made in function of that numb
}


/*!  Define the measurement point
  @param element Array of discretization classes
*/
extern int measurement_point(discretization **&element, TYPE dist_mes)
{
  int dist_mes_points=0;
  TYPE  temp=0;
  for(int k=0;k<numStepsX;k++){
    temp = temp + element[k]->get_dX();
    if(temp >= dist_mes){
      dist_mes_points = k;
      break;
      temp =-1; // Only write one time
    }
  }
  if(dist_mes_points==0){
    dist_mes_points = int(numStepsX/2);
  }
  _INFO("Total number of elements %i",numStepsX);
  _INFO("Measurement_point %i at distance %g",dist_mes_points,dist_mes);

  return dist_mes_points;
}
