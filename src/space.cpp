//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


/*! 
 * \file space.cpp
 * \brief Space discretization.
 *\details This file must be edited by the user. Users give to Neurite the corresponding discretization and it is responsibility of the user to give a correct discretization filling the corresponding arrays properly. The functions shown in this file are only examples with several particular conditions.
*/
#include "space.h"
#include "configuration.h"

extern int numStepsX;
extern unsigned numStepsT;
extern double xspan, tspan, dist_mes;
extern int *terminals, num_terminals;
extern FILE *standard, *study;

int compareints2 (const void * x, const void * y)
{
  return ( *(int*)x - *(int*)y );
}

/*! Loading a neuron from the txt file. In this case, it is a passive dendritic tree.
  @param in Structure with parameters defined in the command line
  @param element Array of discretization classes
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array

  * You have to define the neuron_for_neurite.txt file in the input folder.
 */
extern void loading_neuron_passive_dendritic_tree(input &in, discretization **&element, Cable **&CT_elements, int &num_CT_elements){

  in.input_time = 0.003; // To apply the external intensity
  in.input_current = 0.8E-9; // Amperes. External Intensity
  in.end_input = in.init_time + in.input_time;
  num_terminals = 0;
  xspan = 0; // In this case this value does not make sense
  std::vector<neuron_input> neuron;  
  // Loading the file
  int number_raws=0,i=0;
  char fname[250];
  char *execution_folder = getenv("working_directory");
  if (execution_folder != NULL){
    strcpy (fname, execution_folder);
    strcat (fname, "/neuron_for_neurite.txt");
  }else{
    fprintf(stdout,"ERROR. You did not define the working directory variable\n");
    exit(0);
  }
  //Open file
  std::ifstream inFile (fname);
 
 // Refining each element
  int refining=1; // 1-> No refinement
  // Cont terminals
  int cont_terminals=1; // This is the initial 0 element
  double micrometer=1E-6;
  
  
  if (inFile.is_open()) {
    while(!inFile.eof()){
      neuron.push_back(neuron_input());
      inFile >> neuron[i].me >> neuron[i].mother >> neuron[i].dR >> neuron[i].branching >> neuron[i].dL >> neuron[i].diameter >> neuron[i].length;
      xspan+=neuron[i].length;
      if(neuron[i].me == neuron[i].dR){
	cont_terminals++;
      }
      i++;
    } 
    number_raws = i-1;
  }else{
    // show message:
    _ERROR("Error opening file with the corresponding segmented neuron. See the space.cpp file in the loading_neuron function");
    exit(0);
  }
  num_terminals =  cont_terminals-1; 

  // Loading the neuron structure from the swc format. Pre edited with matlab to take only one branch with several branches.

 if(tspan < 1E-10){
    // Time discretization
    tspan = 0.025; // Total time
  } 

 numStepsX=number_raws*refining+1; // Total number of elements. // +1 is for the initial element
 int k=0;

  element =  new discretization*[numStepsX];
  terminals = new int[num_terminals](); // () initialize to zero
  num_CT_elements = numStepsX;
  // In this case everyone is cable theory
  num_CT_elements = numStepsX;
  CT_elements = new Cable*[num_CT_elements];
  int pos=0,cont_terminal=1;

  // Creating the fictitious 0 element
  terminals[0] = 0; // First terminal is always the fictitious element.
  pos=0;
  k++;
  CT_elements[pos]=new Cable(pos, pos + 1, -1, neuron[0].length/refining*micrometer);
  element[pos] = CT_elements[pos]; 
  element[pos]->set_epsilon(in.epsilon);
  element[pos]->set_DL(0); // No branching
  element[pos]->set_branching(false);
  element[pos]->set_diameter(neuron[0].diameter*micrometer);

 for(i=0;i<number_raws;i++){
   // If you are branching
   if(neuron[i].branching){
     for(int j=1;j<=refining;j++){
       pos = (neuron[i].me-1)*refining + j;
       CT_elements[pos]=new Cable(pos, pos + 1, neuron[i].mother*refining + j - 1, neuron[i].length/refining*micrometer);
       element[pos] = CT_elements[pos]; 
       element[pos]->set_epsilon(in.epsilon);
       if(j==refining){
	 element[pos]->set_DL(neuron[i].dL*refining - refining + 1); // branching
	 element[pos]->set_branching(true);
       }else{
	 element[pos]->set_DL(0); // branching
	 element[pos]->set_branching(false);
       }
       element[pos]->set_diameter(neuron[i].diameter*micrometer);
       k++;
     }   
   }else if(neuron[i].me == neuron[i].dR){   // If you are terminal
     for(int j=1;j<=refining;j++){
       pos = (neuron[i].me-1)*refining + j;
       if(j==refining){
       CT_elements[pos]=new Cable(pos, pos, neuron[i].mother*refining + j - 1, neuron[i].length/refining*micrometer);// The last terminal is yourself
       }else{
       CT_elements[pos]=new Cable(pos, pos+1, neuron[i].mother*refining + j - 1, neuron[i].length/refining*micrometer);
       }
       element[pos] = CT_elements[pos]; 
       element[pos]->set_epsilon(in.epsilon);
       element[pos]->set_DL(0); // No branching
       element[pos]->set_branching(false);
       element[pos]->set_diameter(neuron[i].diameter*micrometer);

       k++;
     }
     terminals[cont_terminal] = pos;
     cont_terminal++;
   }else{
     if(i>0 && neuron[i-1].me == neuron[i-1].dR){ // But the previous one was terminal
       for(int j=1;j<=refining;j++){
	 pos = (neuron[i].me-1)*refining + j;
	 if(j==1){
	   CT_elements[pos]=new Cable(pos, pos + 1, neuron[i].mother*refining + j - 1, neuron[i].length/refining*micrometer);
	 }else{
	   CT_elements[pos]=new Cable(pos, pos + 1, pos-1, neuron[i].length/refining*micrometer);
	 }
	 element[pos] = CT_elements[pos]; 
	 element[pos]->set_epsilon(in.epsilon);
	 element[pos]->set_DL(0); // No branching
	 element[pos]->set_branching(false);
	 element[pos]->set_diameter(neuron[i].diameter*micrometer);

	 k++;
       }
     }else{ // Normal normal
       for(int j=1;j<=refining;j++){
	 pos = (neuron[i].me-1)*refining + j;
	 CT_elements[pos]=new Cable(pos, pos + 1, neuron[i].mother*refining + j - 1, neuron[i].length/refining*micrometer);
	 element[pos] = CT_elements[pos]; 
	 element[pos]->set_epsilon(in.epsilon);
	 element[pos]->set_DL(0); // No branching
	 element[pos]->set_branching(false);
	 element[pos]->set_diameter(neuron[i].diameter*micrometer);
	 k++;
       }
     }
   }
 }
 dist_mes=100E-6;
 fprintf(standard,"Total steps = %i ", numStepsX);
 element[1]->set_input_current(in.input_current);
}


/*! General tree with n jumps of cable theory elements. Passive cable.
  @param in Structure with parameters defined in the command line
  @param element Array of discretization classes
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array

  * Please, keep in your mind that you must define a new function following the way of the discretization. This function is only an example.
 */
extern void general_symmetric_tree_CT(input &in, discretization **&element, Cable **&CT_elements, int &num_CT_elements){

  TYPE dXtemp=0;
  in.input_time = 0.003; // To apply the external intensity
  in.input_current = 0.8E-9; // Nano Amperes. External Intensity
  in.init_time=0;
  in.end_input = in.init_time + in.input_time;   // when you finish to apply the current
  num_terminals = 0;
  //int min_size_branch = 70, max_size_branch = 70;// number of elements in each branch
  int *size_branch; // A pointer to know the size of the tree
  int **flag_jump, *flag_jumpE;
  int k=0;
  xspan=0;

  if(tspan < 1E-10){
    // Time discretization
    tspan = 0.025; // Total time
  } // if not is because you put as an argument.

    // All elements are CT elements. We define the element size and begin the tree in function of the elements at each branch
  dXtemp = 10E-6; //
  numStepsX=0;

  int trams=0, cont_trams=0, initial_sections = 4; // This number defines the number of jumps from the beginning.

  num_terminals = 1*pow(2,initial_sections-1) + 1;  // The last 1 is for the first fictitious element
  trams = (num_terminals-1)*2-1; // Experience

  int temp_numStepsX=0;
  size_branch = new int[trams]();// () initialize to zero
  flag_jump = new int*[trams](); // At most This number must be adjusted. THIS ALGORITHM IS TEMPORAL. flag_jump[x][y]. x is the position in the vector and y=1 if already has branch, and y=0 vice versa.
  flag_jumpE = new int[trams*2]();
   for(int i = 0; i < trams; i++){
	flag_jump[i] = flag_jumpE + (i*2);
      }

   srand(time(NULL)); // To generate random numbers. In this case, the number of elements of each branch
  // for(k=0;k<trams;k++){
  //   size_branch[k] = min_size_branch + rand()%(max_size_branch + 1 - min_size_branch);
  //   temp_numStepsX = temp_numStepsX +  size_branch[k];
  // }

   // For a fixed tree // Initial sections=4;
   size_branch[0]=44;
   size_branch[1]=62;
   size_branch[2]=48;
   size_branch[3]=50;
   size_branch[4]=28;
   size_branch[5]=50;
   size_branch[6]=26;
   size_branch[7]=32;
   size_branch[8]=48;
   size_branch[9]=50;
   size_branch[10]=30;
   size_branch[11]=32;
   size_branch[12]=40;
   size_branch[13]=55;
   size_branch[14]=38;
   temp_numStepsX = 633;

  element =  new discretization*[temp_numStepsX];
  terminals = new int[num_terminals](); // () initialize to zero

  // In this case everyone is cable theory
  num_CT_elements = temp_numStepsX ;
  CT_elements = new Cable*[num_CT_elements];

  terminals[0] = 0; // First terminal is always the fictitious element.
  k=0;
  int j=0;
  int jumps_branch = initial_sections;
  int cont_elem_tram=0, cont_jumps_tram=0, cont_back=0, cont_inside=0, cont_jumps=0;
  numStepsX=0;
  j=1;
  int mom=0;
  while(cont_trams < trams){
    cont_elem_tram = 0;
    while(cont_elem_tram < size_branch[k]){
      if(cont_elem_tram == size_branch[k]-1){
	if(cont_jumps_tram == jumps_branch-1){
	  terminals[j]=numStepsX;
	  j++;
	  CT_elements[numStepsX]=new Cable(numStepsX, numStepsX,numStepsX-1, dXtemp); //
	  element[numStepsX] = CT_elements[numStepsX]; 
	  element[numStepsX]->set_epsilon(in.epsilon);
	  cont_back = 0;
	  cont_inside=1;
	  while(cont_inside<cont_jumps){
	    cont_back++;
	    if(flag_jump[cont_jumps-1-cont_inside][1] == 0){
	      flag_jump[cont_jumps-1-cont_inside][1]=1; 
	      mom = flag_jump[cont_jumps-1-cont_inside][0];
	      element[flag_jump[cont_jumps-1-cont_inside][0]]->set_DL(numStepsX+1);
	      cont_jumps=cont_jumps-cont_inside;
	      cont_inside=10000;
	    }
	    cont_inside++;
	  }
	  cont_jumps_tram=0;
	  jumps_branch = cont_back;
	}else{
	  CT_elements[numStepsX]=new Cable(numStepsX, numStepsX+1,numStepsX-1, dXtemp);
	  element[numStepsX] = CT_elements[numStepsX]; 
	  element[numStepsX]->set_epsilon(in.epsilon);
	  element[numStepsX]->set_branching(true);
	  cont_jumps_tram++;
	}
      }else{
	if(mom == 0){
	  CT_elements[numStepsX]=new Cable(numStepsX, numStepsX+1, numStepsX-1, dXtemp);
	  element[numStepsX] = CT_elements[numStepsX]; 
	  element[numStepsX]->set_epsilon(in.epsilon);
	}else{
	  CT_elements[numStepsX]=new Cable(numStepsX, numStepsX+1, mom, dXtemp);
	  mom=0;
	  element[numStepsX] = CT_elements[numStepsX]; 
	  element[numStepsX]->set_epsilon(in.epsilon);
	}
      }
      numStepsX++;
      cont_elem_tram++;
      if(cont_elem_tram == size_branch[k]-1){
	flag_jump[cont_jumps][0]=numStepsX;
	flag_jump[cont_jumps][1]=0;
	cont_jumps++;
      }
    }
    cont_trams++;
    k++;
  }
  // Measurement distance
  dist_mes = dXtemp*numStepsX/10;

  fprintf(standard,"total elements steps=%i === %i\n", numStepsX, temp_numStepsX);
  element[1]->set_input_current(in.input_current); // Only at the initial element

}


/*! General tree with 2 jumps of cable theory elements. Passive cable. The tree is Y. This function does exactly the same than general_tree_n_terminal_CT with n=2 jumps
  @param in Structure with parameters defined in the command line
  @param element Array of discretization classes
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array

  * Please, keep in your mind that you must define a new function following the way of the discretization. This function is only an example.
 */
extern void branching_Y_CT(input &in, discretization **&element, Cable **&CT_elements, int &num_CT_elements){

  TYPE dXtemp=0;
  in.input_time = 0.003; // To apply the external intensity
  in.input_current = 0.8E-9; // Nano Amperes. External Intensity
  in.end_input = in.init_time + in.input_time;   // when you finish to apply the current
  num_terminals =3; 
  int size_branch = 50;// number of elements for each branch.
  int temp_cont=0;
  int temp_terminals=0;

  int flag_jump=0;
  int temp_path=0; 
  xspan=0;

  if(tspan < 1E-10){
    // Time discretization
    tspan = 0.025; // Total time
  } // if not is because you put as an argument.

  
  dXtemp = 50E-6; // 50 elements for each branch with this size
  numStepsX=0;

  element =  new discretization*[size_branch*num_terminals];
  terminals = new int[num_terminals](); // () initialize to zero

  num_CT_elements = size_branch*num_terminals ;
  CT_elements = new Cable*[num_CT_elements];

  terminals[0] = 0;
  while(temp_terminals < num_terminals){
    temp_cont = 0;
    while(temp_cont<size_branch){
      temp_cont++;
      if(temp_terminals==2){
	terminals[temp_terminals-1]=numStepsX-1;

	CT_elements[numStepsX]=new Cable(numStepsX, numStepsX+1,flag_jump, dXtemp); //my mother is flag_jump
	element[numStepsX] = CT_elements[numStepsX]; 
	element[numStepsX]->set_epsilon(in.epsilon);
	element[flag_jump]->set_DL(numStepsX);
	element[flag_jump]->set_branching(true);
	element[numStepsX-1]->set_DR(numStepsX-1);
	    
	temp_terminals++;

      }else{
	CT_elements[numStepsX]=new Cable(numStepsX, numStepsX+1,numStepsX-1, dXtemp);
	element[numStepsX] = CT_elements[numStepsX]; 
	 element[numStepsX]->set_epsilon(in.epsilon);
      }
      numStepsX++;
    }
    if(temp_path==0){
      flag_jump=temp_cont-1;
    }
    temp_path++;
    temp_terminals++;
  }
  terminals[num_terminals-1] = numStepsX-1;
  element[numStepsX-1]->set_DR(numStepsX-1);
  // Measurement distance
  dist_mes = numStepsX*dXtemp/10;
  element[1]->set_input_current(in.input_current); // Only at the initial element
}

/*! General tree with 2 jumps of CT+HH elements. Active cable. The tree is Y
  @param in Structure with parameters defined in the command line
  @param element Array of discretization classes
  @param HH_elements This is an array of classes that contains the Hodgkin and Huxley elements in the discretization
  @param num_HH_elements Size of the HH_elements array
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array

  * Please, keep in your mind that you must define a new function following the way of the discretization. This function is only an example.
 */
extern void branching_Y_CT_HH(input &in, discretization **&element, Hodgkin_Huxley **&HH_elements, int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements){

  int hh_steps_temp=0, ct_steps_temp=0;
  TYPE ct_dXtemp=0, hh_dXtemp=0;
  TYPE node_size = 2.1E-6; // Node of Ranvier. From Shi papers
  TYPE inter_node = 0.0008; // internodal length. General distance
  num_terminals=3; // Initial + two final elements. Y.

  in.input_time = 0.003; // To apply the external intensity
  in.input_current = 0.8E-9; // Nano Amperes. External Intensity
  in.end_input = in.init_time + in.input_time;   // when you finish to apply the current

  if(in.IRE< 1E-10){
    // Default values for the discretization
    ct_steps_temp = 10; // internodal part
  }else{
    ct_steps_temp = in.IRE; // internodal part
  }
  if(in.NRE< 1E-10){
    // Default values for the discretization
    hh_steps_temp = 1; // internodal part
  }else{
    hh_steps_temp = in.NRE; // internodal part
  }

  if(tspan < 1E-10){
    // Time discretization
    tspan = 0.025; // Total time
  } // if not is because you put as an argument.

  //Calculate dXi, dXn and dT
  ct_dXtemp = (TYPE)inter_node/ct_steps_temp;
  hh_dXtemp = (TYPE)node_size/hh_steps_temp;

  fprintf(standard,"for discretization and simulation cabledx =%g  hhdx=%g, steps hh %i y ct %i\n", ct_dXtemp, hh_dXtemp, hh_steps_temp, ct_steps_temp);

  // Test logical options
  if(hh_dXtemp > node_size){
    _ERROR("Error. Your element size %g must be smaller than node of Ranvier size = %g\n",hh_dXtemp, node_size);
    exit(0);
  }
  if(node_size > inter_node){
    _ERROR("Error. Your node of Ranvier size = %g must be smaller than internodal distance =%g\n",node_size,inter_node);
    exit(0);
  }

  numStepsX = 0;
  fprintf(standard,"total elements steps=%i\n", numStepsX);

  int tempcts=0, temphhs=0;

  int size_branch = 50;// 
  int chunks=0;
  int temp_cont=0;
  int temp_terminals=0;


  int things=0;
  things = num_terminals*(size_branch/(hh_steps_temp+ct_steps_temp)+1);

  int flag_jump=0;
  int temp_path=0; //
  terminals = new int[num_terminals](); // () initialize to zero
  terminals[0] = 0;
  element =  new discretization*[size_branch*num_terminals+things];
  HH_elements = new Hodgkin_Huxley*[things];
  CT_elements = new Cable*[size_branch*num_terminals+things];

  while(temp_terminals < num_terminals){
    temp_cont = 0;
    while(temp_cont<size_branch){
      for(tempcts=0;tempcts<ct_steps_temp;tempcts++){
	temp_cont++;
	if(temp_terminals==2){
	  terminals[temp_terminals-1]=numStepsX-1;

	  CT_elements[num_CT_elements]=new Cable(numStepsX, numStepsX+1, flag_jump, ct_dXtemp);
	  element[numStepsX] = CT_elements[num_CT_elements];
	  num_CT_elements++;
	  element[numStepsX]->set_epsilon(in.epsilon);
	  element[flag_jump]->set_DL(numStepsX);
	  element[flag_jump]->set_branching(true);
	  element[numStepsX-1]->set_DR(numStepsX-1);
	    
	  temp_terminals++;

	}else{

	  CT_elements[num_CT_elements]=new Cable(numStepsX, numStepsX+1,numStepsX-1, ct_dXtemp);
	  element[numStepsX] = CT_elements[num_CT_elements]; 
	  num_CT_elements++;
	  element[numStepsX]->set_epsilon(in.epsilon);
	}
	numStepsX++;
      }
      for(temphhs=0;temphhs<hh_steps_temp;temphhs++){
	HH_elements[num_HH_elements] = new Hodgkin_Huxley(numStepsX, numStepsX+1,numStepsX-1, hh_dXtemp);
	element[numStepsX] = HH_elements[num_HH_elements]; 
	num_HH_elements++;
	element[numStepsX]->set_epsilon(in.epsilon);
	temp_cont++;
	numStepsX++;
      }
      chunks++;
    }
    if(temp_path==0){
      flag_jump=temp_cont-1;
    }
    temp_path++;
    temp_terminals++;
  }
   terminals[num_terminals-1] = numStepsX-1;
  element[numStepsX-1]->set_DR(numStepsX-1);
  dist_mes = numStepsX*ct_dXtemp/10;

 element[4]->set_input_current(in.input_current); // 
 element[90]->set_input_current(in.input_current); // 

}
/*! This is only for the HH model.
  @param in Structure with parameters defined in the command line
  @param element Array of discretization classes
  @param HH_elements This is an array of classes that contains the Hodgkin and Huxley elements in the discretization
  @param num_HH_elements Size of the HH_elements array
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array

  * Please, keep in your mind that you must define a new function following the way of the discretization. This function is only an example.
 */
extern void Only_HH_Model_class(input &in, discretization **&element,  Hodgkin_Huxley **&HH_elements,int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements){
  in.IRE = 1;
  in.NRE =1;
  in.dist_mes= 4.2E-3;
  xspan =6.3E-3;
  in.init_time = 0.003; // seconds. When you apply the current.
  in.input_time = 0.003;
  in.end_input = in.init_time + in.input_time;
  in.input_current = 0.8E-9;
  numStepsX =3;
  in.dt_factor = 0.0001;
  num_terminals=2;
  tspan = 0.025;
    
  element =  new discretization*[numStepsX];
  terminals = new int[num_terminals](); // () initialize to zero
  num_CT_elements = 2;
  CT_elements = new Cable*[num_CT_elements];
  num_HH_elements = 1;
  HH_elements = new Hodgkin_Huxley*[num_HH_elements];
  terminals[0] = 0;
  // First fictitious element
  CT_elements[0]=new Cable(0, 0+1, 0-1,2.1E-3);
  element[0] = CT_elements[0];
  element[0]->set_epsilon(0);
  // Active element
  HH_elements[0] = new Hodgkin_Huxley(1, 1,1-1,2.1E-3 );
  element[1] = HH_elements[0] ;
  element[1]->set_epsilon(0);
  // First fictitious element
  CT_elements[1]=new Cable(2, 1+1, 2-1,2.1E-3);
  element[2] = CT_elements[1];
  element[2]->set_epsilon(0);
  terminals[1] = 2;

 element[1]->set_input_current(in.input_current); // Only at the initial element

}

/*! General cable without branching with cable theory elements. Passive cable. This function does exactly the same than general_tree_n_terminal_CT with n=1 jumps
  @param in Structure with parameters defined in the command line
  @param element Array of discretization classes
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array

  * Please, keep in your mind that you must define a new function following the way of the discretization. This function is only an example.
 */
extern void without_branching_CT(discretization **&element,input &in, Cable **&CT_elements, int &num_CT_elements){
  TYPE dXtemp=0;
  int cable_steps_temp=0;
  in.input_time = 0.003; // To apply the external intensity
  in.input_current = 0.8E-9; // Nano Amperes. External Intensity
  in.end_input = in.init_time + in.input_time;   // when you finish to apply the current
  num_terminals =2; // The first and last elements

  if(xspan < 1E-10){
    xspan = 0.002; // Total length. It is approximate; xspan = xspan*(1+epsilon). We put two times the total distance in order not to have boundary problems
  } // if not is because you put as an argument.

  dist_mes = xspan/10;

  if(tspan < 1E-10){
    // Time discretization
    tspan = 0.025; // Total time 
  } // if not is because you put as an argument.

  if(numStepsX != 0){
    cable_steps_temp = numStepsX;
    dXtemp = (TYPE)xspan/cable_steps_temp;
    fprintf(stdout,"dX at the whole cable =%g\n",dXtemp);
  }else{
    cable_steps_temp = 110; // Total number of elements;
    dXtemp = (TYPE)xspan/cable_steps_temp;
    numStepsX = cable_steps_temp ;
    fprintf(stdout,"dX at the whole cable =%g\n",dXtemp);
  }

  element =  new discretization*[numStepsX];
  terminals = new int[num_terminals](); // () initialize to zero
  terminals[0] = 0;
  num_CT_elements = numStepsX;
  CT_elements = new Cable*[num_CT_elements];

  for(int i=0;i<numStepsX;i++){
 
    CT_elements[i]=new Cable(i, i+1,i-1, dXtemp);
    element[i] = CT_elements[i];	  
    element[i]->set_epsilon(in.epsilon);
  }
  terminals[1] = numStepsX-1;
 element[1]->set_input_current(in.input_current); // Only at the initial element
}

/*! General cable without branching with CT+HH elements. BMMB scheme (strain, refinement, etc)
  @param in Structure with parameters defined in the command line
  @param element Array of discretization classes
  @param HH_elements This is an array of classes that contains the Hodgkin and Huxley elements in the discretization
  @param num_HH_elements Size of the HH_elements array
  @param CT_elements This is an array of classes that contains the cable theory elements in the discretization
  @param num_CT_elements Size of the CT_elements array

  * Please, keep in your mind that you must define a new function following the way of the discretization. This function is only an example.
 */
extern void without_branching_CT_HH(input &in, discretization **&element, Hodgkin_Huxley **&HH_elements,int &num_HH_elements, Cable **&CT_elements, int &num_CT_elements){

  // Temporal values
  int i=0;
  double temp=0;
  int temp2=0;
  int flag = 0, cItr=0; // to know if we are at internodal or nodal section
  int *gPoints;
  int hh_numChannels =0, cable_numCable=0, cable_numStepsX=0, hh_numStepsX=0;
  double hh_epsilon=0, cable_epsilon=0, hh_dX=0, cable_dX, cable_inter_node = 0, hh_node_size = 0;

  num_terminals =2;
  terminals = new int[num_terminals](); // () initialize to zero
  terminals[0] = 0;

  in.init_time = 0.04;
  in.input_time = 0.003; // To apply the external intensity
  in.input_current = 0.04E-9; // Nano Amperes. External Intensity.
  //in.input_current = 0;
  hh_numChannels = 0;  // to initialize the variable
  cable_numCable =0; // to initialize the variable
  in.end_input = in.init_time + in.input_time;   // when you finish to apply the current
  dist_mes = 0.010; // The referential measurement distance. BMMB paper.


  if(in.IRE< 1E-10){
    // Default values for the discretization
    //cable_numStepsX = 20; // internodal part
    cable_numStepsX = 100; // internodal part
  }else{
    cable_numStepsX = in.IRE; // internodal part
  }
  if(in.NRE< 1E-10){
    // Default values for the discretization
    hh_numStepsX = 1; // nodal part
  }else{
    hh_numStepsX = in.NRE; // nodal part
  }
  if(xspan < 1E-10){
    
    xspan = 2;//0.02 // Total length. It is approximate; xspan = xspan*(1+epsilon). We put two times the total distance in order not to have boundary problems
  } // if not is because you put as an argument.

  if(tspan < 1E-10){
    // Time discretization
    tspan = 0.040; // 0.025 Total time
  } // if not is because you put as an argument.

  if(hh_node_size < 1E-10){
    //  hh_node_size = 2.1E-6; // Node of Ranvier. From Shi papers
    hh_node_size = 2.1E-6; // Node of Ranvier. From Shi papers
  }
  if(cable_inter_node < 1E-10){
    //cable_inter_node = 0.0008; // internodal length. General distance
   cable_inter_node = 0.0008; // internodal length. General distance
  }
  //  Total number of elements
  temp2 = int((xspan/( cable_inter_node + hh_node_size))+1); // Number of CT+HH parts
  xspan = temp2*( cable_inter_node + hh_node_size); // Recalculate the total length
 
  hh_epsilon = (in.nu*xspan/(hh_node_size*temp2))*in.epsilon;
  cable_epsilon = ((1-in.nu)*xspan/(temp2*cable_inter_node))*in.epsilon;

  if(numStepsX < 1E-10){
    numStepsX = temp2*(cable_numStepsX+hh_numStepsX); // Calculate the total spatial steps
  }

  //Calculate dXi, dXn and dT
  cable_dX = (double)cable_inter_node/cable_numStepsX;
  hh_dX = (double)hh_node_size/hh_numStepsX;

  fprintf(standard,"for discretization and simulation cabledx =%g  hhdx=%g steps cable %i steps hh=%i\n", cable_dX, hh_dX, cable_numStepsX, hh_numStepsX);
  fprintf(standard,"distances cable=%.20f  hh=%.20f\n", cable_dX*cable_numStepsX, hh_dX*hh_numStepsX);

  gPoints = new int[temp2*hh_numStepsX*2](); // this 2 is to get enough memory. It does not matter, because this vector will be deleted before the simulation

  // Test logical options
  if(hh_dX > hh_node_size){
    _ERROR("Error. Your element size %g must be smaller than node of Ranvier size = %g\n",hh_dX, hh_node_size);
    exit(0);
  }
  if(hh_node_size > cable_inter_node){
    _ERROR("Error. Your node of Ranvier size = %g must be smaller than internodal distance =%g\n",hh_node_size, cable_inter_node);
    exit(0);
  }

  // be careful with this part, it is complicated
  i=1; // First ion channel is always at the first element. The element 0 is a fictitious element (BC).
  num_CT_elements++;

  temp=0;
  double tolhh=hh_dX/2;
  double tolct=cable_dX/2;
  while(i<numStepsX){
    if(flag == 0){
      gPoints[cItr] = i;
      cItr++;
      hh_numChannels++;
      temp = temp + hh_dX;
      i++;
      if(fabs(temp + tolhh) > hh_node_size){
	flag=1; // change to internodal part
	temp=0;
      }
    }
    if(flag == 1){
      i++;
      num_CT_elements++;
      temp = temp + cable_dX;
      if(fabs(temp + tolct) > cable_inter_node){
	flag=0; // change to nodal part
	temp=0;
	cable_numCable++;
      }
    } 
  }
  fprintf(standard, "Total numChannels %i. With nodal size=%f and %f as internodal distance.\n",hh_numChannels,hh_node_size, cable_inter_node);

  element =  new discretization*[numStepsX];
  num_HH_elements = hh_numChannels;
  HH_elements = new Hodgkin_Huxley*[num_HH_elements];
  CT_elements = new Cable*[num_CT_elements];
  fprintf(standard, "Total numChannels %i cable =%i total = %i numStepsX=%i\n",hh_numChannels, num_CT_elements, num_HH_elements + num_CT_elements, numStepsX);
  int tempo=0, tempoct=0;
 
  for(i=0;i<numStepsX;i++){
    if(((int*) bsearch(&i, gPoints, hh_numChannels, sizeof(int), compareints2)) != NULL){
      HH_elements[tempo] = new Hodgkin_Huxley(i, i+1,i-1, hh_dX);
      element[i] = HH_elements[tempo];
      element[i]->set_epsilon(hh_epsilon);
      tempo++; 
    }else{
      CT_elements[tempoct]=new Cable(i, i+1,i-1, cable_dX);
      element[i] = CT_elements[tempoct];
      element[i]->set_epsilon(cable_epsilon);  
      tempoct++;
    }
  }
  element[numStepsX-1]->set_DR(numStepsX-1);
  fprintf(standard,"nu=%g\n",in.nu);
  cable_numCable++;
  terminals[1] = numStepsX-1;
  fprintf(standard,"nu of equilibrium:%g\n",hh_node_size*hh_numChannels/xspan);

  fprintf(standard,"Setting the input current at the elements that you defined %g\n",in.input_current);
  element[1]->set_input_current(in.input_current); // Only at the initial element

  delete gPoints;
}



