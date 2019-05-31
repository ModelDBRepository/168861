#
#
# File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
#
# Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
#
#
#!/bin/bash
if [ $# -ne 3 ]; then
    cat >&2<<EOF

¡¡¡ERROR!!!: The number of parameters is not correct ($#).
==================================================
 $0 'solver' 'type' 'configuration'
----------------------------------
 Example:
   $0 cpu I case_defined_by_the_user
==================================================
EOF
    exit -1
fi

export idGPU=0;  # In case you have more than one graphics card in your machine, you must give the corresponding label to Neurite. Typing nvidia-smi in the terminal you can see the characteristics of your graphics cards

DIRECTORY=$(pwd);
export working_directory=$(pwd);
rm -rf "outputs"
mkdir "outputs";
mkdir "outputs/GUI_captures";

# Type of simulation to run (the functions to create the discretization are in the space.cpp file):
   # simulation=
                 # Only_HH_Model: only one HH element. This simulation is only available with the sequential version of the explicit scheme
                 # HH_axon: axon without branching 
                 # HH_Y_branching: axon with one branching point
                 # CT_cable:  passive structure without branching
                 # CT_Y_branching:  passive structure with one branching point
                 # CT_symmetric_tree: passive dendritic symmetric structure with several branching points
                 # CT_segmented_neuron: segmented neuron structure

# Type of execution in the simulation environment
export simulation=CT_cable; 

# Define the processor to use
# Warning: You have to put the number of "../" to reach the bin folder!!!
case $1 in 
    "cpu")
	COMM=../../bin/Neurite_solver-cpu
	SOLVER="solver-cpu"
	;;
    "gpu")
	COMM=../../bin/Neurite_solver-gpu
	SOLVER="solver-gpu"
	;;
    *)
	echo "ERROR: this processor, ($1), does not exist [cpu or gpu]">&2
	exit -1;
	;;
esac

# The solver: Explicit or Implicit
case $2 in
    "E")
	S_TYPE=e
	dtf=0.8
	;;
    "I")
	S_TYPE=i
	dtf=100
	;;
    *)
	echo "ERROR: this solver ($2) does not exist [E or I]">&2
	exit -1;
	;;
esac

# General parameters to run Neurite. Type man ./neurite.7 in the documentation folder to see the list of parameters
case $3 in
    "example_4")
	# dtf implicit=10 explicit=0.8
	echo "Explicit simulations for one axon"
       	OPT="g Tt 0.1 dtf ${dtf}"
	PATH="$DIRECTORY/outputs"
	;;
    *)
	echo "ERROR: this type of simulation, ($3), does not exist. Please, you should create and configure your own simulation in this file">&2
	exit -1;
	;;
esac
if [ ! -d ${PATH} ]; then
    echo "ERROR: you did not create the folder to place the results of the simulations">&2
    exit -1;
fi
# The command line that you are executing
echo "time ${COMM} ${S_TYPE} ${OPT} out ${PATH}/${SOLVER}_${S_TYPE}_simulation.dat log ${PATH}/${SOLVER}_${S_TYPE}_simulation.log"
time ./${COMM} ${S_TYPE} ${OPT} out ${PATH}/${SOLVER}_${S_TYPE}_simulation.dat log ${PATH}/${SOLVER}_${S_TYPE}_simulation.log 
