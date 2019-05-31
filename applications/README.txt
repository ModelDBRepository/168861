It is strongly recommended to define the simulation in an execution.sh bash file. It is also possible executing directly the executable file defining previously the folders, variables, etc. Following the template given in this folder will simplify dramatically the execution of a simulation. 

In each of the examples exposed here, you have the corresponding bash file to configure and execute the simulation.

Once you define the folder to work and configure the bash file, you can execute the simulation as follows:

./example.sh processor solver case_defined_by_the_user

where the parameters for the bash file are

processor:

cpu    -> This is the sequential version of Neurite
gpu    -> This is the parallel version of Neurite

solver:

I -> For the implicit scheme. Use a dtf bigger than 1
E -> For the explicit scheme. Use a dtf smaller or equal to one

case_defined_by_the_user -> You define inside the bash file a case in which you configure the parameters and options needed by Neurite. Type man ./neurite.7 inside the documentation folder for further information.

This is the list of examples shown here:
example_1 -> Only the HH model. 1 element
example_2 -> Axon myelinated without branching
example_3 -> Axon one branching point
example_4 -> Passive cable without branching
example_5 -> Passive cable with one branching point
example_6 -> Passive symmetric tree synthetically created
example_7 -> Axon under mechanical loading
example_8 -> Passive segmented dendrites
example_9 -> Long myelinated axon: parallel advantages
