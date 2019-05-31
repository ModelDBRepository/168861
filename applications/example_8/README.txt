In this simulation the segmented neuron is loaded. The initial segmented neuron must be modified cutting a branch departing from the soma and adapting the enumeration to Neurite's enumeration. You have an example with a Matlab code in the Others folder of the program. The initial neuron for this example is the ID: NMO_00223 from neuromorpho.org, Ishizuka, Cowan and Amaral, 1995. A quantitative analysis of the dendritic organization of pyramidal cells in the rat hippocampus. The Journal of Comparative Neurology.

neuron_for_neurite.txt->    This file is a neuron and is loaded by Neurite to build the discretization. The information is by columns as follows:
			     
			     element	parent	daughter_right	daughter_left	diameter(micro)	length (micrometers)
			     	
In the output folder, in the *.dat file,you have the potential vs time at a given point along the tree.
