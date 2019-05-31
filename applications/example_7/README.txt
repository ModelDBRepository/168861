This is the example 7, an axon under mechanical loading defined in the input_data.txt file. In the output folder, in the *.dat file, you have the potential vs time at a given point along the axon

Execution line:

./example_7.sh processor solver example_7

with processor = cpu or gpu, and solver = E or I.

input_data.txt->     	     This file is mandatory when you are using the mechanical model. It is loaded by the mechinical_model.cpp file to extract the strain configuration.
	      	 	     The file is composed by only one line with the following values (as columns):
		 	     1-  Case    	      	      	       Integer that represents the specific case in the Shi's paper experiments. 1, 2, 3, 4, 5 or 6.
		 	     2-  Kappa   		      	       See the mechanical model paper (BMMB)
		 	     3-  Epsilon macro at t=0     	       Calculated by the mechanical 
		 	     4-  Epsilon damage at t=0                 Calculated by the mechanical 
		 	     5-  Epsilon damage at t=infinity	       Calculated by the mechanical

			     If you are not using the mechanical model, this file is not necessary.


