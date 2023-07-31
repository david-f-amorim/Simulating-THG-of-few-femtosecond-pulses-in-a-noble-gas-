DESY Summer Programme 2023
============================

GOAL:     simulate THG of a few-femtosecond IR laser beam in a 
----      noble gas; 
	  carry out parameter scans to find optimal input beam 
          conditions for desired UV output beam;   
		 


Stage 1: DONE [Week 1-2]
--------------------------
	
did background reading ; cleaned up existing code ; 
provided visualisation functions ; provided output handling 

Stage 2: IN PROGRESS [Week 2]
----------------------------------

provided option to base input pulse on measured data ; provided 
option to feed in custom gas density profile (untested!) ; 
provided option to overlay measured UV spectrum 
	

Stage 3: TO DO
-----------------------------

	- do pressure scans (Ar & Ne) and try to reproduce signs of filamentation (saturation, pulse shape change, second spectral peak)

			+ consider input pulses at 75mW, 150mw, and 300mW (currently only 150mW data available)
			+ so far: filamentation effects could not be reproduced 
				=> most likely due to gas density model breaking down at high pressures?

				CAREFULLY EVALUATE PRESSURE SCANS AT LOW PRESSURES TO EVALUATE THIS HYPOTHESIS!!

				look into using and improving COMSOL simulations 			


	- make simulation more realistic: 

		+ set up non-Gaussian spatial profile 
		+ input measured CEO phase value ?

