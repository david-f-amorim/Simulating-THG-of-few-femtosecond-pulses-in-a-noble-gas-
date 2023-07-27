DESY Summer Programme 2023
============================

GOAL:     simulate THG of a few-femtosecond IR laser beam in a 
----      noble gas; 
	  carry out parameter scans to find optimal input beam 
          conditions for desired UV output beam;   
		 


Stage 1: DONE [Week 1-2 ]
--------------------------
	
did background reading ; cleaned up existing code ; 
provided visualisation functions ; provided output handling 

	(MAYBE: reproduce Josina's results ?) 


Stage 2: IN PROGRESS [Week 2 - ?]
----------------------------------

	- use measured data as input (replacing "Fields.GaussGaussField") [DONE!]

		-> can the CEO phase be extracted from the data?
		-> is pulse duration still acurate? 


		-> use measured spatial beam data to improve input (if exists)
		
	- compare input and output with measured data 

		+ NIR spectrum at input compared to measured data
			-> extract from FROG Speck.dat  ?
		+ measured NIR pulse at input compared to simulated NIR pulse at input 
			-> from Ek.dat 

		+ UV spectrum at output compared to measured data
			-> 
		+ ....


Stage 3: TO DO
-----------------------------

	- vary pressure, beam intensity, etc. to study effect on THG efficiency [Luna parameter scans...]
		-> try to reproduce Josina's measured pressure scans (saturation, pulse shape change, second spectral peak)

	- set up own density function ? set up non-Gaussian spatial profile ?