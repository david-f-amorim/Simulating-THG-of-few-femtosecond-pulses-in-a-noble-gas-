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

	- set up "Fields.DataField" to represent input pulse (replacing "Fields.GaussGaussField")
		
		+ extract spectral phase from FROG data  (read through Jupyter Notebook)
		+ combine spectral phase with frequency and spectral power density
                  into a data file 

	- compare input and output with measured data 

		+ NIR spectrum at input compared to measured data
		+ UV spectrum at output compared to measured data
		+ ....


Stage 3: TO DO
-----------------------------

	- vary pressure, beam intensity, etc. to study effect on THG efficiency [Luna parameter scans...]
		-> try to reproduce Josina's measured pressure scans (saturation, pulse shape change, second spectral peak)

	- set up own density function