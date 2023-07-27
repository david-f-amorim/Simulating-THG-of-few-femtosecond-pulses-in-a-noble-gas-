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

Stage 2: IN PROGRESS [Week 2 - ?]
----------------------------------

	- use measured data as input (replacing "Fields.GaussGaussField") [DONE!]
		
	- compare input and output with measured data 

		+ NIR spectrum at input compared to measured data
			-> from Speck.dat  
		+ measured NIR pulse at input compared to simulated NIR pulse at input 
			-> from Ek.dat 

		+ UV spectrum at output compared to measured data
			-> from data folder [spectra]
		+ UV pulse energy compared to measured data
			-> from data folder [pressure scans]


Stage 3: TO DO
-----------------------------

	- make simulation more realistic: 

		+ set up own density function
		+ set up non-Gaussian spatial profile 
		+ input measured CEO phase value 

	- vary pressure, beam intensity, etc. to study effect on THG efficiency [Luna parameter scans...] (Ar & Ne)
		-> try to reproduce Josina's measured pressure scans (saturation, pulse shape change, second spectral peak)

